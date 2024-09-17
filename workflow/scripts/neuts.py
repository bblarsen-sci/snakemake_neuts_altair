import pandas as pd
import neutcurve
import altair as alt
import sys
import re
_ = alt.data_transformers.disable_max_rows()
#import altair theme and enable
sys.path.append(str(snakemake.input.theme))
import theme
alt.themes.register('main_theme', theme.main_theme)
alt.themes.enable('main_theme')

# Path to the files for snakemake
#neut_file_path = str(snakemake.input.neutFile)
fitParams_output = str(snakemake.output.fitParams)
output_image_path = str(snakemake.output.neutcurve_img)
output_log_file = str(snakemake.log)

log_file = open(output_log_file, 'w')
log_file.write(f"Input file: {snakemake.input.neutFile}\n")

# Snakemake parameter flags
icvalues = snakemake.params.icvalues
height = snakemake.params.height
width = snakemake.params.width
sample_type = snakemake.params.sample_type

# Set flags based on sample type
if sample_type == 'receptor':
    receptor_flag = True
    sera_flag = False
    antibody_flag = False
elif sample_type == 'sera':
    receptor_flag = False
    sera_flag = True
    antibody_flag = False
elif sample_type == 'antibody':
    receptor_flag = False
    sera_flag = False
    antibody_flag = True
else:
    raise ValueError(f"Unknown sample type: {sample_type}")

# Read in the neutralization data

def load_data(file_path):
    try:
        df = pd.read_csv(file_path)
    
        REQUIRED_COLUMNS = ['serum', 'virus', 'replicate', 'concentration', 'fraction infectivity']
        missing_columns = [col for col in REQUIRED_COLUMNS if col not in df.columns]
        if missing_columns:
            log_file.write(f"Missing required columns: {', '.join(missing_columns)}\n")
            raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")
        return df
    except Exception as e:
        raise ValueError(f"Error loading data from {file_path}: {e}")

df = load_data(snakemake.input.neutFile)


def determine_experiment_type(dataframe):
    number_serum = len(dataframe["serum"].unique())
    number_virus = len(dataframe["virus"].unique())
    if number_serum > number_virus:
        return "serum"
    else:
        return "virus"
    
experiment_type = determine_experiment_type(df)
print(experiment_type)



# Determine if the serum or virus is being varied in the experiment
number_serum = len(df["serum"].unique())
number_virus = len(df["virus"].unique())
if number_serum > number_virus:
    vary_serum_flag = True
    vary_virus_flag = False
else:
    vary_serum_flag = False
    vary_virus_flag = True

# Estimate neutralization curves using the `curvefits` module from `neutcurve` package.
def get_neutcurve(df, replicate="average"):
    #estimate fits
    fits = neutcurve.curvefits.CurveFits(
        data=df,
        serum_col="serum",
        virus_col="virus",
        replicate_col="replicate",
        conc_col="concentration",
        fracinf_col="fraction infectivity",
        fixbottom=0,
    )
    
    fitParams = fits.fitParams(ics=icvalues)

    #get list of different sera and viruses that were tested
    serum_list = list(df["serum"].unique())
    virus_list = list(df["virus"].unique())

    curves = [] #initialize an empty list to store neutralization curve data
    
    # Loop over each serum type and retrieve the curve
    for serum in serum_list:
        for virus in virus_list:
            curve = fits.getCurve(serum=serum, virus=virus, replicate=replicate)
            neut_df = curve.dataframe() #turn into a dataframe
            neut_df["serum"] = serum #assign serum name to a column
            neut_df["virus"] = virus #assign virus name to a column
            curves.append(neut_df)

    # Concatenate all the dataframes into one
    combined_curve = pd.concat(curves, axis=0)
    combined_curve["upper"] = combined_curve["measurement"] + combined_curve["stderr"]
    combined_curve["lower"] = combined_curve["measurement"] - combined_curve["stderr"]
    
    return combined_curve, fitParams


neutcurve_df,fitParams = get_neutcurve(df)

################### Export modified fitParams to csv ###################
#fitParams = fitParams.drop(['replicate','nreplicates'],axis=1)

# Function to rename columns
def rename_ic_columns(col_name):
    match = re.match(r'(ic\d{2}$)', col_name)
    if match:
        return f"{match.group(1)}_ug"
    return col_name

# Rename columns
fitParams.rename(columns=rename_ic_columns, inplace=True)

# Function to create new 'ng' columns
def create_ng_columns(df):
    for col in df.columns:
        match = re.match(r'(ic\d{2})_ug', col)
        if match:
            new_col_name = f"{match.group(1)}_ng"
            df[new_col_name] = df[col].mul(1000).round(1)
    return df

fitParams = create_ng_columns(fitParams)

def filter_columns(df):
    return df.loc[:, df.columns.str.contains('serum|virus|_ng', case=False)]

# Apply the column filter
df = filter_columns(fitParams)
df.to_csv(fitParams_output,index=False)

# Plot the neutralization curves
def plot_neut_curve(df):
    if receptor_flag:
        scale = alt.Scale(type='log')
        axis = alt.Axis(format='.0e',tickCount=3)
        title = 'Concentration (µM)'
        legend_title = 'Receptor'
    elif sera_flag:
        scale = alt.Scale(type='log')
        axis = alt.Axis(format='.0e',tickCount=3)
        title = 'Sera Dilution'
        legend_title = 'Serum'
    elif antibody_flag:
        scale = alt.Scale(type='log')
        axis = alt.Axis(format='.0e',tickCount=3)
        title = 'Concentration (µg/mL)'
        legend_title = 'Antibody'
    else:
        raise ValueError(f"Unknown sample type")
    
    if vary_serum_flag:
        color_variable = 'serum'
    elif vary_virus_flag:
        color_variable = 'virus'
        
    
    chart = (
        alt.Chart(df)
        .mark_line(size=1.5)
        .encode(
            x=alt.X(
                "concentration:Q",
                scale=scale,
                axis=axis,
                title=title,
            ),
            y=alt.Y(
                "fit:Q",
                title="Fraction Infectivity",
            ),
            color=alt.Color(color_variable, title=legend_title),
        )
    )
    circle = (
        alt.Chart(df)
        .mark_circle(size=40,opacity=1)
        .encode(
            x=alt.X(
                "concentration",
                scale=scale,
                axis=axis,
                title=title,
            ),
            y=alt.Y("measurement:Q", title="Fraction Infectivity"),
            color=alt.Color(color_variable, title=legend_title),
        )
    )
    error = (
        alt.Chart(df)
        .mark_errorbar(opacity=1)
        .encode(
            x="concentration",
            y=alt.Y("lower", title="Fraction Infectivity"),
            y2="upper",
            color=color_variable,
        )
    )
    plot = chart + circle + error
    plot = plot.properties(width=width,height=height)
    return plot


ephrin_curve = plot_neut_curve(neutcurve_df)

#save plots
ephrin_curve.save(output_image_path, ppi=300)
log_file.close()