import pandas as pd
import neutcurve
import altair as alt
import sys
# allow more rows for Altair
_ = alt.data_transformers.disable_max_rows()

#import altair theme and enable
sys.path.append(str(snakemake.input.theme))
import theme
alt.themes.register('main_theme', theme.main_theme)
alt.themes.enable('main_theme')

# Path to the files for snakemake
neut_file_path = str(snakemake.input.neutFile)
fitParams_output = str(snakemake.output.fitParams)
output_image_path = str(snakemake.output.neutcurve_img)

# Snakemake parameter flags
icvalues = snakemake.params.icvalues
height = snakemake.params.height
width = snakemake.params.width

# Flags to determine the type of curve to plot
receptor_flag = True
sera_flag = False
antibody_flag = False

df = pd.read_csv(neut_file_path)

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
fitParams.to_csv(fitParams_output,index=False)

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
        print('Please specify the type of curve to plot')
    
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
ephrin_curve.save(output_image_path,ppi=300)
