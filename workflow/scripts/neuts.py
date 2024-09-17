import pandas as pd
import neutcurve
import altair as alt
import sys
import re
import os
_ = alt.data_transformers.disable_max_rows()

#import altair theme and enable
sys.path.append(str(snakemake.input.theme))
import theme
alt.themes.register('main_theme', theme.main_theme)
alt.themes.enable('main_theme')

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

def process_filename(file_path):
    input_file = str(file_path)
    filename = os.path.basename(input_file)
    filename_lower = filename.lower()

    if 'receptor' in filename_lower:
        return 'receptor'
    elif 'sera' in filename_lower:
        return 'sera'
    elif 'antibody' in filename_lower:
        return 'antibody'
    else:
        raise ValueError(f"Unknown sample type in filename: {filename}")
    
# Determine if serum variables or virus variables should be used in the plot
def determine_experiment_type(dataframe):
    number_serum = len(dataframe["serum"].unique())
    number_virus = len(dataframe["virus"].unique())
    return number_serum > number_virus, number_serum <= number_virus
    
# Estimate neutralization curves using the `curvefits` module from `neutcurve` package.
def get_neutcurve(df, icvalues, replicate="average"):
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

# Plot the neutralization curves
def plot_neut_curve(df, vary_serum_flag, vary_virus_flag, sample_type):
    if sample_type == 'receptor':
        scale = alt.Scale(type='log')
        axis = alt.Axis(format='.0e',tickCount=3)
        title = 'Concentration (µM)'
        legend_title = 'Receptor'
    elif sample_type == 'sera':
        scale = alt.Scale(type='log')
        axis = alt.Axis(format='.0e',tickCount=3)
        title = 'Sera Dilution'
        legend_title = 'Serum'
    elif sample_type == 'antibody':
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
                scale=alt.Scale(domain=[0, 1]),
                axis=alt.Axis(values=[0,0.5,1]),
            ),
            color=alt.Color(color_variable, title=legend_title),
        )
    )
    circle = (
        alt.Chart(df)
        .mark_circle(size=30,opacity=1)
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
    plot = plot.properties(width=snakemake.params.width,height=snakemake.params.height)
    return plot

# Main execution
if __name__ == "__main__":
    try:
        log_file = open(str(snakemake.log), 'w')
        log_file.write(f"Input file: {snakemake.input.neutFile}\n")

        # Load data
        df = load_data(snakemake.input.neutFile)
        
        # Determine sample type
        sample_type = process_filename(snakemake.input.neutFile)
        
        # Determine experiment type
        vary_serum_flag, vary_virus_flag = determine_experiment_type(df)

        # Get neutralization curves and IC values
        neutcurve_df, fit_params = get_neutcurve(df, snakemake.params.icvalues)
        
        # Save the fit parameters
        fit_params.to_csv(snakemake.output.fitParams, index=False)
        
        # Plot the neutralization curves
        neut_curve = plot_neut_curve(neutcurve_df, vary_serum_flag, vary_virus_flag, sample_type)
        neut_curve.save(snakemake.output.neutcurve_img, ppi=300)
        neut_curve.save(snakemake.output.neutcurve_svg)
        
    except Exception as e:
        log_file.write(f"An error occurred: {str(e)}")
    log_file.close()