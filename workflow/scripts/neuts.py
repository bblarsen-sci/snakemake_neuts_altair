import pandas as pd
import neutcurve
import altair as alt
import re
import os
import httpimport

_ = alt.data_transformers.disable_max_rows()

# Import custom altair theme from github using httpimport module
def import_theme(github_username, github_repo, github_branch):
    with httpimport.github_repo(github_username, github_repo, github_branch):
        import main_theme
    alt.themes.register("main_theme", main_theme.main_theme)
    alt.themes.enable("main_theme")


# Read in the neutralization data
def load_data(file_path):
    df = pd.read_csv(file_path)

    REQUIRED_COLUMNS = [
        "serum",
        "virus",
        "replicate",
        "concentration",
        "fraction infectivity",
    ]
    missing_columns = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_columns:
        log_file.write(f"Missing required columns: {', '.join(missing_columns)}\n")
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")
    return df

# Determine which type of neut was run by parsing the name of the input file
def process_filename(file_path):
    input_file = str(file_path)
    filename = os.path.basename(input_file)
    filename_lower = filename.lower()

    result = {
        'type': None,
        'facet': False
    }

    if "receptor" in filename_lower:
        result['type'] = "receptor"
    elif "sera" in filename_lower:
        result['type'] = "sera"
    elif "antibody" in filename_lower:
        result['type'] = "antibody"
    else:
        log_file.write(f"Unknown sample type in filename: {filename}\n")
        raise ValueError(f"Unknown sample type in filename: {filename}")

    if "facet" in filename_lower:
        result['facet'] = True

    return result


# Determine if serum variables or virus variables should be used in the plot
def determine_experiment_type(dataframe):
    number_serum = len(dataframe["serum"].unique())
    number_virus = len(dataframe["virus"].unique())
    log_file.write(f"Number of serum variables: {number_serum}\n")
    log_file.write(f"Number of virus variables: {number_virus}\n")
    return number_serum > number_virus, number_serum <= number_virus


# Estimate neutralization curves using the `curvefits` module from `neutcurve` package.
def fit_neutcurve(df):
    # estimate fits
    return neutcurve.curvefits.CurveFits(
        data=df,
        serum_col="serum",
        virus_col="virus",
        replicate_col="replicate",
        conc_col="concentration",
        fracinf_col="fraction infectivity",
    )

def get_ic_values(fits):
    fit_params = fits.fitParams(ics=ICVALUES)
    fit_params.to_csv(snakemake.output.fitParams, index=False)

def make_neutcurve_df(fits):
    curves = []  # initialize an empty list to store neutralization curve data
    # Loop over each serum type and retrieve the curve
    for serum in serum_list:
        for virus in virus_list:
            curve = fits.getCurve(serum=serum, virus=virus, replicate="average")
            neut_df = curve.dataframe()  # turn into a dataframe
            neut_df["serum"] = serum  # assign serum name to a column
            neut_df["virus"] = virus  # assign virus name to a column
            curves.append(neut_df)

    # Concatenate all the dataframes into one
    combined_curve = pd.concat(curves, axis=0)
    combined_curve["upper"] = combined_curve["measurement"] + combined_curve["stderr"]
    combined_curve["lower"] = combined_curve["measurement"] - combined_curve["stderr"]
    return combined_curve

def custom_sort_order(array):
    # Helper function to extract numerical part from mutation strings.
    def extract_number(virus):
        num = re.search(r"\d+", virus)
        return (
            int(num.group()) if num else 0
        )  # Convert digits to integer, or 0 if none found.

    array = sorted(
        array, key=extract_number
    )  # Sort array by the numerical value extracted.

    # Ensure 'WT' (wild type) is the first element in the list if it exists.
    if "WT" in array:
        array.remove("WT")  # Remove 'WT' from its current position.
        array.insert(0, "WT")  # Insert 'WT' at the beginning of the list.
    return array

def custom_sort_key(item):
    # Extract the time value from the item string
    match = re.search(r'(\d+)min', item)
    if match:
        return int(match.group(1))
    return 0  # Return 0 for items without a time value


# Plot the neutralization curves
def plot_neut_curve(df):
    if sample_type == "receptor":
        title = "Concentration (µM)"
        if vary_serum_flag:
            legend_title = "Receptor"
            color_variable = "serum"
        elif vary_virus_flag:
            legend_title = "Virus"
            color_variable = "virus"
    elif sample_type == "sera":
        title = "Sera Dilution"
        if vary_serum_flag:
            legend_title = "Serum"
            color_variable = "serum"
        elif vary_virus_flag:
            legend_title = "Virus"
            color_variable = "virus"
    elif sample_type == "antibody":
        title = "Concentration (µg/mL)"
        if vary_serum_flag:
            legend_title = "Antibody"
            color_variable = "serum"
        elif vary_virus_flag:
            legend_title = "Virus"
            color_variable = "virus"
    else:
        log_file.write(f"Unknown sample type: {sample_type}\n")
        raise ValueError("Unknown sample type")

    axis = alt.Axis(format=".0e", tickCount=3)
    scale = alt.Scale(type="log")



    # Get the unique values of the color variable
    unique_values = df[color_variable].unique()
    
    # Sort the unique values using the custom sorting function
    sorted_values = sorted(unique_values, key=custom_sort_key)
    
    # If 'WT' is present, set the color scale to include 'WT' as the first color and make it black
    if "WT" in df["virus"].unique():
        print("WT is present")
        colors = ["black"] + COLORS_FOR_PLOT[: len(df["virus"].unique()) - 1]
        color_scale = alt.Color(
            color_variable,
            title=legend_title,
            scale=alt.Scale(
                domain=custom_sort_order(df["virus"].unique()), range=colors
            ),
        )
    else:
        color_scale = alt.Color(
            color_variable,
            title=legend_title,
            scale=alt.Scale(domain=sorted_values,range=COLORS_FOR_PLOT),
            sort=sorted_values
        )
        
    chart = (
        alt.Chart(df)
        .mark_line(size=LINE_WIDTH)
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
                axis=alt.Axis(values=[0, 0.5, 1]),
            ),
            color=color_scale,
        )
    )
    circle = (
        alt.Chart(df)
        .mark_circle(size=CIRCLE_SIZE, opacity=1)
        .encode(
            x=alt.X(
                "concentration",
                scale=scale,
                axis=axis,
                title=title,
            ),
            y=alt.Y("measurement:Q", title="Fraction Infectivity"),
            color=color_scale,
        )
    )
    error = (
        alt.Chart(df)
        .mark_errorbar(opacity=ERROR_BAR_OPACITY)
        .encode(
            x="concentration",
            y=alt.Y("lower", title="Fraction Infectivity"),
            y2="upper",
            color=color_scale,
        )
    )
    plot = chart + circle + error
    plot = plot.properties(width=WIDTH, height=HEIGHT)
    
    if use_facet:
        facet_variable = 'serum' if vary_virus_flag else 'virus'
        plot = plot.facet(
            facet=alt.Facet(
                f"{facet_variable}:N",
                header=alt.Header(title=None),
            ),
            columns=N_FACETED_COLUMNS,
        )

    return plot




# Main execution
if __name__ == "__main__":
    try:
        log_file = open(str(snakemake.log), "w")
        log_file.write(f"Input file: {snakemake.input.neutFile}\n")

        ICVALUES = snakemake.params.icvalues
        COLORS_FOR_PLOT = snakemake.params.colors
        HEIGHT = snakemake.params.height
        WIDTH = snakemake.params.width
        N_FACETED_COLUMNS = snakemake.params.n_faceted_columns
        LINE_WIDTH = snakemake.params.line_width
        CIRCLE_SIZE = snakemake.params.circle_size
        ERROR_BAR_OPACITY = snakemake.params.error_bar_opacity

        github_username = snakemake.params.github_username
        github_repo = snakemake.params.github_repo
        github_branch = snakemake.params.github_branch

        # import custom altair theme from github
        import_theme(github_username, github_repo, github_branch)

        # Load data
        df = load_data(snakemake.input.neutFile)

        # Determine sample type and faceting
        filename_info = process_filename(snakemake.input.neutFile)
        sample_type = filename_info['type']
        use_facet = filename_info['facet']

        # Determine experiment type
        vary_serum_flag, vary_virus_flag = determine_experiment_type(df)

        # Get neutralization curves and IC values
        fit = fit_neutcurve(df)
        get_ic_values(fit)

        # get list of different sera and viruses that were tested
        serum_list = list(df["serum"].unique())
        virus_list = list(df["virus"].unique())

        neutcurve_df = make_neutcurve_df(fit)

        # Plot the neutralization curves
        neut_curve = plot_neut_curve(neutcurve_df)
        neut_curve.save(snakemake.output.neutcurve_img, ppi=300)
        neut_curve.save(snakemake.output.neutcurve_svg)

    except Exception as e:
        log_file.write(f"An error occurred: {str(e)}")
    log_file.close()
