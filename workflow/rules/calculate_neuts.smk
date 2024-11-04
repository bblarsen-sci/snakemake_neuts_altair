rule neuts:
    input:
        neutFile="data/{sample}.csv",
    output:
        fitParams="results/fitparams/{sample}.csv",
        neutcurve_img="results/images/{sample}.png",
        neutcurve_svg="results/images/{sample}.svg",
    params:
        icvalues=config["icvalues"],
        colors=config["colors"],
        height=config["height"],
        width=config["width"],
        github_username=config["github_username"],
        github_repo=config["github_repo"],
        github_branch=config["github_branch"],
        n_faceted_columns=config["n_faceted_columns"],
        line_width=config["line_width"],
        circle_size=config["circle_size"],
        error_bar_opacity=config["error_bar_opacity"],

    conda:
        "../envs/environment.yml"
    log:
        "results/logs/{sample}_neuts.txt",
        notebook = "results/logs/{sample}_neuts.py.ipynb"
    notebook:
        "../notebooks/neuts.py.ipynb"
