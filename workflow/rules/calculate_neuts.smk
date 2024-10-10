rule neuts:
    input:
        neutFile="data/{sample}.csv",
        theme="config/",  #altair theme
    output:
        fitParams="results/fitparams/{sample}.csv",
        neutcurve_img="results/images/{sample}.png",
        neutcurve_svg="results/images/{sample}.svg",
    params:
        icvalues=config["icvalues"],
        height=config["height"],
        width=config["width"],
    conda:
        "../envs/environment.yml"
    log:
        "results/logs/{sample}_neuts.txt",
    script:
        "../scripts/neuts.py"
