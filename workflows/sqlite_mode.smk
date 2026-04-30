rule all:
    input:
        "results/repro/README.txt"


rule scaffold_marker:
    output:
        "results/repro/README.txt"
    shell:
        "mkdir -p results/repro && printf 'GasRegNet scaffold workflow placeholder\\n' > {output}"
