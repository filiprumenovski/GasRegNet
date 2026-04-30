rule all:
    input:
        "results/full_discovery/README.txt"


rule scaffold_marker:
    output:
        "results/full_discovery/README.txt"
    shell:
        "mkdir -p results/full_discovery && printf 'GasRegNet full-discovery workflow placeholder\\n' > {output}"
