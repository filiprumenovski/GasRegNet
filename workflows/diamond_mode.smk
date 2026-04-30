rule all:
    input:
        "results/diamond/README.txt"


rule scaffold_marker:
    output:
        "results/diamond/README.txt"
    shell:
        "mkdir -p results/diamond && printf 'GasRegNet DIAMOND workflow placeholder\\n' > {output}"
