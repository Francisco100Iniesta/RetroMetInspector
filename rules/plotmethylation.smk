rule plotREF:
    conda:
        "../metcheck.yaml"
    input:
        bed = "methylartist_locus/BedSecuencias{sample}.bed",
        ref = config["referenceGenome"],
        bam = "alns/{sample}.bam",
        bai = "alns/{sample}.bam.bai",
        gtf = config["humanGTF"]
    output:
        directory("methylartist_plots_ref/{sample}")
    log:
        "logs/process_methylartist/{sample}.log"
    shell:
        r"""
        mkdir -p {output}
        scripts/methPlot_ref.sh {input.bed} {input.ref} {input.bam} {input.gtf} {output} > {log} 2>&1
        """
rule custom_and_plots:
    conda:
        "../metcheck.yaml"
    input:
        bed = "methylartist_locus/BedSecuencias{sample}.bed",
        ref = config["referenceGenome"],
        bamselect = "methylartist_locus/BedCandidatas{sample}.bed",
        bam = "alns/{sample}.bam"
    output:
        directory("customref_run/{sample}")
    log:
        "logs/customref_run/{sample}.log"
    threads:
        max(1, config.get("threads", 1) - 1)
    shell:
        r"""
        mkdir -p {output}
        SLURM_CPUS_PER_TASK={threads} scripts/methPlot_customref.sh {input.bed} {input.ref} {input.bamselect} {input.bam} {output} > {log} 2>&1
        """
