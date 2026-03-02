rule process_compare:
    conda:
        "../metcheck.yaml"
    log:
        "logs/metcheck/compare_{sample1}_vs_{sample2}.log"
    input:
        data=rules.analysis_preparatory.output.insertionsTable
    params:
        met_window=config["MetWindows"],
        min_met_diff=config["MINmetDiff"],
        min_reads=config["MinreadsSupport"],
        samples=list(config["samples"].keys()),
        soi=["{sample1}", "{sample2}"]
    output:
        right_bed="metconstructor_compare/{sample1}_vs_{sample2}_compare.right.bed",
        left_bed="metconstructor_compare/{sample1}_vs_{sample2}_compare.left.bed",
        unshared_insertion="rds/{sample1}_vs_{sample2}_unshared.rds",
        bat_file="igv/{sample1}_vs_{sample2}_compare.bat"  # <-- Añadido en carpeta igv
    script:
        "../scripts/processData_compare.R"


rule methylartist_windows_compare:
    conda:
        "../metcheck.yaml"
    log:
        "logs/methylartist/{sample1}_vs_{sample2}.log"
    threads:
        max(1, int(config.get("threads", 1)) - 1)
    input:
        right_bed="metconstructor_compare/{sample1}_vs_{sample2}_compare.right.bed",
        left_bed="metconstructor_compare/{sample1}_vs_{sample2}_compare.left.bed",
        bam1="alns/{sample1}.bam",
        bai1="alns/{sample1}.bam.bai",
        bam2="alns/{sample2}.bam",
        bai2="alns/{sample2}.bam.bai",
        genome=config["referenceGenome"]
    output:
        left_tsv1="metconstructor_compare/{sample1}_vs_{sample2}_compare_S1.left.tsv",
        right_tsv1="metconstructor_compare/{sample1}_vs_{sample2}_compare_S1.right.tsv",
        left_tsv2="metconstructor_compare/{sample1}_vs_{sample2}_compare_S2.left.tsv",
        right_tsv2="metconstructor_compare/{sample1}_vs_{sample2}_compare_S2.right.tsv"
    shell:
        """
        methylartist segmeth \
            -b {input.bam1} \
            -i {input.left_bed} \
            --ref {input.genome} \
            --motif CG \
            -p {threads} \
            -o {output.left_tsv1} \
            --phased

        methylartist segmeth \
            -b {input.bam1} \
            -i {input.right_bed} \
            --ref {input.genome} \
            --motif CG \
            -p {threads} \
            -o {output.right_tsv1} \
            --phased

        methylartist segmeth \
            -b {input.bam2} \
            -i {input.left_bed} \
            --ref {input.genome} \
            --motif CG \
            -p {threads} \
            -o {output.left_tsv2} \
            --phased

        methylartist segmeth \
            -b {input.bam2} \
            -i {input.right_bed} \
            --ref {input.genome} \
            --motif CG \
            -p {threads} \
            -o {output.right_tsv2} \
            --phased
        """


rule process_Meth_compare:
    conda:
        "../metcheck.yaml"
    log:
        "logs/metcheck/compare_{sample1}_vs_{sample2}_meth.log"
    input:
        left_tsv1="metconstructor_compare/{sample1}_vs_{sample2}_compare_S1.left.tsv",
        right_tsv1="metconstructor_compare/{sample1}_vs_{sample2}_compare_S1.right.tsv",
        left_tsv2="metconstructor_compare/{sample1}_vs_{sample2}_compare_S2.left.tsv",
        right_tsv2="metconstructor_compare/{sample1}_vs_{sample2}_compare_S2.right.tsv",
        unshared_insertion="rds/{sample1}_vs_{sample2}_unshared.rds",
        FinalTable_S1="rds/FinalMetTable{sample1}.rds",
        FinalTable2_S1="rds/FinalMetTable2{sample1}.rds",
        FinalTable_S2="rds/FinalMetTable{sample2}.rds",
        FinalTable2_S2="rds/FinalMetTable2{sample2}.rds"
    params:
        met_window=config["MetWindows"],
        min_met_diff=config["MINmetDiff"],
        min_reads=config["MinreadsSupport"],
        samples=list(config["samples"].keys()),
        soi=["{sample1}", "{sample2}"]
    output:
        S1_ALL="rds/{sample1}_vs_{sample2}_S1_ALL.rds",
        S2_ALL="rds/{sample1}_vs_{sample2}_S2_ALL.rds",
        S1_filtered="rds/{sample1}_vs_{sample2}_S1_filtered.rds",
        S2_filtered="rds/{sample1}_vs_{sample2}_S2_filtered.rds"
    script:
        "../scripts/compare_all_meth.R"


rule GenerateMethReport_compare:
    conda:
        "../metcheck.yaml"
    input:
        S1_ALL="rds/{sample1}_vs_{sample2}_S1_ALL.rds",
        S2_ALL="rds/{sample1}_vs_{sample2}_S2_ALL.rds",
        S1_filtered="rds/{sample1}_vs_{sample2}_S1_filtered.rds",
        S2_filtered="rds/{sample1}_vs_{sample2}_S2_filtered.rds"
    output:
        report="UnsharedReport/Compare_Meth_Report{sample1}_vs_{sample2}.html"
    log:
        "logs/compare_Meth_Report{sample1}_vs_{sample2}.log"
    params:
        met_window=config["MetWindows"],
        min_met_diff=config["MINmetDiff"],
        min_reads=config["MinreadsSupport"],
        samples=list(config["samples"].keys()),
        soi=["{sample1}", "{sample2}"]
    script:
        "../scripts/ReportMethyl_compare.Rmd"