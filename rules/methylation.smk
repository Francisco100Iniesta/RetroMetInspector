rule process_detectedETs:
    conda:
        "../metcheck.yaml"
    log:
        "logs/metcheck/{sample}.log"
    input:
        data = rules.analysis_preparatory.output.insertionsTable
    params:
        met_window = config["MetWindows"],
        min_met_diff = config["MINmetDiff"],
        min_reads = config["MinreadsSupport"],
        samples = list(config["samples"].keys())
    output:
        right_bed = "metconstructor/{sample}.right.bed",
        left_bed = "metconstructor/{sample}.left.bed"
    script:
        "../scripts/processData.R"


rule methylartist_windows:
    conda:
        "../metcheck.yaml"
    log:
        "logs/metcheck/{sample}.log"
    input:
        right_bed = "metconstructor/{sample}.right.bed",
        left_bed  = "metconstructor/{sample}.left.bed",
        bam       = "alns/{sample}.bam",
        bai       = "alns/{sample}.bam.bai",
        genome    = config["referenceGenome"]
    output:
        left_tsv  = "metconstructor/{sample}.left.tsv",
        right_tsv = "metconstructor/{sample}.right.tsv"
    shell:
        """
        methylartist segmeth \
            -b {input.bam} \
            -i {input.left_bed} \
            --ref {input.genome} \
            --motif CG \
            -p 16 \
            -o {output.left_tsv} \
            --phased

        methylartist segmeth \
            -b {input.bam} \
            -i {input.right_bed} \
            --ref {input.genome} \
            --motif CG \
            -p 16 \
            -o {output.right_tsv} \
            --phased
        """


rule process_methylartist:
    conda:
        "../metcheck.yaml"
    log:
        "logs/process_methylartist/{sample}.log"
    input:
        data = rules.analysis_preparatory.output.insertionsTable,
        left_tsv  = "metconstructor/{sample}.left.tsv",
        right_tsv = "metconstructor/{sample}.right.tsv"
    params:
        met_window = config["MetWindows"],
        min_met_diff = config["MINmetDiff"],
        min_reads = config["MinreadsSupport"],
        samples = list(config["samples"].keys())
    output:
        AllInfo = "rds/AllInfo{sample}.rds",
        TOTAL_INS = "rds/TOTAL_INS{sample}.rds",
        Query_bcftools = "metconstructor/Bcftools_query{sample}.bed",
        Query_bcftools2 = "metconstructor/Bcftools_query2{sample}.bed",
        BedCandidatas = "methylartist_locus/BedCandidatas{sample}.bed",
        BedSecuencias = "methylartist_locus/BedSecuencias{sample}.bed"
    script:
        "../scripts/processMethylartist.R"


rule bcftools_query:
    conda:
        "../metcheck.yaml"
    log:
        "logs/bcftools/{sample}.log"
    input:
        Query_bcftools = "metconstructor/Bcftools_query{sample}.bed",
        Query_bcftools2 = "metconstructor/Bcftools_query2{sample}.bed",
        vcf = "variants/sniffles2/{sample}.sniffles2.vcf.gz",
        csi = "variants/sniffles2/{sample}.sniffles2.vcf.gz.csi"
    output:
        Query = "metconstructor/Retro_Query{sample}.tsv",
        Query2 = "metconstructor/Retro_Query2{sample}.tsv"
    shell:
        """
        bcftools query -R {input.Query_bcftools} \
            -f '%CHROM\\t%POS\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%INFO/PHASE\\t[%GT\\t%GQ\\t%DR\\t%DV]\\n' \
            {input.vcf} > {output.Query}

        bcftools query -R {input.Query_bcftools2} \
            -f '%CHROM\\t%POS\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%INFO/PHASE\\t[%GT\\t%GQ\\t%DR\\t%DV]\\n' \
            {input.vcf} > {output.Query2}
        """


rule process_bcftools:
    conda:
        "../metcheck.yaml"
    log:
        "logs/process_bcftools/{sample}.log"
    input:
        AllInfo = "rds/AllInfo{sample}.rds",
        TOTAL_INS = "rds/TOTAL_INS{sample}.rds",
        Query = "metconstructor/Retro_Query{sample}.tsv",
        Query2 = "metconstructor/Retro_Query2{sample}.tsv"
    params:
        met_window = config["MetWindows"],
        min_met_diff = config["MINmetDiff"],
        min_reads = config["MinreadsSupport"]
    output:
        FinalTable = "rds/FinalMetTable{sample}.rds",
        FinalTable2 = "rds/FinalMetTable2{sample}.rds"
    script:
        "../scripts/processBcftoolsQuery.R"


rule GenerateReport:
    conda:
        "../metcheck.yaml"
    input:
        rmd = "scripts/ReportMethyl.Rmd",
        final2 = "rds/FinalMetTable2{sample}.rds",
        final = "rds/FinalMetTable{sample}.rds",
        ref_plots = "methylartist_plots_ref/{sample}",
        custom_plots = "customref_run/{sample}"
    output:
        report = "reports/meth_report.{sample}.html"
    log:
        "logs/GenerateReport/{sample}.log"
    params:
        met_window = config["MetWindows"],
        min_met_diff = config["MINmetDiff"],
        min_reads = config["MinreadsSupport"]
    script:
        "../scripts/ReportMethyl.Rmd"
