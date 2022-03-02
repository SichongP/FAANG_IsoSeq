rule bamCoverage:
    input: "data/SR_bam/{file}.markDup.bam"
    output: "Results/RNAseq/bigwig/{file}.bw"
    params: partition = "bmm"
    conda: "../envs/deeptools.yaml"
    resources:
        cpus = 6, cpus_bmm = 6,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300
    shell:
     """
     bamCoverage -b {input} -o {output} -of "bigwig" -bs 10 --effectiveGenomeSize {config[effective_size]} --normalizeUsing "RPGC" --exactScaling --ignoreDuplicates \
     --minMappingQuality 30 
     """

rule computeMatrix:
    input: 
        bws = expand("Results/RNAseq/bigwig/{file}.bw", file = FILES),
        gff = "Results/Annotated/all_samples.gff3"
    output: "Results/RNAseq/matrix/all_samples.mat.gz"
    params: partition = "bmm"
    conda: "../envs/deeptools.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 70000, mem_mb_bmm = 70000,
        time_min = 900
    shell:
     """
     awk '{{if ($3=="mRNA") print $1"\t"int($4)"\t"int($5)"\ttranscript\t0\t"$7}}' {input.gff} > {input.gff}.gene.bed
     computeMatrix reference-point --referencePoint "TSS" -S {input.bws} -R {input.gff}.gene.bed --beforeRegionStartLength 1000 --nanAfterEnd --afterRegionStartLength 1000 --skipZeros -o {output} -p {resources.cpus} --maxThreshold 3000
     """
     
rule plotProfile:
    input: 
        mat = "Results/RNAseq/matrix/all_samples.mat.gz"
    output:
        png = "Results/figures/profile/all_samples.png",
        tab = "Results/figures/profile/all_samples.tab"
    params: partition = "bmm"
    conda: "../envs/deeptools.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 70000, mem_mb_bmm = 70000,
        time_min = 900
    shell:
     """
     plotProfile -m {input.mat} -o {output.png} --outFileNameData {output.tab} --dpi 300 --plotType=fill
     """