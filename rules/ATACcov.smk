ATAC_SAMPLES = ['AH1_Adipose', 'AH1_Liver', 'AH1_Ovary', 'AH2_Heart', 'AH2_Lung', 'AH2_ParietalCortex', 'AH3_Lamina', 'AH3_Muscle', 'AH4_Adipose', 'AH4_Liver', 'AH4_ParietalCortex', 'AH1_Heart', 'AH1_Lung', 'AH1_ParietalCortex', 'AH2_Lamina', 'AH2_Muscle',  'AH3_Adipose', 'AH3_Liver', 'AH3_ParietalCortex', 'AH4_Heart', 'AH4_Lung', 'AH4_Testes', 'AH1_Lamina', 'AH1_Muscle', 'AH2_Adipose', 'AH2_Liver', 'AH2_Ovary', 'AH3_Heart', 'AH3_Lung', 'AH3_Testes', 'AH4_Lamina', 'AH4_Muscle']

rule ATAC_bamCoverage:
    input: "data/ATACseq_bam/{file}.bam"
    output: "Results/ATACseq/bigwig/{file}.bw"
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

rule ATAC_computeMatrix:
    input: 
        bws = expand("Results/ATACseq/bigwig/{file}.bw", file = ATAC_SAMPLES),
        gff = "Results/Annotated/all_samples.gff3"
    output: "Results/ATACseq/matrix/all_samples.mat.gz"
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
     
rule ATAC_plotProfile:
    input: 
        mat = "Results/ATACseq/matrix/all_samples.mat.gz"
    output:
        png = "Results/figures/ATACseq_profile/all_samples.png",
        tab = "Results/figures/ATACseq_profile/all_samples.tab"
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