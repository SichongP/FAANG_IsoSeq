import glob, os
localrules: rename_files, unzip_fastq, moveFiles, get_counts_post_confirm_collapse
configfile: "config.yaml"

wildcard_constraints:
    file="[a-zA-Z0-9_]+"

SUBREADS = {"AH2_Lamina": "Pool_1", "AH2_Ovary": "Pool_1",
            "AH1_Ovary": "Pool_2", "AH1_Lung": "Pool_2",
            "AH1_Liver": "Pool_3", "AH1_Heart": "Pool_3",
            "AH1_Muscle": "Pool_4", "AH1_Skin": "Pool_4",
            "AH4_Lamina": "Pool_5", "AH4_Liver": "Pool_5",
            "AH4_Testes": "Pool_6", "AH4_Heart": "Pool_6",
            "AH4_Muscle": "Pool_7", "AH3_Testes": "Pool_7",
            "AH3_Lung": "Pool_8", "AH3_Skin": "Pool_8"}

FILES = []
SAMPLES = {}
TISSUES = {}

for file in glob.glob("data/Iso-seq/*ccs*bam"):
    file = os.path.split(file)[1].split('.')[0]
    if file.startswith('Parietal'):
        continue
    FILES.append(file)
    tissue, sample = file.split('_')
    SAMPLES[sample] = [tissue] if sample not in SAMPLES.keys() else SAMPLES[sample] + [tissue]
    TISSUES[tissue] = [sample] if tissue not in TISSUES.keys() else TISSUES[tissue] + [sample]

ALL_RNASEQ = []
for file in glob.glob("data/RNAseq_fq/*R1*gz"):
    file = '_'.join(os.path.split(file)[1].split('.')[0].split('_')[:2])
    ALL_RNASEQ.append(file)

rule all:
    input:
        expand("Results/salmon/{step}/{file}/quant.sf", file = ALL_RNASEQ, step = ['sqanti', 'chained', 'annotated', 'refseq', 'ensembl']),
        expand("Results/salmon/{step}/{file}/quant.sf", file = FILES, step = ['collapsed']),
        #"Results/Annotated/all_samples.gff3"
        "Results/figures/profile/all_samples.tab",
        "Results/figures/ATACseq_profile/all_samples.tab",
        expand("Results/Enrichment/{stage}/{file}/all_samples.tsv", file =  FILES, stage = ['Annotated', 'collapsed', 'Chained', 'SQANTI3', 'RefSeq', 'ensembl']),
        expand("Results/mapping/stats/{ref}_{file}.stats.txt", ref = ['refSeqTranscriptomeOnly', 'IsoSeqTranscriptomeOnly', 'ensembleTranscriptomeOnly'], file = ALL_RNASEQ),
        expand("Results/Rarefaction/{file}.{by}.txt", by = ['pbgene', 'pbid'], file = FILES)
        
rule refine:
    # This rule removes polyA tails and concatemers from ccs reads to generate FLNC (full length non-chimeric) transcripts
    # It also removes any reads without at least 20bp of polyA tails
    input: "data/Iso-seq/{file}.ccs.bam"
    output: temp("Results/FLNC/{file}.flnc.bam")
    resources:
        cpus = 4, cpus_bmm = 4,
        mem_mb = 20000, mem_mb_bmm = 20000,
        time_min = 120
    params: partition = 'bmm'
    conda: "envs/pacbio.yaml"
    shell:
     """
     isoseq3 refine {input} data/barcodes.fa {output} --require-polya -j {resources.cpus}
     """
     
rule cluster:
    # This rule clusters reads from same transcripts into a single clustered transcript
    input: "Results/FLNC/{file}.flnc.bam"
    output:
        bam = "Results/Clustered/{file}.clustered.bam",
        fa = "Results/Clustered/{file}.clustered.hq.fasta.gz",
        report = "Results/Clustered/{file}.clustered.cluster_report.csv"
    resources:
        cpus = 2, cpus_bmm = 2,
        mem_mb = 50000, mem_mb_bmm = 50000,
        time_min = 1800
    params: partition = 'bmm'
    conda: "envs/pacbio.yaml"
    shell:
     """
     isoseq3 cluster {input} {output} --use-qvs --verbose -j {resources.cpus}
     """

rule map_long_read_ec3:
    # Map polished long reads to reference genome
    input: "Results/Clustered/{file}.clustered.hq.fasta.gz"
    output: temp("Results/minimap2/{file}.bam")
    conda: "envs/minimap2.yaml"
    resources:
        cpus = 8, cpus_bmm = 8,
        mem_mb = 200000, mem_mb_bmm = 200000,
        time_min = 4200
    params: partition = 'bmm'
    shell:
     """
     minimap2 -t {resources.cpus} -H -ax splice:hq -uf {config[ref_ec3]} {input} | \
       samtools view -bh -o {output}
     """

rule sort_bam:
    # This rule sorts bam files to prepare for collapsing
    input: "Results/minimap2/{file}.bam"
    output: "Results/minimap2/{file}.sorted.bam"
    conda: "envs/samtools.yaml"
    resources:
        cpus = 1,
        mem_mb = 4000,
        time_min = 30
    params: partition = 'bmm'
    shell:
     """
     samtools sort -o {output} {input}
     samtools index {output}
     """

rule unzip_fastq:
    #This rule unzips fasta file produced by Cluster step as Cupcake only takes unzipped file
    input: "Results/Clustered/{file}.clustered.hq.fasta.gz"
    output: temp("Results/Clustered/{file}.clustered.hq.fasta")
    shell:
     """
     gunzip -k {input}
     """

rule collapse_reads:
    # Collapse similar reads to reduce redundancies
    input: 
        bam = "Results/minimap2/{file}.sorted.bam",
        fa = "Results/Clustered/{file}.clustered.hq.fasta"
    output: "Results/Collapsed/{file}.collapsed.gff", "Results/Collapsed/{file}.collapsed.group.txt"
    params: partition = 'bmm', out_prefix = "Results/Collapsed/{file}"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 150000, mem_mb_bmm = 150000,
        time_min = 2200
    conda: "envs/cupcake.yaml"
    shell:
     """
     collapse_isoforms_by_sam.py --input {input.fa} \
      -b {input.bam} -o {params.out_prefix}
     """

rule getAbundance:
    # This rule counts supporting reads for each collapsed transcript (not merged)
    input:
        gff = "Results/Collapsed/{file}.collapsed.gff",
        counts = "Results/Clustered/{file}.clustered.cluster_report.csv"
    output: "Results/Collapsed/{file}.collapsed.abundance.txt"
    params:
        input_prefix = lambda wildcards: "Results/Collapsed/{file}.collapsed".format(file = wildcards.file),
        partition = 'bmm'
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 10000, mem_mb_bmm = 10000,
        time_min = 30
    conda: "envs/cupcake.yaml"
    shell:
     """
     get_abundance_post_collapse.py {params.input_prefix} {input.counts}
     """

rule filterByAbundance:
    # This rule filters collapsed transcripts by their abundance (minimum 2 counts)
    input:
        gff = "Results/Collapsed/{file}.collapsed.gff",
        counts = "Results/Collapsed/{file}.collapsed.abundance.txt"
    output: "Results/Collapsed/{file}.collapsed.min_fl_2.gff"
    params:
        input_prefix = lambda wildcards: "Results/Collapsed/{file}.collapsed".format(file = wildcards.file),
        partition = 'med2'
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 2048, mem_mb_bmm = 2048,
        time_min = 30
    conda: "envs/cupcake.yaml"
    shell:
     """
     filter_by_count.py --min_count 2 --dun_use_group_count {params.input_prefix}
     """

rule filter5primedegraded:
    # This rule removes reads if 5 prime is deemed degraded (shorter 5' with identical 3' sequences)
    input: "Results/Collapsed/{file}.collapsed.min_fl_2.gff"
    output: "Results/Collapsed/{file}.collapsed.min_fl_2.filtered.gff", "Results/Collapsed/{file}.collapsed.min_fl_2.filtered.abundance.txt", "Results/Collapsed/{file}.collapsed.min_fl_2.filtered.rep.fa"
    params:
        input_prefix = lambda wildcards: "Results/Collapsed/{file}.collapsed.min_fl_2".format(file = wildcards.file),
        partition = 'med2'
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 2048, mem_mb_bmm = 2048,
        time_min = 30
    conda: "envs/cupcake.yaml"
    shell:
     """
     filter_away_subset.py {params.input_prefix}
     """

rule moveFiles:
    # This rule reorganizes files so they can be chained together
    # See https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step#chain for details
    input:
        expand("Results/Collapsed/{file}.collapsed.min_fl_2.filtered.gff", file = FILES),
        expand("Results/Collapsed/{file}.collapsed.min_fl_2.filtered.abundance.txt", file = FILES),
        expand("Results/Collapsed/{file}.collapsed.min_fl_2.filtered.rep.fa", file = FILES),
        expand("Results/Collapsed/{file}.collapsed.group.txt", file = FILES)
    output: 
        expand("Results/Collapsed/{file}/collapsed.min_fl_2.filtered.gff", file = FILES),
        expand("Results/Collapsed/{file}/collapsed.min_fl_2.filtered.abundance.txt", file = FILES),
        expand("Results/Collapsed/{file}/collapsed.min_fl_2.filtered.rep.fa", file = FILES),
        expand("Results/Collapsed/{file}/collapsed.group.txt", file = FILES),
        config = "Results/Chained/chaining.config"

    params: files = lambda wildcards: FILES
    shell:
     """
     for name in {params.files}
     do
         mkdir -p Results/Collapsed/${{name}}
         for file in Results/Collapsed/*${{name}}.collapsed.{{min_fl_2.filtered.gff,min_fl_2.filtered.abundance.txt,min_fl_2.filtered.rep.fa,group.txt}}
         do
             mv $file Results/Collapsed/${{name}}/$(echo $(basename $file) | sed -e "s/${{name}}\.//")
         done
         echo "SAMPLE=${{name}};Results/Collapsed/${{name}}" >> {output.config}
     done
     echo "" >> {output.config}
     echo "GROUP_FILENAME=collapsed.group.txt" >> {output.config}
     echo "GFF_FILENAME=collapsed.min_fl_2.filtered.gff" >> {output.config}
     echo "COUNT_FILENAME=collapsed.min_fl_2.filtered.abundance.txt" >> {output.config}
     echo  "FASTQ_FILENAME=collapsed.min_fl_2.filtered.rep.fastq" >> {output.config}
     """

rule fa2fq:
    # Generate generic fastq files from fasta files
    input: "Results/Collapsed/{file}/collapsed.min_fl_2.filtered.rep.fa"
    output: "Results/Collapsed/{file}/collapsed.min_fl_2.filtered.rep.fastq"
    params: partition = 'bmm'
    conda: "envs/cupcake.yaml"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 4000, mem_mb_bmm = 4000,
        time_min = 600
    shell:
     """
     fa2fq.py {input}
     """

rule chainSamples:
    # This rule chains transcripts from replicates 
    input:
        group = expand("Results/Collapsed/{file}/collapsed.group.txt", file = FILES),
        gff = expand("Results/Collapsed/{file}/collapsed.min_fl_2.filtered.gff", file = FILES),
        counts = expand("Results/Collapsed/{file}/collapsed.min_fl_2.filtered.abundance.txt", file = FILES),
        fq = expand("Results/Collapsed/{file}/collapsed.min_fl_2.filtered.rep.fastq", file = FILES),
        config = "Results/Chained/chaining.config"
    output:
        gff = "Results/Chained/all_samples.gff",
        fq = "Results/Chained/all_samples.rep.fastq",
        id = "Results/Chained/all_samples_ids.txt",
        counts = "Results/Chained/all_samples_count.txt"
    params: partition = 'bmm'
    conda: "envs/cupcake.yaml"
    resources:
        cpus = 10, cpus_bmm = 10,
        mem_mb = 350000, mem_mb_bmm = 350000,
        time_min = 3000
    shell:
     """
     chain_samples.py --cpus {resources.cpus} --dun-merge-5-shorter {input.config} count_fl 
     mv all_samples.chained_count.txt {output.counts}
     mv all_samples.chained.gff {output.gff}
     mv all_samples.chained_ids.txt {output.id}
     mv all_samples.chained.rep.fq {output.fq}
     """

rule validateTranscripts:
    # This rule looks through collapsed transcripts and only retain those present in at least 2 of the 4 samples
    input: 
        counts = "Results/Chained/all_samples_count.txt",
        gff = "Results/Chained/all_samples.gff",
        fq = "Results/Chained/all_samples.rep.fastq"        
    output: 
        gff = "Results/Confirmed/all_samples.gff",
        fa = "Results/Confirmed/all_samples.fa"
    conda: "envs/pandas.yaml"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 150000, mem_mb_bmm = 150000,
        time_min = 120
    params: partition = 'bmh'
    script: "scripts/validateTranscripts.py"


rule map_long_read_ec3_confirmed:
    # Map confirmed long reads to reference genome
    input: "Results/Confirmed/all_samples.fa"
    output: "Results/Confirmed/all_samples.bam"
    conda: "envs/minimap2.yaml"
    resources:
        cpus = 4, cpus_bmm = 4,
        mem_mb = 450000, mem_mb_bmm = 450000,
        time_min = 3200
    params: partition = 'bmm'
    shell:
     """
     minimap2 -t {resources.cpus} -H -ax splice:hq -uf {config[ref_ec3]} {input} | \
       samtools view -bh -o {output}
     """
    
rule sort_bam_confirmed:
    # This rule sorts confirmed bam files to prepare for collapsing
    input: "Results/Confirmed/all_samples.bam"
    output: "Results/Confirmed/all_samples.sorted.bam"
    conda: "envs/samtools.yaml"
    resources:
        cpus = 1,
        mem_mb = 4000,
        time_min = 30
    params: partition = 'bmm'
    shell:
     """
     samtools sort -o {output} {input}
     samtools index {output}
     """

rule collapse_reads_confirmed:
    # Collapse reads again after chaining and filtering since chaining does not merge 5' degraded transcripts
    input: 
        bam = "Results/Confirmed/all_samples.sorted.bam",
        fq = "Results/Confirmed/all_samples.fa"
    output: "Results/Confirmed/all_samples.collapsed.gff.unfuzzy", "Results/Confirmed/all_samples.collapsed.rep.fa", "Results/Confirmed/all_samples.collapsed.group.txt.unfuzzy"
    params: partition = 'bmm', out_prefix = "Results/Confirmed/all_samples"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 80000, mem_mb_bmm = 80000,
        time_min = 1200
    conda: "envs/cupcake.yaml"
    shell:
     """
     collapse_isoforms_by_sam.py --input {input.fq} \
      -b {input.bam} -o {params.out_prefix} --max_3_diff 0 --max_5_diff 100
     """

rule get_counts_post_confirm_collapse:
    # Get FL counts of reads after chaining
    input:
        chained_count = "Results/Chained/all_samples_count.txt",
        collapsed_group = "Results/Confirmed/all_samples.collapsed.group.txt.unfuzzy"
    output: collapsed_count = "Results/Confirmed/all_samples.collapsed.count.txt"
    conda: "envs/pandas.yaml"
    script: "scripts/get_counts_post_collapse.py"

rule kallisto_index:
    # Build kallisto index for quantification
    input: "Results/Confirmed/all_samples.collapsed.rep.fa"
    output: "Results/Kallisto_index/all_samples.confirmed.collapsed.idx"
    conda: "envs/kallisto.yaml"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 10000, mem_mb_bmm = 10000,
        time_min = 120
    params: partition = 'bmm'
    shell:
     """
     kallisto index -i {output} {input}
     """

rule kallisto:
    # Quantify chained transcriptome using Kallisto and short read RNAseq data
    input:
        r1 = "data/RNAseq_fq/{file}_R1.fastq.gz",
        r2 = "data/RNAseq_fq/{file}_R2.fastq.gz",
        idx = "Results/Kallisto_index/all_samples.confirmed.collapsed.idx",
        gff = "Results/Confirmed/all_samples.collapsed.gff.unfuzzy"
    output: "Results/Kallisto/Abundance/{file}.abundance.tsv"
    conda: "envs/kallisto.yaml"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 10000, mem_mb_bmm = 10000,
        time_min = 60
    params:
        out_dir = lambda wildcards: "Results/Kallisto/{}/".format(wildcards.file),
        partition = "bmm"
    shell:
     """
     kallisto quant -i {input.idx} -o {params.out_dir} --bias -t {resources.cpus} --gtf {input.gff} {input.r1} {input.r2}
     cp {params.out_dir}abundance.tsv Results/Kallisto/Abundance/{wildcards.file}.abundance.tsv
     """

rule sqanti3:
    # This rule runs sqanti3 on confirmed transcripts
    input:
        gff = "Results/Confirmed/all_samples.collapsed.gff.unfuzzy",
        ref_gtf = "data/ensembl.ref.gtf",
        counts = "Results/Confirmed/all_samples.collapsed.count.txt",
        rnaseq = "data/RNAseq_data.fofn",
        exp_data = expand("Results/Kallisto/Abundance/{file}.abundance.tsv", file = FILES)
        #sj = "data/SJ/*ParietalCortex*tab"
    output: "Results/SQANTI3/all_samples_SQANTI3_report.pdf", "Results/SQANTI3/all_samples_classification.txt", "Results/SQANTI3/all_samples_corrected.gtf", "Results/SQANTI3/all_samples_junctions.txt", "Results/SQANTI3/all_samples_corrected.fasta"
    params:
        partition = 'bmm',
        sj = "data/SJ/",
        sr_bam = "data/SR_bam/",
        exp = "Results/Kallisto/Abundance/"
    conda: "envs/sqanti3.yaml"
    resources:
        cpus = 22, cpus_bmm = 22,
        mem_mb = 120000, mem_mb_bmm = 120000,
        time_min = 2800
    shell:
     """
     #gunzip -c {config[ref_ec3]} > data/sqanti3_genome.fa
     export PYTHONPATH=/home/pengsc/projects/Iso-seq/tools/cDNA_Cupcake/sequence/
     mkdir -p Results/SQANTI3/
     ../tools/SQANTI3/sqanti3_qc.py {input.gff} {input.ref_gtf} data/sqanti3_genome.fa \
     --fl_count {input.counts} -c {params.sj} --SR_bam {params.sr_bam} --expression {params.exp} \
     -t {resources.cpus} -o all_samples -d Results/SQANTI3/ --report pdf
     """
     
rule annotate_sq_tx:
    # Filter nonsense-mediated decay transcripts, and annotate final transcripts based on known and homologous genes 
    input:
        classification = "Results/SQANTI3/all_samples_classification.txt",
        sq_gtf = "Results/SQANTI3/all_samples_corrected.gtf",
        junctions = "Results/SQANTI3/all_samples_junctions.txt",
        salmon_cluster = "Results/salmon/sqanti/salmon_index/duplicate_clusters.tsv",
        salmon_counts = expand("Results/salmon/sqanti/{file}/quant.sf", file = FILES)
    output: gff = "Results/Annotated/all_samples.gff3"
    conda: "envs/pandas.yaml"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 20000, mem_mb_bmm = 20000,
        time_min = 300
    params: partition = 'bmm'
    script: "scripts/filter_consolidate_sq_tx.py"
    
rule plotEnrichment_filtered_combined:
    # Plot fraction of reads in transcripts (FRiT) against final transcriptome
    input: 
        gff = "Results/Annotated/all_samples.gff3",
        bam = "data/SR_bam/{file}.markDup.bam"
    output: 
        pdf = "Results/Enrichment/Annotated/{file}/all_samples.pdf",
        tsv = "Results/Enrichment/Annotated/{file}/all_samples.tsv"
    conda: "envs/deeptools.yaml"
    resources:
        cpus = 20, cpus_bmm = 20,
        mem_mb = 20000, mem_mb_bmm = 20000,
        time_min = 900
    params: partition = 'bmm'
    shell:
     """
     gffread -T -F -o Results/Enrichment/Annotated/temp_{wildcards.file}_all_samples.gtf {input.gff}
     plotEnrichment -b {input.bam} --BED Results/Enrichment/Annotated/temp_{wildcards.file}_all_samples.gtf -o {output.pdf} --outRawCounts {output.tsv} --smartLabels -p {resources.cpus} --samFlagExclude 1024 --minMappingQuality 20
     rm Results/Enrichment/Annotated/temp_{wildcards.file}_all_samples.gtf
     """
     
rule plotEnrichment_sqanti:
    # Plot fraction of reads in transcripts (FRiT) against pre-filtered transcriptome
    input: 
        gff = "Results/SQANTI3/all_samples_corrected.gtf",
        bam = "data/SR_bam/{file}.markDup.bam"
    output: 
        pdf = "Results/Enrichment/SQANTI3/{file}/all_samples.pdf",
        tsv = "Results/Enrichment/SQANTI3/{file}/all_samples.tsv"
    conda: "envs/deeptools.yaml"
    resources:
        cpus = 20, cpus_bmm = 20,
        mem_mb = 20000, mem_mb_bmm = 20000,
        time_min = 900
    params: partition = 'med2'
    shell:
     """
     plotEnrichment -b {input.bam} --BED {input.gff} -o {output.pdf} --outRawCounts {output.tsv} --smartLabels -p {resources.cpus} --samFlagExclude 1024 --minMappingQuality 20
     """
     
rule plotEnrichment_chained:
    # Plot fraction of reads in transcripts (FRiT) against chained transcriptome
    # We assess FRiT in three steps to measure losses of signal during filtering
    input: 
        gff = "Results/Chained/all_samples.gff",
        bam = "data/SR_bam/{file}.markDup.bam"
    output: 
        pdf = "Results/Enrichment/Chained/{file}/all_samples.pdf",
        tsv = "Results/Enrichment/Chained/{file}/all_samples.tsv"
    conda: "envs/deeptools.yaml"
    resources:
        cpus = 20, cpus_bmm = 20,
        mem_mb = 20000, mem_mb_bmm = 20000,
        time_min = 900
    params: partition = 'high2'
    shell:
     """
     gffread -T -F -o Results/Enrichment/Chained/temp_{wildcards.file}.gtf {input.gff}
     plotEnrichment -b {input.bam} --BED Results/Enrichment/Chained/temp_{wildcards.file}.gtf -o {output.pdf} --outRawCounts {output.tsv} --smartLabels -p {resources.cpus} --samFlagExclude 1024 --minMappingQuality 20
     rm Results/Enrichment/Chained/temp_{wildcards.file}.gtf
     """
     
rule plotEnrichment_RefSeq:
    # Plot fraction of reads in transcripts (FRiT) against RefSeq transcriptome
    input: 
        gff = "/home/pengsc/reference/equcab3/RefSeq/AnnotationRelease103/GCF_002863925.1_EquCab3.0_genomic.renamed.gff",
        bam = "data/SR_bam/{file}.markDup.bam"
    output: 
        pdf = "Results/Enrichment/RefSeq/{file}/all_samples.pdf",
        tsv = "Results/Enrichment/RefSeq/{file}/all_samples.tsv"
    conda: "envs/deeptools.yaml"
    resources:
        cpus = 20, cpus_bmm = 20,
        mem_mb = 20000, mem_mb_bmm = 20000,
        time_min = 900
    params: partition = 'high2'
    shell:
     """
     gffread -T -F -o Results/Enrichment/RefSeq/temp_{wildcards.file}.gtf {input.gff}
     plotEnrichment -b {input.bam} --BED Results/Enrichment/RefSeq/temp_{wildcards.file}.gtf -o {output.pdf} --outRawCounts {output.tsv} --smartLabels -p {resources.cpus} --samFlagExclude 1024 --minMappingQuality 20
     rm Results/Enrichment/RefSeq/temp_{wildcards.file}.gtf
     """
     
rule plotEnrichment_ensembl:
    # Plot fraction of reads in transcripts (FRiT) against Ensembl transcriptome
    input: 
        gff = "/home/pengsc/reference/equcab3/ensemble/Equus_caballus.EquCab3.0.102.gff3.renamed",
        bam = "data/SR_bam/{file}.markDup.bam"
    output: 
        pdf = "Results/Enrichment/ensembl/{file}/all_samples.pdf",
        tsv = "Results/Enrichment/ensembl/{file}/all_samples.tsv"
    conda: "envs/deeptools.yaml"
    resources:
        cpus = 20, cpus_bmm = 20,
        mem_mb = 20000, mem_mb_bmm = 20000,
        time_min = 900
    params: partition = 'high2'
    shell:
     """
     gffread -T -F -o Results/Enrichment/ensembl/temp_{wildcards.file}.gtf {input.gff}
     plotEnrichment -b {input.bam} --BED Results/Enrichment/ensembl/temp_{wildcards.file}.gtf -o {output.pdf} --outRawCounts {output.tsv} --smartLabels -p {resources.cpus} --samFlagExclude 1024 --minMappingQuality 20
     rm Results/Enrichment/ensembl/temp_{wildcards.file}.gtf
     """

rule plotEnrichment_collapsed:
    input: 
        gff = "Results/Collapsed/{file}/collapsed.min_fl_2.filtered.gff",
        bam = "data/SR_bam/{file}.markDup.bam"
    output: 
        pdf = "Results/Enrichment/collapsed/{file}/all_samples.pdf",
        tsv = "Results/Enrichment/collapsed/{file}/all_samples.tsv"
    conda: "envs/deeptools.yaml"
    resources:
        cpus = 20, cpus_bmm = 20,
        mem_mb = 20000, mem_mb_bmm = 20000,
        time_min = 900
    params: partition = 'high2'
    shell:
     """
     gffread -T -F -o Results/Enrichment/collapsed/temp_{wildcards.file}.gtf {input.gff}
     plotEnrichment -b {input.bam} --BED Results/Enrichment/collapsed/temp_{wildcards.file}.gtf -o {output.pdf} --outRawCounts {output.tsv} --smartLabels -p {resources.cpus} --samFlagExclude 1024 --minMappingQuality 20
     rm Results/Enrichment/collapsed/temp_{wildcards.file}.gtf
     """
    
    
include: "rules/salmon.smk"
include: "rules/RNAcov.smk"
include: "rules/ATACcov.smk"
include: "rules/transcriptome_alignment.smk"
include: "rules/Rarefaction.smk"