import glob, os
configfile: "config.yaml"
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
    input: expand("Results/salmon/merged/{file}/quant.sf", file = ALL_RNASEQ)

rule sqanti3:
    input:
        gtf = "Results/Stringtie_transcripts/stranded/faang.vs.refseq.annotated.stranded.gtf",
        ref_gtf = "data/refSeq.gtf",
        rnaseq = "data/RNAseq_data.fofn",
        exp_data = expand("Results/Kallisto_stringtie/Abundance/{file}.abundance.tsv", file = FILES)
        #sj = "data/SJ/*ParietalCortex*tab"
    output: "Results/SQANTI3_stringtie/all_samples_SQANTI3_report.pdf", "Results/SQANTI3_stringtie/all_samples_classification.txt", "Results/SQANTI3_stringtie/all_samples_corrected.gtf", "Results/SQANTI3_stringtie/all_samples_junctions.txt", "Results/SQANTI3_stringtie/all_samples_corrected.fasta"
    params:
        partition = 'bmm',
        sj = "data/SJ/",
        sr_bam = "data/SR_bam/",
        exp = "Results/Kallisto/Abundance/"
    conda: "envs/sqanti3.yaml"
    resources:
        cpus = 22, cpus_bmm = 22,
        mem_mb = 120000, mem_mb_bmm = 120000,
        time_min = 2800, partition = 'bmh'
    shell:
     """
     #gunzip -c {config[ref_ec3]} > data/sqanti3_genome.fa
     export PYTHONPATH=/home/pengsc/projects/Iso-seq/tools/cDNA_Cupcake/sequence/
     mkdir -p Results/SQANTI3_stringtie/
     ../tools/SQANTI3/sqanti3_qc.py {input.gtf} {input.ref_gtf} data/sqanti3_genome.fa \
     -c {params.sj} --SR_bam {params.sr_bam} --expression {params.exp} \
     -t {resources.cpus} -o all_samples -d Results/SQANTI3_stringtie/ --report pdf
     """
     
rule kallisto_index:
    # Build kallisto index for quantification
    input: "Results/Stringtie_transcripts/stranded/faang.vs.refseq.stranded.fa"
    output: "Results/Stringtie_transcripts/stranded/faang.vs.refseq.stranded.idx"
    conda: "envs/kallisto.yaml"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 4000, mem_mb_bmm = 4000,
        time_min = 120, partition = 'bmh'
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
        idx = "Results/Stringtie_transcripts/confirmed_all_samples.idx",
        gff = "Results/Stringtie_transcripts/confirmed_all_samples.gtf"
    output: "Results/Kallisto_stringtie/Abundance/{file}.abundance.tsv"
    conda: "envs/kallisto.yaml"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 10000, mem_mb_bmm = 10000,
        time_min = 60, partition = 'high2'
    params:
        out_dir = lambda wildcards: "Results/Kallisto_stringtie/{}/".format(wildcards.file),
        partition = "bmm"
    shell:
     """
     kallisto quant -i {input.idx} -o {params.out_dir} --bias -t {resources.cpus} --gtf {input.gff} {input.r1} {input.r2}
     cp {params.out_dir}abundance.tsv Results/Kallisto_stringtie/Abundance/{wildcards.file}.abundance.tsv
     """

rule salmon_idx:
    input: tx_fa = "Results/Merged/iterative/faang_refseq_ensembl.annotated.fa"
    output: "Results/salmon/merged/gentrome.fa.gz"
    params:
        dir = "Results/salmon/merged",
        partition = "med2"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300, partition = "med2"
    shell:
     """
     grep "^>" <(cat ~/reference/refgenie/equcab3/alias/equcab3/fasta/default/equcab3.fa) | cut -d " " -f 1 > {params.dir}/decoys.txt
     sed -i.bak -e 's/>//g' {params.dir}/decoys.txt
     cat {input.tx_fa} ~/reference/refgenie/equcab3/alias/equcab3/fasta/default/equcab3.fa > {params.dir}/gentrome.fa
     gzip {params.dir}/gentrome.fa
     salmon index -t {output} -d {params.dir}/decoys.txt -p {resources.cpus} -i {params.dir}/salmon_index --gencode -k 29
     """
     
rule salmon:
    input:
        r1 = "data/RNAseq_fq/{file}_R1.fastq.gz",
        r2 = "data/RNAseq_fq/{file}_R2.fastq.gz",
        gentrome = "Results/salmon/merged/gentrome.fa.gz"
    output: "Results/salmon/merged/{file}/quant.sf"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 20000, mem_mb_bmm = 20000,
        time_min = 60, partition = "med2"
    params: dir = "Results/salmon/merged/", partition = "med2"
    shell:
     """
     salmon quant -i {params.dir}/salmon_index -l A -1 {input.r1} -2 {input.r2} --validateMappings -o {params.dir}/{wildcards.file}/ -p {resources.cpus}
     """
