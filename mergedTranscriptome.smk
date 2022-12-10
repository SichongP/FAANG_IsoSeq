import glob, os
configfile: "config.yaml"
wildcard_constraints:
    file="[a-zA-Z0-9_]+"
    
ALL_RNASEQ = []
for file in glob.glob("data/RNAseq_fq/*R1*gz"):
    file = '_'.join(os.path.split(file)[1].split('.')[0].split('_')[:2])
    ALL_RNASEQ.append(file)
    
rule all:
    input: expand("Results/Merged/iterative/mapping/stats/{file}.stats.txt", file = ALL_RNASEQ)

rule bam_stat:
    input: bam = "Results/Merged/iterative/mapping/{file}.sorted.bam", bai = "Results/Merged/iterative/mapping/{file}.sorted.bam.bai"
    output: "Results/Merged/iterative/mapping/stats/{file}.stats.txt"
    params: partition='low2'
    resources: mem_mb=8000, time_min=120, cpus=4, partition='low2'
    conda: "../envs/samtools.yaml"
    shell:
     """
     samtools stats -@ {resources.cpus} {input.bam} 1>{output}
     """

rule indexIsoseqTranscriptome:
    input: fa = "Results/Merged/iterative/faang_refseq_ensembl.annotated.fa"
    output: "Results/Merged/iterative/faang_refseq_ensembl.annotated.fa.amb"
    conda: "../envs/bwa.yaml"
    params: partition='med2'
    resources: mem_mb=4000, mem_mb_med2=8000, time_min=100, cpus=2, cpus_med2=2, partition='med2'
    shell:
     """
     bwa index {input.fa}
     """

rule IsoseqTranscriptome:
    input: fq = expand("Results/trimmed/{{file}}_{read}.fq.gz", read = ["val_1", "val_2"]), ref = "Results/Merged/iterative/faang_refseq_ensembl.annotated.fa", idx = "Results/Merged/iterative/faang_refseq_ensembl.annotated.fa.amb"
    output: temp("Results/Merged/iterative/mapping/{file}.sam")
    params: partition='med2'
    resources: mem_mb=16000, mem_mb_bmm=16000, time_min=900, cpus=8, cpus_bmm=8, partition='med2'
    conda: "../envs/bwa.yaml"
    shell:
     """
     bwa mem -t {resources.cpus} -R "@RG\\tID:{wildcards.file}SM:{wildcards.file}PL:illumina" -V -o {output} {input.ref} {input.fq} 
     """

rule TranscriptomeToBAM:
    input: "Results/Merged/iterative/mapping/{file}.sam"
    output: bam = "Results/Merged/iterative/mapping/{file}.sorted.bam", bai = "Results/Merged/iterative/mapping/{file}.sorted.bam.bai"
    params: partition='med2'
    resources: mem_mb=8000, mem_mb_med2=8000, time_min=900, cpus=4, cpus_med2=4, partition='med2'
    conda: "../envs/sambamba.yaml"
    shell:
     """
     sambamba view -S -o {output.bam}.temp.bam -t {resources.cpus} -f bam {input}
     sambamba sort -m 7G -o {output.bam} -t {resources.cpus} {output.bam}.temp.bam
     rm {output.bam}.temp.bam
     sambamba index -t {resources.cpus} {output.bam}
     """