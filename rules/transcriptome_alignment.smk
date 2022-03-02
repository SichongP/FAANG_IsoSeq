wildcard_constraints:
    ref = "[a-zA-Z]+"

rule trim:
    input: expand("data/RNAseq_fq/{{file}}_{read}.fastq.gz", read = ["R1", "R2"])
    output: expand("Results/trimmed/{{file}}_{read}.fq.gz", read = ["val_1", "val_2"])
    conda: "../envs/trim_galore.yaml"
    resources: time_min=8000, mem_mb=8000
    params: partition='bmm'
    shell:
     """
     trim_galore --paired -o Results/trimmed/  --basename {wildcards.file} {input}
     """
     
rule mappingEnsembleTranscriptome:
    input: fq = expand("Results/trimmed/{{file}}_{read}.fq.gz", read = ["val_1", "val_2"]), ref = "/home/pengsc/reference/equcab3/ensemble/cdna/Equus_caballus.EquCab3.0.cdna.all.fa.gz"
    output: temp("Results/mapping/ensembleTranscriptomeOnly/{file}.sam")
    params: partition='bmm'
    resources: mem_mb=20000, mem_mb_bmm=20000, time_min=900, cpus=8, cpus_bmm=8
    conda: "../envs/bwa.yaml"
    shell:
     """
     bwa mem -t {resources.cpus} -R "@RG\\tID:{wildcards.file}SM:{wildcards.file}PL:illumina" -V -o {output} {input.ref} {input.fq} 
     """

rule TranscriptomeToBAM:
    input: "Results/mapping/{ref}/{file}.sam"
    output: bam = "Results/mapping/{ref}/{file}.sorted.bam", bai = "Results/mapping/{ref}/{file}.sorted.bam.bai"
    params: partition='bmm'
    resources: mem_mb=8000, mem_mb_bmm=8000, time_min=900, cpus=3, cpus_bmm=3
    conda: "../envs/sambamba.yaml"
    shell:
     """
     sambamba view -S -o {output.bam}.temp.bam -t {resources.cpus} -f bam {input}
     sambamba sort -m 7G -o {output.bam} -t {resources.cpus} {output.bam}.temp.bam
     rm {output.bam}.temp.bam
     sambamba index -t {resources.cpus} {output.bam}
     """
     
rule mappingRefSeqTranscriptome:
    input: fq = expand("Results/trimmed/{{file}}_{read}.fq.gz", read = ["val_1", "val_2"]), ref = "/home/pengsc/reference/equcab3/RefSeq/rna/equcab3_RefSeq_cdna_from_genomic.fna"
    output: temp("Results/mapping/refSeqTranscriptomeOnly/{file}.sam")
    params: partition='bmm'
    resources: mem_mb=20000, mem_mb_bmm=20000, time_min=900, cpus=8, cpus_bmm=8
    conda: "../envs/bwa.yaml"
    shell:
     """
     bwa mem -t {resources.cpus} -R "@RG\\tID:{wildcards.file}SM:{wildcards.file}PL:illumina" -V -o {output} {input.ref} {input.fq} 
     """

rule indexIsoseqTranscriptome:
    input: fa = "Results/Annotated/all_samples.fa"
    output: "Results/Annotated/all_samples.fa.amb"
    conda: "../envs/bwa.yaml"
    params: partition='bmm'
    resources: mem_mb=8000, mem_mb_bmm=8000, time_min=100, cpus=1, cpus_bmm=1
    shell:
     """
     bwa index all_samples.fa
     """

rule IsoseqTranscriptome:
    input: fq = expand("Results/trimmed/{{file}}_{read}.fq.gz", read = ["val_1", "val_2"]), ref = "Results/Annotated/all_samples.fa", idx = "Results/Annotated/all_samples.fa.amb"
    output: temp("Results/mapping/IsoSeqTranscriptomeOnly/{file}.sam")
    params: partition='bmm'
    resources: mem_mb=20000, mem_mb_bmm=20000, time_min=900, cpus=8, cpus_bmm=8
    conda: "../envs/bwa.yaml"
    shell:
     """
     bwa mem -t {resources.cpus} -R "@RG\\tID:{wildcards.file}SM:{wildcards.file}PL:illumina" -V -o {output} {input.ref} {input.fq} 
     """

rule bam_stat:
    input: bam = "Results/mapping/{ref}/{file}.sorted.bam", bai = "Results/mapping/{ref}/{file}.sorted.bam.bai"
    output: "Results/mapping/stats/{ref}_{file}.stats.txt"
    params: partition='med2'
    resources: mem_mb=8000, time_min=120, cpus=4
    conda: "../envs/samtools.yaml"
    shell:
     """
     samtools stats -@ {resources.cpus} {input.bam} 1>{output}
     """