

rule salmon_sqanti_idx:
    input: tx_fa = "Results/SQANTI3/all_samples_corrected.fasta"
    output: fa = "Results/salmon/sqanti/gentrome.fa.gz", cluster = "Results/salmon/sqanti/salmon_index/duplicate_clusters.tsv"
    params:
        dir = "Results/salmon/sqanti",
        partition = "bmm"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300
    shell:
     """
     grep "^>" <(cat {config[ref_ec3]}) | cut -d " " -f 1 > {params.dir}/decoys.txt
     sed -i.bak -e 's/>//g' {params.dir}/decoys.txt
     cat {input.tx_fa} {config[ref_ec3]} > {params.dir}/gentrome.fa
     gzip {params.dir}/gentrome.fa
     salmon index -t {output.fa} -d {params.dir}/decoys.txt -p {resources.cpus} -i {params.dir}/salmon_index --gencode -k 29
     """
     
rule salmon_sqanti:
    input:
        r1 = "data/RNAseq_fq/{file}_R1.fastq.gz",
        r2 = "data/RNAseq_fq/{file}_R2.fastq.gz",
        gentrome = "Results/salmon/sqanti/gentrome.fa.gz"
    output: "Results/salmon/sqanti/{file}/quant.sf"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 60
    params: dir = "Results/salmon/sqanti", partition = "bmm"
    shell:
     """
     salmon quant -i {params.dir}/salmon_index -l A -1 {input.r1} -2 {input.r2} --validateMappings -o {params.dir}/{wildcards.file}/ -p {resources.cpus}
     """

rule salmon_chained_idx:
    input: tx_fa = "Results/Chained/all_samples.rep.fastq"
    output: "Results/salmon/chained/gentrome.fa.gz"
    params:
        dir = "Results/salmon/chained",
        partition = "bmm"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300
    shell:
     """
     sed -n '1~4s/^@/>/p;2~4p' {input} > {params.dir}/temp_tx.fa
     grep "^>" <(cat {config[ref_ec3]}) | cut -d " " -f 1 > {params.dir}/decoys.txt
     sed -i.bak -e 's/>//g' {params.dir}/decoys.txt
     cat {params.dir}/temp_tx.fa {config[ref_ec3]} > {params.dir}/gentrome.fa
     gzip {params.dir}/gentrome.fa
     salmon index -t {output} -d {params.dir}/decoys.txt -p {resources.cpus} -i {params.dir}/salmon_index --gencode -k 29
     """
     
rule salmon_chained:
    input:
        r1 = "data/RNAseq_fq/{file}_R1.fastq.gz",
        r2 = "data/RNAseq_fq/{file}_R2.fastq.gz",
        gentrome = "Results/salmon/chained/gentrome.fa.gz"
    output: "Results/salmon/chained/{file}/quant.sf"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300
    params: dir = "Results/salmon/chained", partition = "bmm"
    shell:
     """
     salmon quant -i {params.dir}/salmon_index -l A -1 {input.r1} -2 {input.r2} --validateMappings -o {params.dir}/{wildcards.file}/ -p {resources.cpus}
     """

rule salmon_collapsed_idx:
    input: tx_fa = "Results/Collapsed/{file}/collapsed.min_fl_2.filtered.rep.fastq"
    output: "Results/salmon/collapsed/{file}/gentrome.fa.gz"
    params:
        dir = lambda wildcards: "Results/salmon/collapsed/{file}".format(file = wildcards.file),
        partition = "bmm"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300
    shell:
     """
     sed -n '1~4s/^@/>/p;2~4p' {input} > {params.dir}/temp_tx.fa
     grep "^>" <(cat {config[ref_ec3]}) | cut -d " " -f 1 > {params.dir}/decoys.txt
     sed -i.bak -e 's/>//g' {params.dir}/decoys.txt
     cat {params.dir}/temp_tx.fa {config[ref_ec3]} > {params.dir}/gentrome.fa
     gzip {params.dir}/gentrome.fa
     salmon index -t {output} -d {params.dir}/decoys.txt -p {resources.cpus} -i {params.dir}/salmon_index --gencode -k 29
     """
     
rule salmon_collapsed:
    input:
        r1 = "data/RNAseq_fq/{file}_R1.fastq.gz",
        r2 = "data/RNAseq_fq/{file}_R2.fastq.gz",
        gentrome = "Results/salmon/collapsed/{file}/gentrome.fa.gz"
    output: "Results/salmon/collapsed/{file}/quant.sf"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300
    params:
        dir = lambda wildcards: "Results/salmon/collapsed/{file}".format(file = wildcards.file),
        partition = "bmm"
    shell:
     """
     salmon quant -i {params.dir}/salmon_index -l A -1 {input.r1} -2 {input.r2} --validateMappings -o {params.dir}/ -p {resources.cpus}
     """

rule get_fasta:
    input: "Results/Annotated/all_samples.gff3"
    output: "Results/Annotated/all_samples.fa"
    conda: "../envs/gffread.yaml"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 4000, mem_mb_bmm = 4000,
        time_min = 60
    params: partition = "high2"
    shell:
     """
     awk '{{if ($3=="exon") print $0}}' {input} | gffread - -g /home/pengsc/reference/equcab3/UCSC/equCab3.fa -w {output}
     """

rule salmon_annotated_idx:
    input: tx_fa = "Results/Annotated/all_samples.fa"
    output: "Results/salmon/annotated/gentrome.fa.gz"
    params:
        dir = "Results/salmon/annotated",
        partition = "bmm"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300
    shell:
     """
     grep "^>" <(cat {config[ref_ec3]}) | cut -d " " -f 1 > {params.dir}/decoys.txt
     sed -i.bak -e 's/>//g' {params.dir}/decoys.txt
     cat {input.tx_fa} {config[ref_ec3]} > {params.dir}/gentrome.fa
     gzip {params.dir}/gentrome.fa
     salmon index -t {output} -d {params.dir}/decoys.txt -p {resources.cpus} -i {params.dir}/salmon_index --gencode -k 29
     """
     
rule salmon_annotated:
    input:
        r1 = "data/RNAseq_fq/{file}_R1.fastq.gz",
        r2 = "data/RNAseq_fq/{file}_R2.fastq.gz",
        gentrome = "Results/salmon/annotated/gentrome.fa.gz"
    output: "Results/salmon/annotated/{file}/quant.sf"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 60
    params: dir = "Results/salmon/annotated/", partition = "bmm"
    shell:
     """
     salmon quant -i {params.dir}/salmon_index -l A -1 {input.r1} -2 {input.r2} --validateMappings -o {params.dir}/{wildcards.file}/ -p {resources.cpus}
     """
     
rule salmon_refseq_idx:
    input: tx_fa = "/home/pengsc/reference/equcab3/RefSeq/AnnotationRelease103/GCF_002863925.1_EquCab3.0_cds_from_genomic.fna.gz"
    output: "Results/salmon/refseq/gentrome.fa.gz"
    params:
        dir = "Results/salmon/refseq",
        partition = "med2"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300
    shell:
     """
     grep "^>" <(cat {config[ref_ec3]}) | cut -d " " -f 1 > {params.dir}/decoys.txt
     sed -i.bak -e 's/>//g' {params.dir}/decoys.txt
     cat <(gunzip -c {input.tx_fa}) {config[ref_ec3]} > {params.dir}/gentrome.fa
     gzip {params.dir}/gentrome.fa
     salmon index -t {output} -d {params.dir}/decoys.txt -p {resources.cpus} -i {params.dir}/salmon_index -k 29
     """
     
rule salmon_refseq:
    input:
        r1 = "data/RNAseq_fq/{file}_R1.fastq.gz",
        r2 = "data/RNAseq_fq/{file}_R2.fastq.gz",
        gentrome = "Results/salmon/refseq/gentrome.fa.gz"
    output: "Results/salmon/refseq/{file}/quant.sf"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 60
    params: dir = "Results/salmon/refseq/", partition = "bmm"
    shell:
     """
     salmon quant -i {params.dir}/salmon_index -l A -1 {input.r1} -2 {input.r2} --validateMappings -o {params.dir}/{wildcards.file}/ -p {resources.cpus}
     """
     
rule salmon_ensembl_idx:
    input: tx_fa = "/home/pengsc/reference/equcab3/ensemble/cdna/Equus_caballus.EquCab3.0.cdna.all.fa.gz"
    output: "Results/salmon/ensembl/gentrome.fa.gz"
    params:
        dir = "Results/salmon/ensembl",
        partition = "med2"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 300
    shell:
     """
     grep "^>" <(cat {config[ref_ec3]}) | cut -d " " -f 1 > {params.dir}/decoys.txt
     sed -i.bak -e 's/>//g' {params.dir}/decoys.txt
     cat <(gunzip -c {input.tx_fa}) {config[ref_ec3]} > {params.dir}/gentrome.fa
     gzip {params.dir}/gentrome.fa
     salmon index -t {output} -d {params.dir}/decoys.txt -p {resources.cpus} -i {params.dir}/salmon_index --gencode -k 29
     """
     
rule salmon_ensembl:
    input:
        r1 = "data/RNAseq_fq/{file}_R1.fastq.gz",
        r2 = "data/RNAseq_fq/{file}_R2.fastq.gz",
        gentrome = "Results/salmon/ensembl/gentrome.fa.gz"
    output: "Results/salmon/ensembl/{file}/quant.sf"
    conda: "../envs/salmon.yaml"
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 30000, mem_mb_bmm = 30000,
        time_min = 60
    params: dir = "Results/salmon/ensembl/", partition = "bmm"
    shell:
     """
     salmon quant -i {params.dir}/salmon_index -l A -1 {input.r1} -2 {input.r2} --validateMappings -o {params.dir}/{wildcards.file}/ -p {resources.cpus}
     """
