rule prepFile:
    input: "Results/Collapsed/{file}.collapsed.gff"
    output: "Results/Rarefaction/counts/{file}.all.txt"
    conda: "../envs/cupcake.yaml"
    params:
        partition = 'bmm',
        out_prefix = "Results/Rarefaction/counts/{file}",
        in_prefix = "Results/Collapsed/{file}.collapsed"
    resources:
        cpus = 1, cpus_bmm = 1,
        mem_mb = 50000, mem_mb_bmm = 50000,
        time_min = 2200
    shell:
     """
     make_file_for_subsampling_from_collapsed.py -i {params.in_prefix} -o {params.out_prefix}
     """
     
rule subSampleByGene:
    input: "Results/Rarefaction/counts/{file}.all.txt"
    output: "Results/Rarefaction/{file}.{by}.txt"
    conda: "../envs/cupcake.yaml"
    params: partition = 'bmm'
    resources:
        cpus = 12, cpus_bmm = 12,
        mem_mb = 10000, mem_mb_bmm = 10000,
        time_min = 2200
    shell:
     """
     /home/pengsc/bin/miniconda3/envs/jupyter/bin/subsample.py --by {wildcards.by} --min_fl_count 2 --step 1000 {input} --ncores {resources.cpus}  > {output}
     """
