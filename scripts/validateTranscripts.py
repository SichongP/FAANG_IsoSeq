import pandas as pd
from Bio import SeqIO
import numpy as np

input_gff = snakemake.input['gff']
input_count = snakemake.input['counts']
input_fq = snakemake.input['fq']

output_fa = snakemake.output['fa']
output_gff = snakemake.output['gff']


fl_count = pd.read_csv(input_count, sep = '\t', index_col = 0)
fl_count = pd.melt(fl_count.reset_index(), id_vars = 'superPBID', value_name = 'fl_count', var_name = 'name').fillna(0)
fl_count[['sample', 'tissue']] = fl_count['name'].str.split('_', expand = True)
detected_in = fl_count.groupby(['superPBID', 'tissue']).apply(lambda x: (x['fl_count']>0).sum())
IDs = list(detected_in.reset_index().groupby('superPBID').sum()[(detected_in.reset_index().groupby('superPBID').sum()>1).values].index)

with open(output_fa, 'w') as out:
    for record in SeqIO.parse(input_fq, "fastq"):
        if record.id.split('|')[0] in IDs:
            SeqIO.write(record, out, "fasta")
            

merged_gff = pd.read_csv(input_gff, sep = '\t', names = 
                        ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
confirmed_gff = merged_gff[merged_gff['attribute'].str.split(';', expand = True)[0].str.split('"', expand = True)[1].isin(IDs)]
confirmed_gff = confirmed_gff[~confirmed_gff['chr'].str.startswith('chrUn')]
confirmed_gff.to_csv(output_gff, index = False, header = False, sep = '\t', doublequote = False, quotechar = "'")


