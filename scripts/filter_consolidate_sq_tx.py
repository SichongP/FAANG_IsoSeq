import pandas as pd
import numpy as np
import mygene, glob

input_class = snakemake.input['classification']
input_sq_gtf = snakemake.input['sq_gtf']
input_junctions = snakemake.input['junctions']
salmon_cluster = snakemake.input['salmon_cluster']
salmon_counts = snakemake.input['salmon_counts']

output_gff = snakemake.output['gff']


classification = pd.read_csv(input_class, sep = '\t')
classification['associated_transcript'] = classification['associated_transcript'].str.split(":",expand = True)[1].fillna('novel')

mg = mygene.MyGeneInfo()
gene_names = mg.querymany(classification['associated_gene'].unique(), fields = ['name','symbol'], species = 9796, as_dataframe = True, verbose = False, df_index = False)
gene_names = gene_names.groupby('query').apply(lambda x: pd.Series([','.join(x['name'].dropna().to_list()), ','.join(x['symbol'].dropna().to_list())])).replace('', np.nan).reset_index()
gene_names.columns = ['query', 'name', 'symbol']
classification = classification.merge(gene_names.drop_duplicates().reset_index(), left_on = 'associated_gene', right_on = 'query', how = 'left', validate = 'm:1')

ref_annot = pd.read_csv("data/ensembl.ref.gtf", sep = '\t', names = ['seqID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
attributes = ref_annot['attributes'].str.split(';', expand = True).iloc[:,:2]
ref_annot['transcript_id'] = attributes[0].str.split('"', expand = True)[1].str.split(":", expand = True)[1]
ref_annot['gene_id'] = attributes[1].str.split('"', expand = True)[1]

sq_gtf = pd.read_csv(input_sq_gtf, sep = '\t', names = ['seqID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
attributes = sq_gtf['attributes'].str.split(';', expand = True).iloc[:,:2]
sq_gtf['PBtranscript_id'] = attributes[0].str.split('"', expand = True)[1]
sq_gtf['PBgene_id'] = attributes[1].str.split('"', expand = True)[1]

read_count_df = pd.DataFrame(index = ['Name'])
tpm_df = pd.DataFrame(columns = ['Name', 'TPM'])
for file in salmon_counts:
    sample, tissue = file.split("/")[3].split('_')
    read_count = pd.read_csv(file, sep = '\t', index_col = 0)[['NumReads']]
#    tpm = pd.read_csv(i, sep = '\t')[['TPM']]
    read_count.columns = [sample + '_' + tissue]
#    tpm.columns = [sample]
    read_count_df = read_count.join(read_count_df, how = 'left')
read_count_df = read_count_df.reset_index()
sample_names = read_count_df.columns[1:]

salmon_cluster = pd.read_csv(salmon_cluster, sep = '\t')
has_duplicate = salmon_cluster['RetainedRef'].to_list() + salmon_cluster['DuplicateRef'].to_list()
classification['has_duplicate'] = classification['isoform'].isin(has_duplicate)
classification = classification.merge(salmon_cluster, left_on = 'isoform', right_on = 'DuplicateRef', how = 'left')
classification['salmonIsoformID'] = classification.apply(lambda x: x['isoform'] if x['RetainedRef'] is np.nan else x['RetainedRef'], axis = 1)
classification = classification.merge(read_count_df, left_on = 'salmonIsoformID', right_on = 'Name', how = 'left')
tx_with_cov = classification[((classification[sample_names]>2).sum(axis = 1)>=2) & (classification[sample_names].mean(axis = 1) > 5)]['isoform'].to_list()

junctions = pd.read_csv(input_junctions, sep = '\t')
junctions['log_mean_unique_coverage'] = junctions.loc[:,junctions.columns.str.contains('unique')].mean(axis = 1)
no_junction_cov = junctions[junctions['log_mean_unique_coverage']==0]['isoform'].unique()

#intrapriming = classification[['isoform', 'perc_A_downstream_TTS', 'seq_A_downstream_TTS']]
#intrapriming['intrapriming'] = (intrapriming['perc_A_downstream_TTS'] >= 60) | intrapriming['seq_A_downstream_TTS'].str.startswith('A'*6)
#intrapriming_tx = intrapriming[intrapriming['intrapriming']]['isoform']
NMD_tx = classification[classification['predicted_NMD'].fillna(False).astype(bool)]['isoform']
tx_to_remove = NMD_tx.to_list() + list(no_junction_cov)

classification_filtered = classification[(~classification['isoform'].isin(tx_to_remove)) & (classification['isoform'].isin(tx_with_cov))]
sq_gtf_filtered = sq_gtf[sq_gtf['PBtranscript_id'].isin(classification_filtered.isoform)]

merged = sq_gtf_filtered.merge(classification_filtered, left_on = 'PBtranscript_id', right_on = 'isoform', how = 'left')

def parse_transcript(row):
    attributes = ""
    annotated_classes = ['full-splice_match', 'incomplete-splice_match']
    tx_id, gene_id, symbol = row['associated_transcript'] if row['structural_category'] in annotated_classes else row['PBtranscript_id'], row['associated_gene'], row['symbol']
    if row['type'] == 'transcript':
        incomplete_3_prime = True if row['perc_A_downstream_TTS'] >= 60 or str(row['seq_A_downstream_TTS']).startswith('A'*6) else False
        TX = [row['seqID'], 'PacBio', 'mRNA', row['start'], row['end'], '.', row['strand_x'], '.', 
              "ID={tx};Parent={gene};Name={sym};incomplete_3_prime={inc3prime};structural_category={structural_category};num_exons={num_exons}; FAANG_ID={faang_id}".format(
                  tx = tx_id, gene = gene_id, sym = symbol,
                  inc3prime = incomplete_3_prime, num_exons = row['exons'],
                  structural_category = row['structural_category'],
                  faang_id = row['isoform']
              )
             ]
        if np.isnan(row['CDS_length']):
            CDS, UTR5, UTR3 = [], [], []
        else:
            CDS = [row['seqID'], 'PacBio', 'CDS', min(row['CDS_genomic_start'], row['CDS_genomic_end']), max(row['CDS_genomic_start'], row['CDS_genomic_end']), '.', row['strand_x'], '.', "ID={tx}_cds;Parent={gene};Name={sym}".format(tx = tx_id, gene = gene_id, sym = symbol)]
            if row['strand_x'] == '+':
                utr5_start, utr5_end, utr3_start, utr3_end = row['start'], int(row['CDS_genomic_start'] - 1), int(row['CDS_genomic_end'] + 1), row['end']
            elif row['strand_x'] == '-':
                utr5_start, utr5_end, utr3_start, utr3_end = int(row['CDS_genomic_start'] + 1), row['end'], row['start'], int(row['CDS_genomic_end'] - 1), 
            UTR5 = [row['seqID'], 'PacBio', 'five_prime_UTR', utr5_start, utr5_end, '.', row['strand_x'], '.', "Parent={tx}".format(tx = tx_id)] if utr5_start < utr5_end else None
            UTR3 = [row['seqID'], 'PacBio', 'three_prime_UTR', utr3_start, utr3_end, '.', row['strand_x'], '.', "Parent={tx}".format(tx = tx_id)] if utr3_start < utr3_end else None
        return [TX, CDS, UTR5, UTR3]
    elif row['type'] == 'exon':
        EXONS = [row['seqID'], 'PacBio', 'exon', row['start'], row['end'], '.', row['strand_x'], '.', "Parent={tx}".format(tx = tx_id)]
        return pd.Series(EXONS, index = ['seqID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    
parsed_TX = []
for i in merged[merged['type']=='transcript'].apply(parse_transcript, axis = 1):
    i = [j for j in i if j is not None]
    parsed_TX.extend(i)
parsed_TX = pd.DataFrame(parsed_TX, columns = ['seqID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
parsed_TX = pd.concat([parsed_TX, merged[merged['type']=='exon'].apply(parse_transcript, axis = 1)])
parsed_TX = parsed_TX.sort_values(['seqID', 'start', 'type', 'end']).dropna()
parsed_TX[['start', 'end']] = parsed_TX[['start', 'end']].astype(int)
parsed_TX.to_csv(output_gff, index = False, sep = '\t', header = False)
