import pandas as pd

count_file = snakemake.input['chained_count']
group_file = snakemake.input['collapsed_group']
out_file = snakemake.output['collapsed_count']

chained_count = pd.read_csv(count_file, sep = '\t')
collapsed_group = pd.read_csv(group_file, sep = '\t', names = ['groupID', 'indID'])
mapping = []
for i in collapsed_group.apply(lambda row: [row['groupID'], row['indID'].split(',')], axis = 1):
    mapping.extend([(i[0],n) for n in i[1]])
mapping = pd.DataFrame(mapping, columns = ['groupID', 'indID'])
chained_count = chained_count.merge(mapping, left_on = 'superPBID', right_on = 'indID')
chained_count = chained_count.groupby('groupID').sum().astype(int)
chained_count.index.names = ['superPBID']
chained_count = chained_count.replace(0, 'NA')
chained_count.to_csv(out_file, sep = '\t')
