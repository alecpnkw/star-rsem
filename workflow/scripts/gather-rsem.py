import pandas as pd

counts = pd.DataFrame({
    'sample': pd.Series(dtype='str'),
    snakemake.params['key']: pd.Series(dtype='str'),
    snakemake.params['value']: pd.Series(dtype='float'),
})

for (i,fhandle) in enumerate(snakemake.input):
    new_counts = pd.read_table(fhandle, sep = '\t')
    new_counts.loc[:, 'sample'] = snakemake.params['samples'][i]
    new_counts = new_counts[:, counts.colums]
    counts = pd.concat([counts, new_counts])

#pivot
counts_mat = counts.pivot(index = snakemake.params['key'], columns = 'sample', values = snakemake.params['value'])
counts_mat.to_csv(snakemake.output[0], sep = '\t')