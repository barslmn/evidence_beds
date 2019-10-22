import pandas as pd

df = pd.read_csv('refSeq.Exons.hg19.bed', sep='\t', names=[
                 'chr', 'start', 'end', 'NM_exon', 'something', 'orientation'])

df['NMname'], df['NM'], df['cdsname'], df['exon'], df['rest'] = df.NM_exon.str.split(
    '_', 4).str

counts = df['NM'].value_counts().to_frame()
counts = counts.reset_index()
counts.columns = ['NM', 'counts']
merged = pd.merge(df, counts, on='NM')

# ++++++++++++ #
fmerged = merged[merged['orientation'] == '+']
fe1 = fmerged[fmerged['counts'] == 1]
fmerged = fmerged[fmerged['counts'] != 1]

fex_ultima = fmerged[(fmerged['counts'] - fmerged['exon'].astype(int)) == 1]
fex_penultima = fmerged[(fmerged['counts'] - fmerged['exon'].astype(int)) == 2]
fex_penultima['start'] = fex_penultima['end'].astype('int') - 50

fmerged = pd.merge(fex_penultima[['chr', 'start', 'NM']], fex_ultima[['end', 'NM']], on='NM')

# ---------- #
rmerged = merged[merged['orientation'] == '-']
re1 = rmerged[rmerged['counts'] == 1]
rmerged = rmerged[rmerged['counts'] != 1]

rex_ultima = rmerged[rmerged['exon'] == '0']
rex_penultima = rmerged[rmerged['exon'] == '1']
rex_penultima['end'] = rex_penultima['start'].astype('int') + 50

rmerged = pd.merge(rex_ultima[['chr', 'start', 'NM']], rex_penultima[['end', 'NM']], on='NM')

# ±1±1±1±1±1±1±1±1 #
dfs = [fe1, re1, fmerged, rmerged]
dfs = [df[['chr', 'start', 'end', 'NM']] for df in dfs]
concated = pd.concat(dfs)

# filter by canonical
df_canonical = pd.read_csv('knownCanonical.bed', sep='\t', names=['0', '1', '2', '3', '4', '5', 'NM'])
df_canonical.dropna(inplace=True)
df_canonical['NMname'], df_canonical['NMid'] = df_canonical.NM.str.split('_', 1).str

concated['NMid'], concated['NMversion'] = concated.NM.str.split('.', 1).str
print(concated.head())
canonicals = pd.merge(concated, df_canonical[['NMid']], on='NMid')


# add anotation
canonicals['PVS1M'] = 'PVS1M50nt'
canonicals[['chr', 'start', 'end', 'PVS1M']].to_csv('PVS1M.bed', index=False, sep='\t', header=False)

# fff = pd.merge(fex_ultima, fex_penultima[['NM']], how='outer', indicator=True, on='NM')
# diff_ultima = fff[fff['_merge'] == 'left_only']
# print(diff_ultima)
# print(len(diff_ultima))
# print(diff_ultima['counts'].value_counts())
