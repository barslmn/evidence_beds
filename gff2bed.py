import pandas as pd


def gff2bed(gff_file, output_file, exon_start_minus, exon_end_plus):
    # read
    gff = pd.read_table(gff_file,  # .gz
                        comment='#',
                        compression='gzip',
                        # column names
                        names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    # only exons
    gff = gff[gff['feature'] == 'gene']
    # only CHR START END cols
    gff['attribute_'] = gff.attribute.str.split('gene_name=').str[1]
    gff['gene_symbol'] = gff.attribute_.str.split(';').str[0]
    gff = gff[['seqname', 'start', 'end', 'gene_symbol']]
    # add splice sites
    gff['start'] = gff['start'] + exon_start_minus
    gff['end'] = gff['end'] + exon_end_plus
    gff.to_csv(output_file, sep='\t', index=False)


# example
gff2bed('gencode.v31lift37.annotation.gff3.gz', 'AllGenes.bed', 0, 0)
