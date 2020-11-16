import pandas as pd

Ortholog_df = pd.read_csv(
    'https://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data',
    sep='\t',header=None)

Ortholog_df = Ortholog_df.iloc[:,[0,1,3]]
Ortholog_df.columns = ['Ortholog_id','Tax_id','Symbol']
valid_orthologs = Ortholog_df.Ortholog_id.value_counts().index[
    Ortholog_df.Ortholog_id.value_counts()>1]
Ortholog_df = Ortholog_df[Ortholog_df.Ortholog_id.isin(valid_orthologs)]

m = Ortholog_df[Ortholog_df.Tax_id == 10090]
h = Ortholog_df[Ortholog_df.Tax_id == 9606]
m2h = m.merge(h, on='Ortholog_id')
m2h = pd.DataFrame(m2h.set_index('Symbol_x')['Symbol_y'])
m2h.name = 'Mouse_symbol'
m2h.columns = ['Human_symbol']
m2h.to_csv(
    '/project/shared/xiao_wang/software/rnaseqDE//scripts/M2H_symbol_conversion.txt',
    sep='\t')