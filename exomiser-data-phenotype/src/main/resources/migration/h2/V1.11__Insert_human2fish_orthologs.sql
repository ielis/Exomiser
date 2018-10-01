INSERT INTO human2fish_orthologs SELECT *
                                 FROM CSVREAD('${import.path}/human2fishOrthologs.pg',
                                              'zfin_gene_id|zfin_gene_symbol|human_gene_symbol|entrez_id',
                                              'charset=UTF-8 fieldDelimiter='' fieldSeparator=| nullString=NULL');