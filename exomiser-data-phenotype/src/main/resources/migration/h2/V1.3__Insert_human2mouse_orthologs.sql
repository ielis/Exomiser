INSERT INTO human2mouse_orthologs SELECT *
                                  FROM CSVREAD('${import.path}/human2mouseOrthologs.pg',
                                               'mgi_gene_id|mgi_gene_symbol|human_gene_symbol|entrez_id',
                                               'charset=UTF-8 fieldDelimiter='' fieldSeparator=| nullString=NULL');