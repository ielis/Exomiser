INSERT INTO disease_hp SELECT *
                       FROM CSVREAD('${import.path}/diseaseHp.pg', 'disease_id|hp_id',
                                    'charset=UTF-8 fieldDelimiter='' fieldSeparator=| nullString=NULL');