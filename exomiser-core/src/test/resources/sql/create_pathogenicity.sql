DROP TABLE pathogenicity IF EXISTS;

CREATE TABLE pathogenicity (
    chromosome smallint NOT NULL,
    "position" integer NOT NULL,
    ref character(1) NOT NULL,
    alt character(1) NOT NULL,
    sift double precision,
    polyphen double precision,
    mut_taster double precision,
    cadd double precision,
    cadd_raw double precision 
);
