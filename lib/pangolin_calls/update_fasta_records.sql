CREATE TEMPORARY TABLE pangolin_tmp (taxon varchar NOT NULL, lineage varchar NOT NULL, conflict varchar, ambiguity_score float4, scorpio_call varchar, scorpio_support numeric, scorpio_conflict numeric, scorpio_notes text, version varchar NOT NULL, pangolin_version date NOT NULL, scorpio_version varchar NOT NULL, constellation_version varchar NOT NULL, is_designated text NOT NULL, qc_status varchar, qc_notes varchar, note varchar); 

COPY pangolin_tmp FROM :pangolin_csv_path WITH (format csv, header on);
UPDATE fasta_records SET pangolin_lineage=pangolin_tmp.pangolin_lineage FROM pangolin_tmp WHERE pangolin_tmp.strain=fasta_records.strain;
DROP TABLE pangolin_tmp;