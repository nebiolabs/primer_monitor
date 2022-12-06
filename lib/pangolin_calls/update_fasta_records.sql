CREATE TEMPORARY TABLE pangolin_tmp (strain text, pangolin_lineage text, conflict real, ambiguity_score text, scorpio_call text, scorpio_support real, scorpio_conflict real, scorpio_notes text, version text, pangolin_version text, scorpio_version text, constellation_version text, is_designated text, qc_status text, qc_notes text, note text);
COPY pangolin_tmp FROM :pangolin_csv_path WITH (FORMAT csv, HEADER on);
UPDATE fasta_records SET pangolin_lineage=pangolin_tmp.pangolin_lineage FROM pangolin_tmp WHERE pangolin_tmp.strain=fasta_records.strain;
DROP TABLE pangolin_tmp;