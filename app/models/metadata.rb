class FastaRecord < ApplicationRecord
  
    def self.parse(metadata_tsv)
      raise "Unable to find counts file #{metadata_tsv}" unless File.exist?(metadata_tsv)
  
      metadata = []
      record_count = 0
      # ActiveRecord::Base.logger.debug "Processing nonconverted_read_file #{metadata_tsv}"
      File.readlines(metadata_tsv).each do |line|
        parsed_fields = line.chomp.split("\t")
        (strain, _virus, gisaid_epi_isl, genbank_accession, _date, region, country, division, _location, _region_exposure, _country_exposure, _division_exposure, _segment, _length, _host, _age, _sex, _pangolin_lineage, _GISAID_clade, _originating_lab, _submitting_lab, _authors, _url, _title, _paper_url, date_submitted) = parsed_fields
        if strain && gisaid_epi_isl && genbank_accession && region && country && division && date_submitted
          record_count += 1
          if not FastaRecord.exists?(strain: strain)
            metadata << FastaRecord.new(strain: strain, gisaid_epi_isl: gisaid_epi_isl, genbank_accession: genbank_accession, region: region, country: country, division: division, date_submitted: date_submitted)
          end
        else
          next
          # ActiveRecord::Base.logger.debug "#{line} was not in the expected format\n"
        end
      end
  
      if record_count == 0
        raise "Unable to parse any records from #{metadata_tsv}"
      end
  
      metadata
    end
  end
  