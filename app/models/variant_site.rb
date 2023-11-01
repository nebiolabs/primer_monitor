# frozen_string_literal: true

class VariantSite < ApplicationRecord
  self.ignored_columns = %w[usable_del_or_snp usable_insertion]
  def self.parse(variants_tsv, organism)
    raise "Unable to find counts file #{variants_tsv}" unless File.exist?(variants_tsv)

    variants = []

    # ActiveRecord::Base.logger.debug "Processing variants_tsv #{variants_tsv}"
    File.readlines(variants_tsv).each do |line|
      next if line.blank?

      accession, ref_pos, variant_type, variant = line.chomp.split("\t")
      next if variant.nil? # big insertions can push a variant off the end of the genome e.g. USA/VA-SU-SC_65/2021

      variants << (build_variant_site accession, ref_pos, variant_type, variant, organism)
    end

    raise "Unable to parse any records from #{variants_tsv}" if variants.empty?

    variants
  end

  def self.build_variant_site(accession, ref_pos, variant_type, variant, organism)
    fasta_record_id = FastaRecord.existing_fasta_accession_ids[accession]
    raise "Failed to find fasta record for accession: \"#{accession}\"" unless fasta_record_id

    ref_pos = Integer(ref_pos) - 1 # convert 1-based to 0-based

    ref_end = ref_pos + variant.length
    variant = "#{variant.length}-" if variant_type.include? 'D'
    variant = "#{variant.length}N" if variant.include? 'N'
    organism_id = Organism.find_by(slug: organism).id
    VariantSite.new(ref_start: ref_pos.to_s, ref_end: ref_end, variant_type: variant_type,
                    variant: variant, fasta_record_id: fasta_record_id, organism_id: organism_id)
  end
end
