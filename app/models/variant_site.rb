# frozen_string_literal: true

class VariantSite < ApplicationRecord
  self.ignored_columns = %w[usable_del_or_snp usable_insertion]
  def self.parse(variants_tsv)
    raise "Unable to find counts file #{variants_tsv}" unless File.exist?(variants_tsv)

    variants = []

    # ActiveRecord::Base.logger.debug "Processing variants_tsv #{variants_tsv}"
    File.readlines(variants_tsv).each do |line|
      next if line.blank?

      strain, ref_pos, variant_type, variant = line.chomp.split("\t")
      next if variant.nil? # big insertions can push a variant off the end of the genome e.g. USA/VA-SU-SC_65/2021

      fasta_record_id = FastaRecord.existing_fasta_strain_ids[strain]
      raise "Failed to find fasta record for strain: #{strain}" unless fasta_record_id

      ref_pos = Integer(ref_pos) - 1 # convert 1-based to 0-based

      ref_end = ref_pos + variant.length
      variant = "#{variant.length}-" if variant_type.include? 'D'
      variant = "#{variant.length}N" if variant.include? 'N'
      variants << VariantSite.new(ref_start: ref_pos.to_s, ref_end: ref_end, variant_type: variant_type,
                                  variant: variant, fasta_record_id: fasta_record_id)
    end

    raise "Unable to parse any records from #{variants_tsv}" if variants.empty?

    variants
  end
end
