# frozen_string_literal: true

class VariantSite < ApplicationRecord
  def self.parse(variants_tsv)
    raise "Unable to find counts file #{variants_tsv}" unless File.exist?(variants_tsv)

    variants = []

    # ActiveRecord::Base.logger.debug "Processing variants_tsv #{variants_tsv}"
    File.readlines(variants_tsv).each do |line|
      parsed_fields = line.chomp.split("\t")
      (strain, ref_pos, variant_type, variant) = parsed_fields
      next unless strain && ref_pos && variant_type && variant

      fasta_record_id = FastaRecord.existing_fasta_strain_ids[strain]
      raise "Failed to find fasta record for stain #{strain}" unless fasta_record_id

      ref_end = Integer(ref_pos) + Integer(variant.length.to_s)
      variant = "#{variant.length}-" if variant_type.include? 'D'
      variant = "#{variant.length}N" if variant.include? 'N'
      variants << VariantSite.new(ref_start: ref_pos, ref_end: ref_end, variant_type: variant_type,
                                  variant: variant, fasta_record_id: fasta_record_id)
    end

    raise "Unable to parse any records from #{variants_tsv}" if variants.empty?

    variants
  end
end
