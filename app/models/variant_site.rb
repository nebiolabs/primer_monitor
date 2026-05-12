# frozen_string_literal: true

class VariantSite < ApplicationRecord
  self.ignored_columns = %w[usable_del_or_snp usable_insertion]
  def self.parse(variants_tsv, taxon)
    raise "Unable to find counts file #{variants_tsv}" unless File.exist?(variants_tsv)

    variants = []
    ref_seq = load_reference_sequence(taxon)

    # ActiveRecord::Base.logger.debug "Processing variants_tsv #{variants_tsv}"
    File.readlines(variants_tsv).each do |line|
      next if line.blank?

      accession, ref_pos, variant_type, variant = line.chomp.split("\t")
      next if variant.nil? # big insertions can push a variant off the end of the genome e.g. USA/VA-SU-SC_65/2021

      variants << (build_variant_site accession, ref_pos, variant_type, variant, taxon, ref_seq)
    end

    raise "Unable to parse any records from #{variants_tsv}" if variants.empty?

    variants
  end

  def self.build_variant_site(accession, ref_pos, variant_type, variant, taxon, ref_seq = nil)
    fasta_record_id = FastaRecord.existing_fasta_accession_ids[accession]
    raise "Failed to find fasta record for accession: \"#{accession}\"" unless fasta_record_id

    ref_pos = Integer(ref_pos) - 1 # convert 1-based to 0-based

    ref_allele = lookup_ref_allele(ref_seq, ref_pos, variant_type, variant)

    ref_end = ref_pos + variant.length
    variant = "#{variant.length}-" if variant_type.include? 'D'
    variant = "#{variant.length}N" if variant.include? 'N'
    VariantSite.new(ref_start: ref_pos.to_s, ref_end:, variant_type:,
                    variant:, ref: ref_allele, fasta_record_id:, organism_taxon_id: taxon.id)
  end

  def self.load_reference_sequence(taxon)
    fasta_path = Rails.root.join('igvstatic', taxon.organism.slug, 'ref', "#{taxon.reference_accession}.fasta")
    return nil unless File.exist?(fasta_path)

    seq = +''
    File.foreach(fasta_path) do |line|
      next if line.start_with?('>')

      seq << line.chomp
    end
    seq
  end

  def self.lookup_ref_allele(ref_seq, ref_pos, variant_type, variant)
    return nil if ref_seq.nil?

    case variant_type
    when 'X' # SNP: ref base at the position
      ref_seq[ref_pos]
    when 'I' # insertion: anchor base just before the insertion point
      ref_pos > 0 ? ref_seq[ref_pos - 1] : nil
    when 'D' # deletion: the deleted reference bases (before normalization, variant is "---")
      ref_seq[ref_pos, variant.length]
    end
  end
end
