class Variants < ApplicationRecord
  
    def self.parse(variants_tsv)
      raise "Unable to find counts file #{variants_tsv}" unless File.exist?(variants_tsv)
  
      variants = []
  
      ActiveRecord::Base.logger.debug "Processing variants_tsv #{variants_tsv}"
      File.readlines(variants_tsv).each do |line|
        parsed_fields = line.chomp.split("\t")
        (strain, ref_pos, variant_type, variant) = parsed_fields
        if strain && ref_pos && variant_type && variant
          variants << variants.new(strain: strain, position: ref_pos, variant_type: variant_type, variant: variant)
        else
          ActiveRecord::Base.logger.debug "#{line} was not in the expected format\n"
        end
      end
  
      if variants.empty?
        raise "Unable to parse any records from #{variants_tsv}"
      end
  
      variants
    end
  end
  