class CorrectUniqueConstraintOnVariantSites < ActiveRecord::Migration[6.1]
  def self.up
    execute "ALTER TABLE variant_sites DROP CONSTRAINT ensure_variant_unique;
                 ALTER TABLE variant_sites ADD CONSTRAINT ensure_variant_unique UNIQUE(ref_start, fasta_record_id, variant_type);"
  end

  def self.down
    execute "ALTER TABLE variant_sites DROP CONSTRAINT ensure_variant_unique;
                 ALTER TABLE variant_sites ADD CONSTRAINT ensure_variant_unique UNIQUE(ref_start, fasta_record_id);"
  end
end