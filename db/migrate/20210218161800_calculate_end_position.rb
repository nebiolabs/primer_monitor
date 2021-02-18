class CalculateEndPosition < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      UPDATE variant_sites
      SET ref_end = ref_start + SPLIT_PART(variant, 'N', 1)::int
      WHERE variant LIKE '%N%';

      UPDATE variant_sites
      SET ref_end = ref_start + SPLIT_PART(variant, '-', 1)::int
      WHERE variant LIKE '%-%';

      UPDATE variant_sites
      SET ref_end = ref_start + LENGTH(variant)
      WHERE variant NOT LIKE '%N%' AND variant NOT LIKE '%-%';
    
    SQL
      
  end

  def down
    puts "Irreversible :("
  end
end
