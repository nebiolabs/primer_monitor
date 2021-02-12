class UpdateVariantsWithDashes < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      UPDATE variant_sites 
      SET variant = replace(variant, 'N', '-')
      WHERE variant_type = 'D' AND variant LIKE '%N%';
      
      
      UPDATE variant_sites 
      SET variant = concat(LENGTH(variant), '-')
      WHERE variant_type = 'D' AND NOT variant LIKE '%-%';
    SQL
  end
  def down
    puts "Irreversible :("
  end
end
