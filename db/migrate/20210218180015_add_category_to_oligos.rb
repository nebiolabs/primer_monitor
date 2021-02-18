class AddCategoryToOligos < ActiveRecord::Migration[6.1]
  def up
    add_column :oligos, :locus, :string
    add_column :oligos, :category, :string
    execute <<-SQL
      ALTER TABLE oligos ALTER COLUMN category DROP DEFAULT;
      ALTER TABLE oligos ALTER COLUMN category SET DATA TYPE VARCHAR(255);
      DROP TYPE IF EXISTS oligo_category;
      CREATE TYPE oligo_category AS ENUM ('BIP', 'FIP', 'LB', 'LF', 'B3', 'F3', 'Reverse', 'Forward','Probe'); 
      ALTER TABLE oligos ALTER COLUMN category SET DATA TYPE oligo_category USING (category::oligo_category);
    SQL
  end
  def down
    execute <<-SQL
      ALTER TABLE primer_sets ALTER COLUMN category SET DATA TYPE VARCHAR(255);
      DROP TYPE IF EXISTS oligo_category;
    SQL
    remove_column :oligos, :category
    add_column :oligos, :locus

  end
end
