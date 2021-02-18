class AddStatusToPrimerSet < ActiveRecord::Migration[6.1]
  def up
    add_column :primer_sets, :status, :string, default: 'complete'
    execute <<-SQL
      ALTER TABLE primer_sets ALTER COLUMN status DROP DEFAULT;
      ALTER TABLE primer_sets ALTER COLUMN status SET DATA TYPE VARCHAR(255);
      DROP TYPE IF EXISTS primer_set_status;
      CREATE TYPE primer_set_status AS ENUM ('pending', 'complete', 'failed'); 
      ALTER TABLE primer_sets ALTER COLUMN status SET DATA TYPE primer_set_status USING (status::primer_set_status);
      ALTER TABLE primer_sets ALTER COLUMN status SET DEFAULT 'pending'::primer_set_status;
    SQL
  end
  def down
    execute <<-SQL
      ALTER TABLE primer_sets ALTER COLUMN status DROP DEFAULT;
      ALTER TABLE primer_sets ALTER COLUMN status SET DATA TYPE VARCHAR(255);
      DROP TYPE IF EXISTS fa_tools_status;
    SQL
    remove_column :primer_sets, :status

  end


end
