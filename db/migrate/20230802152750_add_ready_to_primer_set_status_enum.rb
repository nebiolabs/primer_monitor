class AddReadyToPrimerSetStatusEnum < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      ALTER TYPE primer_set_status ADD VALUE 'ready';
    SQL
  end

  def down
    raise ActiveRecord::IrreversibleMigration
  end
end
