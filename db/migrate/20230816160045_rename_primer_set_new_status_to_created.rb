class RenamePrimerSetNewStatusToCreated < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      ALTER TYPE primer_set_status RENAME VALUE 'new' TO 'created';
    SQL
  end

  def down
    execute <<-SQL
      ALTER TYPE primer_set_status RENAME VALUE 'created' TO 'new';
    SQL
  end
end
