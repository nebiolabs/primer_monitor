class RenamePrimerSetStatusEnumValues < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      ALTER TYPE primer_set_status RENAME VALUE 'pending' TO 'new';
      ALTER TYPE primer_set_status RENAME VALUE 'ready' TO 'processing';
    SQL
  end

  def down
    execute <<-SQL
      ALTER TYPE primer_set_status RENAME VALUE 'new' TO 'pending';
      ALTER TYPE primer_set_status RENAME VALUE 'processing' TO 'ready';
    SQL
  end
end
