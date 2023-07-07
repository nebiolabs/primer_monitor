class MakePangolinLineageForeignKey < ActiveRecord::Migration[6.1]
  def change
    rename_column :pangolin_calls, :lineage, '_lineage_name'
    change_column_null :pangolin_calls, :_lineage_name, true
    add_reference :pangolin_calls, :lineage
  end
end
