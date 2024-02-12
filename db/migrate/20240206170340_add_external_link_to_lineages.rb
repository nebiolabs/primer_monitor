class AddExternalLinkToLineages < ActiveRecord::Migration[6.1]
  def change
    add_column :lineages, :external_link, :string, null: true
  end
end
