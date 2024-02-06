class AddUniqueIndexToOrganismSlug < ActiveRecord::Migration[6.1]
  def change
    add_index :organisms, :slug, unique: true
  end
end
