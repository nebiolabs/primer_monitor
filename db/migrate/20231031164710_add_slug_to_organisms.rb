class AddSlugToOrganisms < ActiveRecord::Migration[6.1]
  def change
    add_column :organisms, :slug, :string
  end
end
