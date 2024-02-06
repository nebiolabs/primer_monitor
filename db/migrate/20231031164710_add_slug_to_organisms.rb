class AddSlugToOrganisms < ActiveRecord::Migration[6.1]
  def change
    # a human-readable identifier for an organism
    add_column :organisms, :slug, :string
  end
end
