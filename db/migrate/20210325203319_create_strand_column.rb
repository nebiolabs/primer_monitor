class CreateStrandColumn < ActiveRecord::Migration[6.1]
  def change
    add_column :oligos, :strand, :string
  end
end