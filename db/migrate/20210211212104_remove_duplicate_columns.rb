class RemoveDuplicateColumns < ActiveRecord::Migration[6.1]
  def change
    remove_column :oligos, :start_pos, :integer
    remove_column :oligos, :end_pos, :integer
  end
end
