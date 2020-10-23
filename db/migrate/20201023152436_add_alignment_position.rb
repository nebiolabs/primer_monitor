class AddAlignmentPosition < ActiveRecord::Migration[6.0]
  def change
    add_column :oligos, :start_pos, :integer
    add_column :oligos, :end_pos, :integer
  end
end
