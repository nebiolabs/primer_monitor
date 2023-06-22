class RenamePangolinCallFields < ActiveRecord::Migration[6.1]
  def change
    rename_column :fasta_records, :pending_pangolin_call_id_id, :pending_pangolin_call_id
    rename_column :fasta_records, :pangolin_call_id_id, :pangolin_call_id
  end
end