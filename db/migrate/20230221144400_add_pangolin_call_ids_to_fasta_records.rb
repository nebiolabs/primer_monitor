class AddPangolinCallIdsToFastaRecords < ActiveRecord::Migration[6.1]
  def change
    add_reference :fasta_records, :pangolin_call_id, null: true, foreign_key: {to_table: :pangolin_calls}
    add_reference :fasta_records, :pending_pangolin_call_id, null: true, foreign_key: {to_table: :pangolin_calls}
  end
end
