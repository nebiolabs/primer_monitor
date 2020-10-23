class CreateFastaRecord < ActiveRecord::Migration[6.0]
  def change
    create_table :fasta_records do |t|
      t.string :strain
      t.string :genbank_accession
      t.string :gisaid_epi_isl
      t.string :region
      t.string :country
      t.string :division
      t.date :date_submitted

      t.timestamps
      end
  end
end
