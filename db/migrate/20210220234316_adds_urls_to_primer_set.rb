class AddsUrlsToPrimerSet < ActiveRecord::Migration[6.1]
  def change
    add_column :primer_sets, :citation_url, :string
    add_column :primer_sets, :doi, :string
    add_column :primer_sets, :amplification_method_id, :bigint
    add_foreign_key :primer_sets, :amplification_methods

  end
end
