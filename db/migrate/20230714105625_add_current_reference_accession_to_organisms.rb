class AddCurrentReferenceAccessionToOrganisms < ActiveRecord::Migration[6.1]
  def change
    add_column :organisms, 'reference_accession', :string
  end
end
