class RenameToFeatureType < ActiveRecord::Migration[6.1]
  def change
    rename_column :genomic_features, :type, :feature_type
  end
end
