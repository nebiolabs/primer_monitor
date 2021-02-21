class CreateAmplificationMethods < ActiveRecord::Migration[6.1]
  def change
    create_table :amplification_methods do |t|
      t.string :name
      t.string :description_url

      t.timestamps
    end
  end
end
