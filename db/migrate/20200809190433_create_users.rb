class CreateUsers < ActiveRecord::Migration[6.0]
  def change
    create_table :users do |t|
      t.string :first,  null: false
      t.string :last, null: false
      t.string :email, null: false, unique: true
      t.boolean :activated, default: false, null: false

      t.timestamps
    end
  end
end
