class AddDefaultAmpData < ActiveRecord::Migration[6.1]
  def up
    qpcr = AmplificationMethod.find_or_create_by(name: 'qPCR', description_url: 'https://www.neb.com/applications/dna-amplification-pcr-and-qpcr/qpcr-and-rt-qpcr' )
    lamp = AmplificationMethod.find_or_create_by(name: 'LAMP', description_url: 'https://www.neb.com/applications/dna-amplification-pcr-and-qpcr/isothermal-amplification/loop-mediated-isothermal-amplification-lamp' )

    lamp_sets = PrimerSet.where("name ilike '%lamp%'")
    lamp_sets.each do |ps|
      ps.amplification_method_id = lamp.id
      ps.save!
    end

    qpcr_sets = PrimerSet.where("not name ilike '%lamp%'")
    qpcr_sets.each do |ps|
      ps.amplification_method_id = qpcr.id
      ps.save!
    end
  end

  def down
    PrimerSet.all.each do |ps|
      ps.amplification_method_id = nil
      ps.save!
    end
    AmplificationMethod.delete_all()
  end
end