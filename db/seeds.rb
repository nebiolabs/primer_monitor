# frozen_string_literal: true

# This file should contain all the record creation needed to seed the database with its default values.
# The data can then be loaded with the rails db:seed command (or created alongside the database with db:setup).

# password should be changed...
admin = User.create_with(first: 'Admin', last: 'User',
                         password: Rails.application.credentials.admin_password,
                         password_confirmation: Rails.application.credentials.admin_password,
                         active: true, approved: true, confirmed: true)
            .find_or_create_by!(email: ENV['ADMIN_EMAIL'])

%w[administrator pi operator participant].each do |role_name|
  Role.find_or_create_by!(name: role_name)
end

admin.roles << Role.find_by(name: 'administrator') if admin.roles.empty?

LineageCaller.find_or_create_by!(name: 'default', script_name: 'default')

pangolin = LineageCaller.find_or_create_by!(name: 'pangolin', script_name: 'pangolin')

sars = Organism.find_or_create_by!(name: 'SARS-CoV-2', slug: 'sars-cov-2')
OrganismTaxon.find_or_create_by!(ncbi_taxon_id: 2_697_049, organism_id: sars.id, lineage_caller_id: pangolin.id,
                                 name: 'SARS-CoV-2', reference_accession: 'NC_045512.2')

nextclade_rsva = LineageCaller.find_or_create_by!(name: 'nextclade_rsva', script_name: 'nextclade')
nextclade_rsvb = LineageCaller.find_or_create_by!(name: 'nextclade_rsvb', script_name: 'nextclade')

rsv = Organism.find_or_create_by!(name: 'RSV', slug: 'rsv')
OrganismTaxon.find_or_create_by!(ncbi_taxon_id: 208_893, organism_id: rsv.id, lineage_caller_id: nextclade_rsva.id,
                                 name: 'RSV-A', reference_accession: 'NC_038235.1')
OrganismTaxon.find_or_create_by!(ncbi_taxon_id: 208_895, organism_id: rsv.id, lineage_caller_id: nextclade_rsvb.id,
                                 name: 'RSV-B', reference_accession: 'NC_001781.1')

lamp = AmplificationMethod.find_or_create_by!(name: 'LAMP')
qpcr = AmplificationMethod.find_or_create_by!(name: 'qPCR')

lamp_primerset = PrimerSet.find_or_create_by(name: 'NEB LAMP (E2019)', organism: sars, user: admin,
                                             amplification_method: lamp)

lamp_primers = {
  'E1-LF': 'CGCTATTAACTATTAACG',
  'E1-LB': 'GCGCTTCGATTGTGTGCGT',
  'E1-BIP': 'TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT',
  'E1-FIP': 'ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG',
  'E1-B3': 'TTCAGATTTTTAACACGAGAGT',
  'E1-F3': 'TGAGTACGAACTTATGTACTCAT',
  'N2-F3': 'ACCAGGAACTAATCAGACAAG',
  'N2-B3': 'GACTTGATCTTTGAAATTTGGATCT',
  'N2-FIP': 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC',
  'N2-BIP': 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA',
  'N2-LF': 'GGGGGCAAATTGTGCAATTTG',
  'N2-LB': 'CTTCGGGAACGTGGTTGACC'
}
lamp_primers.each_pair do |name, seq|
  Oligo.find_or_create_by(name: name, sequence: seq, primer_set: lamp_primerset)
end

luna_primerset = PrimerSet.find_or_create_by(name: 'NEB Luna qPCR (E3019)', organism: sars, user: admin,
                                             amplification_method: qpcr)

luna_primers = {
  'N1-F': 'GACCCCAAAATCAGCGAAAT',
  'N1-R': 'TCTGGTTACTGCCAGTTGAATCTG',
  'N2-F': 'TTACAAACATTGGCCGCAAA',
  'N2-R': 'GCGCGACATTCCGAAGAA',
  'N2-Probe': 'ACAATTTGCCCCCAGCGCTTCAG',
  'N1-Probe': 'ACCCCGCATTACGTTTGGTGGACC'
}

luna_primers.each_pair do |name, seq|
  Oligo.find_or_create_by(name: name, sequence: seq, primer_set: luna_primerset)
end

features = [['N', 28_274, 29_533],
            ['ORF10', 29_558, 29_674],
            ['orf1ab', 266, 21_555],
            ['ORF3a', 25_393, 26_220],
            ['S', 21_563, 25_384],
            ['M', 26_523, 27_191],
            ['E', 26_245, 26_472],
            ['ORF7a', 27_394, 27_759],
            ['ORF6', 27_202, 27_387],
            ['ORF8', 27_894, 28_259],
            ['ORF7b', 27_756, 27_887],
            ['5\'UTR', 1, 265],
            ['3\'UTR', 29_675, 29_903]]
organism = Organism.find_or_create_by(name: 'SARS-CoV-2')
features.each do |feature|
  GenomicFeature.find_or_create_by!(name: feature[0], ref_start: feature[1], ref_end: feature[2],
                                    organism_id: organism.id)
end

# Creates a list of views to operate on (views/*.sql). Order is defined by numeric prefix to view name
view_files = Dir["#{__dir__}/views/*.sql"].sort

view_defs = view_files.each_with_object({}) do |view_file, h| # depends on hash being ordered
  view_def = File.readlines(view_file)
  mat_view_matches = view_def.grep(/CREATE/)[0]&.match(/CREATE MATERIALIZED VIEW ([\w.]+) AS/)
  view_name_matches = view_def.grep(/CREATE/)[0]&.match(/CREATE VIEW ([\w.]+) AS/)
  if mat_view_matches
    h[mat_view_matches[1]] = { mat_view: view_def }
  elsif view_name_matches
    h[view_name_matches[1]] = { view: view_def }
  else
    raise "Expected VIEW or MATERIALIZED VIEW in #{view_file}"
  end
end
# puts view_defs.keys
conn = ActiveRecord::Base.connection
# drops all views in views directory, reversed to deal with view interdependencies
view_defs.keys.reverse_each do |v|
  Rails.logger.info("dropping #{v}")
  conn.execute(
    view_defs[v][:view] && "drop view if exists #{v} CASCADE;" ||
      view_defs[v][:mat_view] && "drop materialized view if exists #{v} CASCADE;"
  )
end

# creates all views
view_defs.each_key do |v|
  Rails.logger.info("Creating #{v}")
  conn.execute(
    view_defs[v][:view] && view_defs[v][:view].join.to_s ||
      view_defs[v][:mat_view] && view_defs[v][:mat_view].join.to_s
  )
end
