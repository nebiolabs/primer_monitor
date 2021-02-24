# frozen_string_literal: true
# This file should contain all the record creation needed to seed the database with its default values.
# The data can then be loaded with the rails db:seed command (or created alongside the database with db:setup).


#password should be changed...
admin = User.create_with(first: 'Admin', last:'User',
                         password: Rails.application.credentials.admin_password,
                         password_confirmation: Rails.application.credentials.admin_password,
                         active: true, approved: true, confirmed: true)
            .find_or_create_by!(email: 'primer-monitor@neb.com')

%w[administrator pi operator participant].each do |role_name|
  Role.find_or_create_by!(name: role_name)
end

if admin.roles.size == 0
  admin.roles << Role.find_by(name: 'administrator')
end

sars = Organism.find_or_create_by!(name: "SARS-CoV-2", ncbi_taxon_id:2697049)

lamp_primerset = PrimerSet.find_or_create_by(name: "NEB LAMP (E2019)", organism: sars, user: admin )

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
  'N2-LB': 'CTTCGGGAACGTGGTTGACC',
}
lamp_primers.each_pair do |name,seq|
  Oligo.find_or_create_by(name: name, sequence: seq, primer_set: lamp_primerset)
end

luna_primerset = PrimerSet.find_or_create_by(name: "NEB Luna qPCR (E3019)", organism: sars, user: admin )

luna_primers = {
  'N1-F': 'GACCCCAAAATCAGCGAAAT',
  'N1-R': 'TCTGGTTACTGCCAGTTGAATCTG',
  'N2-F': 'TTACAAACATTGGCCGCAAA',
  'N2-R': 'GCGCGACATTCCGAAGAA',
  'N2-Probe': 'ACAATTTGCCCCCAGCGCTTCAG',
  'N1-Probe': 'ACCCCGCATTACGTTTGGTGGACC'
}

luna_primers.each_pair do |name,seq|
  Oligo.find_or_create_by(name: name, sequence: seq, primer_set: luna_primerset)
end



features = [['N',28274,29533],
['ORF10',29558,29674],
['orf1ab',266,21555],
['ORF3a',25393,26220],
['S',21563,25384],
['M',26523,27191],
['E',26245,26472],
['ORF7a',27394,27759],
['ORF6',27202,27387],
['ORF8',27894,28259],
['ORF7b',27756,27887],
['5\'UTR',1,265],
['3\'UTR',29675,29903],
]
organism = Organism.find_or_create_by(name: 'SARS-CoV-2')
features.each do |feature|
  GenomicFeature.create!(name: feature[0], ref_start: feature[1], ref_end: feature[2], organism_id: organism.id)
end



