# frozen_string_literal: true

class Lineage < ApplicationRecord
  belongs_to :organism
  has_many :pangolin_calls
  has_many :fasta_records, through: :pangolin_calls

  def self.parse(pangolin_csv)
    raise "Unable to find calls file #{pangolin_csv}" unless File.exist?(pangolin_csv)

    new_lineage_names = Set[] # new empty set
    record_count = 0

    File.readlines(pangolin_csv).each do |line|
      next if line.start_with?('taxon,')

      record_count += 1
      if line.chomp.split("\t")[1].blank? || line.chomp.split(',')[1].nil?
        ActiveRecord::Base.logger.info line.chomp.split(',')
      end

      # prevent trying to re-add lineages that are already in the database
      next if Lineage.exists? name: line.chomp.split(',')[1]

      new_lineage_names.add line.chomp.split(',')[1]
    end

    # if there are no records (not only pre-existing lineages, but nothing at all)
    raise "Unable to parse any records from #{pangolin_csv}" if record_count.zero?

    new_lineages = []

    sarscov2_id = Organism.find_by(ncbi_taxon_id: 2_697_049).id # hardcoded SARS-CoV-2 id for now

    new_lineage_names.each do |lineage_name|
      new_lineages << Lineage.new(name: lineage_name, caller_name: 'pangolin', organism_id: sarscov2_id)
    end

    new_lineages
  end
end
