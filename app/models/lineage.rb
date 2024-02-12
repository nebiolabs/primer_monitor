# frozen_string_literal: true

class Lineage < ApplicationRecord
  belongs_to :organism
  has_many :lineage_calls, dependent: :destroy
  has_many :fasta_records, through: :lineage_calls

  def self.parse(calls_csv, organism, caller_id)
    raise "Unable to find calls file #{calls_csv}" unless File.exist?(calls_csv)

    new_lineage_names = Set[] # new empty set
    external_links = {}
    record_count = 0

    File.readlines(calls_csv).each do |line|
      next if line.start_with?('taxon,')

      record_count += 1
      new_lineage, external_link = parse_record(line)
      unless new_lineage.nil?
        new_lineage_names.add new_lineage
        external_links[new_lineage] = external_link if external_link.present?
      end
    end

    # if there are no records (not only pre-existing lineages, but nothing at all)
    raise "Unable to parse any records from #{calls_csv}" if record_count.zero?

    new_lineages = []

    new_lineage_names.each do |lineage_name|
      new_lineages << Lineage.new(name: lineage_name, caller_name: LineageCaller.find_by(id: caller_id).name,
                                  organism_id: organism.id, external_link: external_links[lineage_name])
    end

    new_lineages
  end

  def self.parse_record(line)
    fields = line.chomp.split(',')
    ActiveRecord::Base.logger.info fields if line.chomp.split("\t")[1].blank? || fields.nil?

    # prevent trying to re-add lineages that are already in the database
    if Lineage.exists? name: fields[1]
      nil
    else
      [fields[2], fields[1]]
    end
  end
end
