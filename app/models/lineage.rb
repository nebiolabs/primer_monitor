# frozen_string_literal: true

class Lineage < ApplicationRecord
  belongs_to :organism
  has_many :lineage_calls, dependent: :destroy
  has_many :fasta_records, through: :lineage_calls

  def self.external_link_url(organism, lineage_name)
    # expand this with new organisms as support is added for them
    case organism.slug
    when 'sars-cov-2'
      "https://outbreak.info/situation-reports?pango=#{lineage_name}"
    end
  end

  # Returns unsaved Lineage records for the lineage names in calls_csv that are not yet in the database.
  def self.parse(calls_csv, organism, caller_id)
    names = lineage_names(calls_csv)

    # Look up every name in one query. A per-row lookup costs one database round trip per call.
    new_lineage_names = names - Lineage.where(name: names.to_a).distinct.pluck(:name)

    caller_name = LineageCaller.find_by(id: caller_id).name
    new_lineage_names.map do |lineage_name|
      Lineage.new(name: lineage_name, caller_name:, organism_id: organism.id,
                  external_link: external_link_url(organism, lineage_name))
    end
  end

  # Returns the set of lineage names (second column) in calls_csv, skipping the header.
  # Raises if the file holds no records at all.
  def self.lineage_names(calls_csv)
    raise "Unable to find calls file #{calls_csv}" unless File.exist?(calls_csv)

    names = Set[]
    record_count = 0
    File.foreach(calls_csv) do |line|
      next if line.start_with?('taxon,', 'seqName,')

      record_count += 1
      name = line.chomp.split(',')[1]
      if name.blank?
        logger.warn("No lineage name in #{calls_csv} line: #{line.chomp}")
      else
        names.add(name)
      end
    end
    raise "Unable to parse any records from #{calls_csv}" if record_count.zero?

    names
  end

  # Maps each of the given lineage names to its id, using one query.
  # If a name exists more than once, the lowest id wins.
  def self.ids_by_name(names)
    where(name: names.to_a).order(:id).pluck(:name, :id).each_with_object({}) do |(name, id), ids|
      ids[name] ||= id
    end
  end
end
