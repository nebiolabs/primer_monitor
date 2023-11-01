#!/usr/bin/env ruby
# frozen_string_literal: true

require_relative 'config/environment' # loads rails env, models, etc.

# aggregates results from the passed files and stores them in the configured database

@log = ActiveRecord::Base.logger = Logger.new($stderr)

def setup_db_connection
  db_config = YAML.safe_load(
    File.open(
      File.join(
        File.dirname(__FILE__), 'config', 'database.yml'
      )
    ), aliases: true
  )[ENV['DATABASE_ENV'] || 'production']
  ActiveRecord::Base.establish_connection(db_config)
end

def define_options
  Slop.parse(ARGV.map(&:strip)) do |o|
    o.string '--metadata_tsv', 'TSV file containing metadata information'
    o.string '--variants_tsv', 'TSV file containing variants information'
    o.string '--lineage_csv', 'CSV file containing lineage calls'
    o.string '--organism', 'Organism slug for finding organism ID'
    o.bool '--pending', 'Save calls as pending', default: false
    o.bool '--import_calls', 'Performs lineage data import', default: false
    o.bool '--import_seqs', 'Performs sequence/variant data import', default: false
    o.bool '--rebuild_views', 'Rebuilds materialized views', default: false

    # The available log levels are: :debug, :info, :warn, :error, and
    # :fatal, corresponding to the log level numbers from 0 up to 4
    # respectively. See rails docs.
    o.string '--verbose', 'Verbosity level of ActiveRecord logger', default: 'INFO'

    o.on '--help', '-h' do
      puts o
      exit
    end
  end
end

def parse_options
  begin
    slop_opts = define_options
  rescue Slop::MissingArgument => e
    @log.error "fatal: #{e}"
    exit(1)
  rescue Slop::UnknownOption => e
    @log.error "fatal: #{e}"
    exit(1)
  end
  slop_opts.to_hash
end

def import_metadata(metadata_file, organism)
  return unless metadata_file

  @log.info("starting import: #{metadata_file}: ")

  fasta_records = FastaRecord.parse(metadata_file, organism)
  result = FastaRecord.import(fasta_records, validate: false, on_duplicate_key_ignore: true)
  result.failed_instances.each { |rec| @log.error("Failed to insert \"#{rec}\"") }
  @log.info("Loaded #{result.ids.size}/#{fasta_records.size} new fasta records")
end

def import_variants(variants_file, organism)
  return unless variants_file

  variant_records = VariantSite.parse(variants_file, organism)
  # fasta_ids = variant_records.map(&:fasta_record_id).uniq
  # VariantSite.where(fasta_record_id: fasta_ids).delete_all
  VariantSite.import(variant_records, validate: false, on_duplicate_key_ignore: true)
end

def create_new_notifications!
  new_proposed_notifications = ProposedNotification.new_proposed_notifications()
  return if new_proposed_notifications.empty?

  ProposedNotification.import(new_proposed_notifications, validate: false)
end

def import_lineages(lineage_csv, pending, organism)
  return unless lineage_csv

  @log.info("starting import: #{lineage_csv}: ")

  lineages = Lineage.parse(lineage_csv, organism)
  ActiveRecord::Base.transaction do
    ActiveRecord::Base.connection.execute('LOCK lineages IN ROW EXCLUSIVE MODE')
    result_lineages = Lineage.import(lineages, on_duplicate_key_ignore: true)
    result_lineages.failed_instances.each { |rec| @log.error("Failed to insert lineage \"#{rec}\"") }
    @log.info("Loaded #{result_lineages.ids.size}/#{lineages.size} new lineages")
  end

  calls = LineageCall.parse(lineage_csv)
  result_calls = LineageCall.import(calls, validate: false)
  result_calls.failed_instances.each { |rec| @log.error("Failed to insert lineage call \"#{rec}\"") }
  @log.info("Loaded #{result_calls.size}/#{calls.size} new lineage calls")

  LineageCall.update_fasta_recs pending
end

def main
  opts = parse_options
  @log.level = Logger.const_get(opts[:verbose])

  setup_db_connection
  ActiveRecord::Base.transaction do
    if opts[:import_seqs]
      import_metadata(opts[:metadata_tsv], opts[:organism])
      import_variants(opts[:variants_tsv], opts[:organism])
    end
    import_lineages(opts[:lineage_csv], opts[:pending], opts[:organism]) if opts[:import_calls]
    if opts[:rebuild_views]
      %w[variant_overlaps counts time_counts oligo_variant_overlaps
         identify_primers_for_notifications initial_score lineage_info].each do |view|
        ActiveRecord::Base.connection.execute("REFRESH MATERIALIZED VIEW #{view}")
        @log.info("refreshing #{view}")
      end
      create_new_notifications!
      # cannot send notifications here since this may be running on a cluster node.
    end
  end
end
main if $PROGRAM_NAME.end_with?('upload.rb')
