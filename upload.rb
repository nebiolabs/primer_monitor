#!/usr/bin/env ruby
# frozen_string_literal: true

# aggregates results from the passed files and stores them in the configured database
require 'rubygems'
require 'bundler/setup'
require 'active_record'
require 'activerecord-import'
require 'active_record/associations'
require 'logger'
require 'pg'
require 'yaml'
require 'pp'
# require 'slop'
# require 'digest/md5'
# require 'date'
# require 'tempfile'
# require 'tmpdir'
# require 'csv'

# Use 'safe' in order to prevent rails from failing due to monkey
# patching Enumerable:
require 'descriptive_statistics/safe'

# need application record first since others depend on it
require_relative "#{__dir__}/app/models/application_record.rb"

# requires all the model files
Dir["#{__dir__}/app/models/*.rb"].each do |f|
  require_relative f
end

def setup_db_connection
  db_config = YAML.safe_load(
    File.open(
      File.join(
        File.dirname(__FILE__), 'config', 'database.yml'
      )
    )
  )[ENV['DATABASE_ENV'] || 'production']
  ActiveRecord::Base.establish_connection(db_config)
end

def parse_options
    # ADDING A NEW OPTION? TODO CHECKLIST:

    # Does the option specify an input file? Add it to input_file_params
    # hash/method.

    # Is the option required for read_group (library)? Add it to
    # required_read_group_params hash/method, otherwise to
    # optional_read_group_params.

    # Does the option processing change the data for the view in
    # dna_production_quality_metrics (the view is in Tableau for DNA
    # production QC table)? See, for example,
    # db/migrate/20190530213731_chg_dna_production_quality_metrics_to_matviews.rb. Add
    # to matview_for_option_str multiline string/table both the option
    # and the corresponding materialized view(s) to be refreshed when
    # the option is processed.

    begin
      slop_opts = Slop.parse(ARGV.map(&:strip)) do |o|
        o.string '--metadata_tsv', 'Tsv file containing metadata information'
        o.string '--variants_tsv', 'Tsv file containing variants information'

        # The available log levels are: :debug, :info, :warn, :error, and
        # :fatal, corresponding to the log level numbers from 0 up to 4
        # respectively. See rails docs.
        o.string '--verbose', 'Verbosity level of ActiveRecord logger', default: 'INFO'

        o.on '--help' do
          puts o
          exit
        end
      end
    rescue Slop::MissingArgument => e
      @log.error "fatal: #{e}"
      exit(1)
    rescue Slop::UnknownOption => e
      @log.error "fatal: #{e}"
      exit(1)
    rescue Slop::UnknownOption => e
      @log.error 'fatal: ' + String(e)
      exit(1)
    end
    slop_opts.to_hash
  end

def import_metadata(metadata_file)
  return unless metadata_file

  metadata_records = FastaRecord.parse(metadata_file)
  FastaRecord.import(metadata_records, validate: false)
end

def import_variants(variants_file)
  return unless variants_file

  variant_records = VariantSite.parse(variants_file)
  fasta_ids = variant_records.map(&:fasta_record_id).uniq
  fasta_ids.each do |fasta_id|
    VariantSite.where(fasta_record_id: fasta_id).delete_all
  end
  VariantSite.import(variant_records, validate: false)
end

def main
  opts = parse_options
  @log.level = Logger.const_get(opts[:verbose])

  # @log.debug 'checking samtools is installed properly...'
  setup_db_connection

  # import_metadata('./test/fixtures/metadata.tsv')
  # import_variants('./test/fixtures/variants.tsv')
  import_metadata(opts[:metadata])
  import_variants(opts[:variants])

end
main if $PROGRAM_NAME.end_with?('upload.rb')
