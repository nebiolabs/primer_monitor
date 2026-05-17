# frozen_string_literal: true

module VariantRefBackfillTask
  module_function

  def backfill_ref_for_taxon(taxon)
    ref_seq = VariantSite.load_reference_sequence(taxon)
    unless ref_seq
      puts "Skipping #{taxon.name}: FASTA not found at expected path"
      return
    end

    puts "Processing #{taxon.name} (#{ref_seq.length} bp)..."

    conn = ActiveRecord::Base.connection
    raw = conn.raw_connection

    with_ref_backfill_table(raw) do
      copy_ref_backfill_rows(raw, taxon, ref_seq)
      update_variant_site_refs(conn, taxon)
      update_lineage_overlap_refs(conn, taxon)
    end
  end

  def with_ref_backfill_table(raw)
    raw.exec('CREATE TEMP TABLE IF NOT EXISTS _ref_backfill ' \
             '(ref_start integer, ref_end integer, variant_type varchar, ref varchar, ' \
             'PRIMARY KEY (ref_start, ref_end, variant_type))')
    raw.exec('TRUNCATE _ref_backfill')
    yield
  ensure
    raw.exec('DROP TABLE _ref_backfill')
  end

  def copy_ref_backfill_rows(raw, taxon, ref_seq)
    raw.copy_data('COPY _ref_backfill (ref_start, ref_end, variant_type, ref) FROM STDIN') do
      copy_snp_rows(raw, ref_seq)
      copy_variant_rows(raw, taxon, ref_seq, 'I')
      copy_variant_rows(raw, taxon, ref_seq, 'D')
    end
  end

  def copy_snp_rows(raw, ref_seq)
    ref_seq.length.times do |pos|
      ref = VariantSite.reference_for_range(ref_seq, 'X', pos, pos + 1)
      raw.put_copy_data("#{pos}\t#{pos + 1}\tX\t#{ref}\n")
    end
  end

  def copy_variant_rows(raw, taxon, ref_seq, variant_type)
    VariantSite.pending_ref_ranges(taxon.id, variant_type).each do |ref_start, ref_end|
      ref = VariantSite.reference_for_range(ref_seq, variant_type, ref_start.to_i, ref_end.to_i)
      next unless ref

      raw.put_copy_data("#{ref_start}\t#{ref_end}\t#{variant_type}\t#{ref}\n")
    end
  end

  def update_variant_site_refs(conn, taxon)
    result = conn.execute(<<~SQL)
      UPDATE variant_sites
      SET ref = b.ref
      FROM _ref_backfill b
      WHERE variant_sites.ref_start    = b.ref_start
        AND variant_sites.ref_end      = b.ref_end
        AND variant_sites.variant_type = b.variant_type
        AND variant_sites.organism_taxon_id = #{taxon.id}
        AND variant_sites.ref IS NULL
    SQL
    puts "  Updated #{result.cmd_tuples} variant_sites rows"
  end

  def update_lineage_overlap_refs(conn, taxon)
    result = conn.execute(<<~SQL)
      UPDATE lineage_variant_primer_overlaps lvpo
      SET ref = b.ref
      FROM _ref_backfill b
      WHERE lvpo.ref_start    = b.ref_start
        AND lvpo.ref_end      = b.ref_end
        AND lvpo.variant_type = b.variant_type
        AND lvpo.organism_id  = #{taxon.organism_id}
        AND lvpo.ref IS NULL
    SQL
    puts "  Updated #{result.cmd_tuples} lineage_variant_primer_overlaps rows"
  end
end

namespace :variants do
  desc 'Backfill ref column in variant_sites from reference FASTA files. Safe to re-run.'
  task backfill_ref: :environment do
    OrganismTaxon.includes(:organism).find_each { |taxon| VariantRefBackfillTask.backfill_ref_for_taxon(taxon) }
  end
end
