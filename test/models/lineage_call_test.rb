# frozen_string_literal: true

require 'test_helper'

class LineageCallTest < ActiveSupport::TestCase
  setup do
    @caller = lineage_callers(:default)
    @lineage = Lineage.create!(name: 'BA.2', caller_name: @caller.name, organism: organisms(:sars_cov2))
    @rec_one = fasta_records(:one)
    @rec_two = fasta_records(:two)
    @rec_one.update!(genbank_accession: 'OX000001.1')
    @rec_two.update!(genbank_accession: 'OX000002.1')
  end

  def new_call(accession)
    LineageCall.create!(taxon: accession, lineage: @lineage, lineage_caller_id: @caller.id)
  end

  test 'pending update overwrites a stale pending call, only for the given calls' do
    stale = new_call('OX000001.1')
    @rec_one.update!(pending_lineage_call_id: stale.id)
    fresh = new_call('OX000001.1')
    other = new_call('OX000002.1') # not in this batch

    LineageCall.update_fasta_recs(true, [fresh.id])

    assert_equal fresh.id, @rec_one.reload.pending_lineage_call_id
    assert_nil @rec_two.reload.pending_lineage_call_id, "call #{other.id} was outside the batch"
  end

  test 'current update only fills records without a call' do
    existing = new_call('OX000001.1')
    @rec_one.update!(lineage_call_id: existing.id)
    newer_one = new_call('OX000001.1')
    new_two = new_call('OX000002.1')

    LineageCall.update_fasta_recs(false, [newer_one.id, new_two.id])

    assert_equal existing.id, @rec_one.reload.lineage_call_id
    assert_equal new_two.id, @rec_two.reload.lineage_call_id
  end

  test 'empty batch is a no-op' do
    assert_nothing_raised { LineageCall.update_fasta_recs(true, []) }
  end

  def calls_csv(*rows)
    @csv = Tempfile.new(['calls', '.csv'])
    @csv.write((['taxon,lineage,conflict,ambiguity_score'] + rows).join("\n") << "\n")
    @csv.close
    @csv.path
  end

  teardown { @csv&.unlink }

  test 'parse builds a call per row with the lineage id' do
    calls = LineageCall.parse(calls_csv('OX000001.1,BA.2,0.0,', 'OX000002.1,BA.2,0.5,0.9'), @caller.id)

    assert_equal %w[OX000001.1 OX000002.1], calls.map(&:taxon)
    assert_equal [@lineage.id], calls.map(&:lineage_id).uniq
    assert_equal [@caller.id], calls.map(&:lineage_caller_id).uniq
  end

  test 'parse uses one query regardless of row count' do
    path = calls_csv(*Array.new(200) { |i| "OX#{i}.1,BA.2,0.0," })

    # upload.rb runs without the query cache that tests have, so turn it off to count real queries.
    Lineage.uncached do
      assert_queries_count(1) { LineageCall.parse(path, @caller.id) }
    end
  end

  test 'parse raises on a lineage missing from the database' do
    error = assert_raises(RuntimeError) { LineageCall.parse(calls_csv('OX000001.1,NOPE.1,0.0,'), @caller.id) }
    assert_match 'Failed to find lineage NOPE.1', error.message
  end
end
