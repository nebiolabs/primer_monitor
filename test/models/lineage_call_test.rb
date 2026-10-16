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
end
