# frozen_string_literal: true

require 'test_helper'

class LineageTest < ActiveSupport::TestCase
  HEADER = 'taxon,lineage,conflict,ambiguity_score'

  setup do
    @caller = lineage_callers(:default)
    @organism = organisms(:sars_cov2)
    @existing = Lineage.create!(name: 'BA.2', caller_name: @caller.name, organism: @organism)
  end

  def calls_csv(*rows)
    file = Tempfile.new(['calls', '.csv'])
    file.write(([HEADER] + rows).join("\n") << "\n")
    file.close
    @files = (@files || []) << file
    file.path
  end

  teardown { @files&.each(&:unlink) }

  test 'parse returns only lineages not already in the database, once each' do
    path = calls_csv('OX1.1,BA.2,0.0,', 'OX2.1,XBB.1.5,0.0,', 'OX3.1,XBB.1.5,0.0,', 'OX4.1,JN.1,0.0,')

    new_lineages = Lineage.parse(path, @organism, @caller.id)

    assert_equal %w[JN.1 XBB.1.5], new_lineages.map(&:name).sort
    assert(new_lineages.all?(&:new_record?))
    assert_equal ['default'], new_lineages.map(&:caller_name).uniq
  end

  test 'parse uses a fixed number of queries regardless of row count' do
    rows = Array.new(200) { |i| "OX#{i}.1,BA.2,0.0," }
    path = calls_csv(*rows)

    # One lookup of existing names plus one for the caller. Tests run with the query cache on, but
    # upload.rb does not, so turn it off to count the queries the importer actually sends.
    Lineage.uncached do
      assert_queries_count(2) { Lineage.parse(path, @organism, @caller.id) }
    end
  end

  test 'parse raises when the file has only a header' do
    assert_raises(RuntimeError) { Lineage.parse(calls_csv, @organism, @caller.id) }
  end

  test 'ids_by_name maps known names and omits unknown ones' do
    ids = Lineage.ids_by_name(Set['BA.2', 'NOPE.1'])

    assert_equal({ 'BA.2' => @existing.id }, ids)
  end
end
