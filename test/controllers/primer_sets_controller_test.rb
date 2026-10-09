# frozen_string_literal: true

require 'test_helper'

class PrimerSetsControllerTest < ActionDispatch::IntegrationTest
  setup do
    @primer_set = primer_sets(:one)
    sign_in(users(:admin_user))
  end

  test 'should get index' do
    Organism.any_instance.stubs(:primer_sets_config).returns([{}, {}])
    get organism_primer_sets_url(organisms(:sars_cov2))

    assert_response :success
  end

  test 'should get new' do
    get new_primer_set_url

    assert_response :success
  end

  test 'should create primer_set' do
    assert_difference('PrimerSet.count') do
      post primer_sets_url, params: { primer_set: {
        user_id: @primer_set.user_id,
        name: 'test',
        organism_id: organisms(:sars_cov2).id,
        amplification_method_id: amplification_methods(:qPCR).id,
        oligos_attributes: [{ name: 'F', sequence: 'ATCG' },
                            { name: 'R', sequence: 'CGTA' }]
      } }
    end

    assert_redirected_to edit_primer_set_url(PrimerSet.last)
  end

  test 'should create primer_set with oligos as the form submits them (keyed by row)' do
    assert_difference('Oligo.count', 2) do
      post primer_sets_url, params: { primer_set: {
        user_id: @primer_set.user_id,
        name: 'test',
        organism_id: organisms(:sars_cov2).id,
        amplification_method_id: amplification_methods(:qPCR).id,
        oligos_attributes: { '0' => { name: 'F', sequence: 'ATCG' }, '1' => { name: 'R', sequence: 'CGTA' } }
      } }
    end
  end

  test 'should show primer_set' do
    Organism.any_instance.stubs(:primer_sets_config).returns([{}, {}])
    get primer_set_url(@primer_set)

    assert_response :success
  end

  test 'show marks FASTA and BED pending until the primer set is on the data server' do
    @primer_set.update!(status: :complete)
    Organism.any_instance.stubs(:primer_sets_config).returns([{ data_server: 'http://data' }, {}])
    get primer_set_url(@primer_set)

    assert_select 'strong', text: 'FASTA:' do |strong|
      assert_includes strong.first.parent.text, 'Pending'
    end
    assert_select 'strong', text: 'BED:' do |strong|
      assert_includes strong.first.parent.text, 'Pending'
    end
    assert_select '#igv', count: 0
  end

  test 'show links FASTA and BED once the primer set is on the data server' do
    @primer_set.update!(status: :complete)
    Organism.any_instance.stubs(:primer_sets_config)
            .returns([{ data_server: 'http://data', organism_slug: 'sars-cov-2' }, { @primer_set.name => 'Charite' }])
    get primer_set_url(@primer_set)

    assert_select 'a[href=?]', 'http://data/sars-cov-2/primer_sets_fasta/Charite.fasta'
    assert_select 'a[href=?]', 'http://data/sars-cov-2/primer_sets_bed/Charite.bed'
    assert_select '#igv'
  end

  test 'should get edit' do
    get edit_primer_set_url(@primer_set)

    assert_response :success
  end

  test 'should update primer_set' do
    patch primer_set_url(@primer_set), params:
      { primer_set: { user_id: @primer_set.user_id, name: "#{@primer_set.name}∆" } }

    assert_redirected_to edit_primer_set_url(@primer_set)
  end

  test 'should destroy primer_set' do
    assert_difference('PrimerSet.count', -1) do
      delete primer_set_url(@primer_set)
    end

    assert_redirected_to primer_sets_url
  end
end
