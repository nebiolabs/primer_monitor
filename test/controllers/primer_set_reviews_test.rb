# frozen_string_literal: true

require 'test_helper'

# Approving and rejecting primer sets (PrimerSets::ReviewsController) and the review panel on the show page
class PrimerSetReviewsTest < ActionDispatch::IntegrationTest
  setup do
    @primer_set = primer_sets(:one)
    sign_in(users(:admin_user))
  end

  test 'an administrator approves a set awaiting review from its page' do
    Organism.any_instance.stubs(:primer_sets_config).returns([{}, {}])
    get primer_set_url(@primer_set)

    assert_select '#review_status', text: /waiting for your review/
    assert_select 'form[action=?] button', primer_set_review_path(@primer_set), text: 'Approve and publish'

    patch primer_set_review_url(@primer_set), params: { decision: 'approve' }

    assert_redirected_to primer_set_url(@primer_set)
    assert_predicate @primer_set.reload, :processing?
  end

  test 'approving a set with unaligned oligos runs the alignment again' do
    OligoAlignmentPosition.delete_all
    PrimerSet.any_instance.expects(:align_primers).once
    patch primer_set_review_url(@primer_set), params: { decision: 'approve' }

    assert_predicate @primer_set.reload, :processing?
  end

  test 'an administrator rejects a set' do
    patch primer_set_review_url(@primer_set), params: { decision: 'reject' }

    assert_redirected_to primer_set_url(@primer_set)
    assert_predicate @primer_set.reload, :failed?
  end

  test 'a submitter can neither approve their own set nor set its status through the form' do
    sign_out(users(:admin_user))
    sign_in(users(:one))

    patch primer_set_review_url(@primer_set), params: { decision: 'approve' }

    assert_predicate @primer_set.reload, :created?

    patch primer_set_url(@primer_set), params: { primer_set: { status: 'processing', user_id: users(:two).id } }

    assert_predicate @primer_set.reload, :created?
    assert_equal users(:one), @primer_set.user
  end

  test 'a submitter sees that their set is waiting for review, without the review buttons' do
    sign_out(users(:admin_user))
    sign_in(users(:one))
    Organism.any_instance.stubs(:primer_sets_config).returns([{}, {}])
    get primer_set_url(@primer_set)

    assert_select '#review_status', text: /An administrator will review/
    assert_select 'form[action=?]', primer_set_review_path(@primer_set), count: 0
  end

  test 'an unknown decision is a bad request and leaves the status alone' do
    patch primer_set_review_url(@primer_set), params: { decision: 'publish' }

    assert_response :bad_request
    assert_predicate @primer_set.reload, :created?
  end
end
