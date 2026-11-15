# frozen_string_literal: true

module PrimerSets
  # The administrator's review decision on a primer set: PATCH /primer_sets/:primer_set_id/review with decision
  # approve or reject. Approving sets processing, which the backend publishes once every oligo is aligned.
  class ReviewsController < ApplicationController
    DECISIONS = {
      'approve' => { status: :processing, notice: 'Primer set approved. It will be published within a few minutes.' },
      'reject' => { status: :failed, notice: 'Primer set rejected. Only its submitter and administrators can see it.' }
    }.freeze

    def update
      primer_set = PrimerSet.find(params.expect(:primer_set_id))
      authorize! :confirm, primer_set
      decision = DECISIONS[params.expect(:decision)]
      return head :bad_request unless decision

      unless primer_set.update(status: decision[:status])
        return redirect_to primer_set,
                           alert: "The primer set could not be saved: #{primer_set.errors.full_messages.to_sentence}"
      end

      # a set approved with missing alignments would never be published, so align it again
      primer_set.align_primers if primer_set.processing? && primer_set.unaligned_oligos_count.positive?
      redirect_to primer_set, notice: decision[:notice]
    end
  end
end
