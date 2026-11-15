# frozen_string_literal: true

module PrimerSetsHelper
  STATUS_COLORS = { 'created' => 'is-warning', 'processing' => 'is-info', 'complete' => 'is-success',
                    'failed' => 'is-danger' }.freeze

  def primer_set_status_tag(primer_set)
    tag.span(primer_set.status_label, class: ['tag', 'is-light', STATUS_COLORS[primer_set.status.to_s]])
  end

  def primer_set_status_notification_class(primer_set)
    ['notification', 'is-light', STATUS_COLORS[primer_set.status.to_s]]
  end
end
