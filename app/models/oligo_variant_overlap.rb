# frozen_string_literal: true

# read-only view of each usable variant overlapping a primer in the last 168 days, refreshed nightly; the same rows
# the notification alerts are built from
class OligoVariantOverlap < ApplicationRecord
  self.table_name = 'oligo_variant_overlaps'

  def readonly?
    true
  end
end
