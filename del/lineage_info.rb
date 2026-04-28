# frozen_string_literal: true

# this model represents the lineage info materialized view that is updated by a
# cron job as new sequences are added to the database
class LineageInfo < ApplicationRecord
  self.table_name = 'lineage_info'

  def readonly?
    true
  end
end
