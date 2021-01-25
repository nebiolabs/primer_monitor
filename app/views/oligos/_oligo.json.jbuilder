# frozen_string_literal: true

json.extract! oligo, :id, :name, :sequence, :created_at, :updated_at
json.url oligo_url(oligo, format: :json)
