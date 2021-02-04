# frozen_string_literal: true

json.extract! organism, :id, :name, :ncbi_taxon_id, :created_at, :updated_at
json.url organism_url(organism, format: :json)
