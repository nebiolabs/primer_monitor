# frozen_string_literal: true

class LineagesController < ApplicationController
  def index
    authorize! :index, LineagesController
    if params.key? :organism_name
      @organism = Organism.find_by(slug: params[:organism_name])

      query = <<-SQL
        SELECT name, times_seen, last_seen FROM lineage_info WHERE organism_id=?;
      SQL
      @lineages = ActiveRecord::Base.connection.execute(ActiveRecord::Base.sanitize_sql([query, @organism.id])).to_a
    else
      # hardcoded legacy redirect
      redirect_to organism_lineage_variants_url(Organism.find_by(slug: 'sars-cov-2').slug), status: :moved_permanently
    end
  end

  def show
    authorize! :show, LineagesController
    @organism = Organism.find_by(name: params[:organism_name])
    @lineage = Lineage.find_by(name: params[:name])
    query = <<-SQL
      SELECT times_seen, first_seen, last_seen FROM lineage_info WHERE organism_id=? AND name=?;
    SQL
    # first because there should only be one record
    @lineage_info = ActiveRecord::Base.connection.execute(ActiveRecord::Base.sanitize_sql([query, @organism.id,
                                                                                           @lineage.name])).to_a.first
  end

  def to_param
    name
  end
end
