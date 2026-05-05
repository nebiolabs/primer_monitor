# frozen_string_literal: true

class LineageVariantsController < ApplicationController
  def index
    authorize! :index, LineageVariantsController

    @organism = Organism.find_by(slug: params[:organism_slug])
    return head :not_found unless @organism

    @config = {
      "data_server": ENV['IGV_DATA_SERVER'],
      "organism_slug": @organism.slug,
      "organism_name": @organism.name
    }

    @track_data = @organism.lineage_variants_data(@config[:data_server], @config[:organism_slug])

    return unless @track_data[:data_fetched]

    @url_lineage = params[:lineage].presence
    @url_primer_sets = params[:primer_sets]&.split(',')

    prepare_lineage_display
  end

  private

  def prepare_lineage_display
    individual_frequencies = @organism.recent_lineage_frequencies
    lineage_sets = @track_data[:lineage_sets]

    @lineage_frequencies = lineage_sets.keys.to_h do |k|
      prefix = "#{k}."
      [k, individual_frequencies.filter { |name, _| name == k || name.start_with?(prefix) }.sum { |_, count| count }]
    end

    @sorted_lineages = sorted_lineage_entries(lineage_sets)
    @effective_lineage = @url_lineage || dominant_lineage || @track_data[:default_lineage]
    @effective_primer_sets = @url_primer_sets ||
                             @track_data[:default_tracks].presence ||
                             @track_data[:primer_sets].keys.select { |k| k.start_with?('NEB') }
  end

  def dominant_lineage
    @lineage_frequencies.reject { |k, v| k == 'all' || v.zero? }.max_by { |_, v| v }&.first
  end

  def sorted_lineage_entries(lineage_sets)
    keys_set = lineage_sets.keys.to_set
    children_of, roots = build_lineage_tree(lineage_sets.keys, keys_set)
    result = flatten_lineage_tree(roots, children_of)
    result << ['all', 0] if lineage_sets.key?('all')
    result
  end

  def build_lineage_tree(keys, keys_set)
    children_of = Hash.new { |h, k| h[k] = [] }
    roots = []
    keys.reject { |k| k == 'all' }.each do |k|
      parent = closest_ancestor_key(k, keys_set)
      (parent ? children_of[parent] : roots) << k
    end
    children_of.each_value { |v| v.sort_by! { |k| -@lineage_frequencies.fetch(k, 0) } }
    roots.sort_by! { |k| -@lineage_frequencies.fetch(k, 0) }
    [children_of, roots]
  end

  def flatten_lineage_tree(roots, children_of, depth = 0)
    roots.flat_map { |k| [[k, depth]] + flatten_lineage_tree(children_of.fetch(k, []), children_of, depth + 1) }
  end

  def closest_ancestor_key(key, keys_set)
    parts = key.split('.')
    (parts.length - 1).downto(1) do |i|
      prefix = parts[0, i].join('.')
      return prefix if keys_set.include?(prefix)
    end
    nil
  end
end
