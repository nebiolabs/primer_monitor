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
    @config[:initial_lineage] = @effective_lineage
    @config[:initial_primer_sets] = @effective_primer_sets
  end

  def variant_overlaps
    authorize! :index, LineageVariantsController

    @organism = Organism.find_by(slug: params[:organism_slug])
    return head :not_found unless @organism

    lineage_param = params[:lineage].presence
    primer_set_params = Array(params[:primer_sets]).reject(&:blank?)

    return head :bad_request if lineage_param.nil? || primer_set_params.empty?

    rows = LineageVariantPrimerOverlap
      .joins(oligo: [:primer_set, { oligo_alignment_positions: :organism_taxon }])
      .where(
        lineage_variant_primer_overlaps: {
          organism_id: @organism.id,
          lineage_group_key: lineage_param
        },
        primer_sets: { name: primer_set_params },
        organism_taxa: { organism_id: @organism.id }
      )
      .select(
        'lineage_variant_primer_overlaps.ref_start',
        'lineage_variant_primer_overlaps.ref_end',
        'lineage_variant_primer_overlaps.variant_type',
        'lineage_variant_primer_overlaps.variant',
        'lineage_variant_primer_overlaps.frequency_pct',
        'lineage_variant_primer_overlaps.first_seen_date',
        'lineage_variant_primer_overlaps.first_seen_lineage',
        'lineage_variant_primer_overlaps.first_seen_location',
        'lineage_variant_primer_overlaps.last_seen_date',
        'lineage_variant_primer_overlaps.last_seen_lineage',
        'lineage_variant_primer_overlaps.last_seen_location',
        'oligos.name AS oligo_name',
        'oligos.sequence AS oligo_sequence',
        'oligos.strand AS oligo_strand',
        'oligo_alignment_positions.ref_start AS oligo_start',
        'oligo_alignment_positions.ref_end AS oligo_end',
        'primer_sets.name AS primer_set_name'
      )

    variants_map = {}
    rows.each do |row|
      key = [row.ref_start, row.ref_end, row.variant_type, row.variant]
      variants_map[key] ||= {
        ref_start: row.ref_start,
        ref_end: row.ref_end,
        variant_type: row.variant_type,
        variant: row.variant,
        frequency_pct: row.frequency_pct,
        first_seen: { date: row.first_seen_date, lineage: row.first_seen_lineage, location: row.first_seen_location },
        last_seen:  { date: row.last_seen_date,  lineage: row.last_seen_lineage,  location: row.last_seen_location  },
        oligos: []
      }
      variants_map[key][:oligos] << {
        name: row.oligo_name,
        sequence: row.oligo_sequence,
        oligo_start: row.oligo_start,
        oligo_end: row.oligo_end,
        strand: row.oligo_strand,
        primer_set: row.primer_set_name
      }
    end

    render json: { variants: variants_map.values }
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
    @effective_lineage = @url_lineage || dominant_lineage || valid_lineage_default
    @effective_primer_sets = @url_primer_sets ||
                             @track_data[:default_tracks].presence ||
                             @track_data[:primer_sets].keys.select { |k| k.start_with?('NEB') }.presence ||
                             @track_data[:primer_sets].keys
  end

  def dominant_lineage
    @lineage_frequencies.reject { |k, v| k == 'all' || v.zero? }.max_by { |_, v| v }&.first
  end

  def valid_lineage_default
    default = @track_data[:default_lineage]
    return default if default && @track_data[:lineage_sets].key?(default)
    return 'all' if @track_data[:lineage_sets].key?('all')

    @sorted_lineages.first&.first
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
