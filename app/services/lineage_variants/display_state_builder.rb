# frozen_string_literal: true

module LineageVariants
  class DisplayStateBuilder
    def self.call(organism:, track_data:, url_lineage:, url_primer_sets:)
      new(organism:, track_data:, url_lineage:, url_primer_sets:).build
    end

    def initialize(organism:, track_data:, url_lineage:, url_primer_sets:)
      @organism = organism
      @track_data = track_data
      @url_lineage = url_lineage
      @url_primer_sets = url_primer_sets
    end

    def build
      {
        lineage_frequencies:,
        sorted_lineages:,
        effective_lineage: url_lineage || dominant_lineage || valid_lineage_default,
        effective_primer_sets:
      }
    end

    private

    attr_reader :organism, :track_data, :url_lineage, :url_primer_sets

    def lineage_sets
      track_data[:lineage_sets]
    end

    def lineage_frequencies
      @lineage_frequencies ||= organism.grouped_lineage_frequencies(lineage_sets)
    end

    def sorted_lineages
      @sorted_lineages ||= begin
        keys_set = lineage_sets.keys.to_set
        children_of, roots = build_lineage_tree(lineage_sets.keys, keys_set)
        result = flatten_lineage_tree(roots, children_of)
        result << ['all', 0] if lineage_sets.key?('all')
        result
      end
    end

    def effective_primer_sets
      url_primer_sets ||
        track_data[:default_tracks].presence ||
        track_data[:primer_sets].keys.select { |key| key.start_with?('NEB') }.presence ||
        track_data[:primer_sets].keys
    end

    def dominant_lineage
      lineage_frequencies.reject { |key, value| key == 'all' || value.zero? }.max_by { |_, value| value }&.first
    end

    def valid_lineage_default
      default = track_data[:default_lineage]
      return default if default && lineage_sets.key?(default)
      return 'all' if lineage_sets.key?('all')

      sorted_lineages.first&.first
    end

    def build_lineage_tree(keys, keys_set)
      children_of = Hash.new { |hash, key| hash[key] = [] }
      roots = []
      keys.reject { |key| key == 'all' }.each do |key|
        parent = closest_ancestor_key(key, keys_set)
        (parent ? children_of[parent] : roots) << key
      end
      children_of.each_value { |values| values.sort_by! { |key| -lineage_frequencies.fetch(key, 0) } }
      roots.sort_by! { |key| -lineage_frequencies.fetch(key, 0) }
      [children_of, roots]
    end

    def flatten_lineage_tree(roots, children_of, depth = 0)
      roots.flat_map { |key| [[key, depth]] + flatten_lineage_tree(children_of.fetch(key, []), children_of, depth + 1) }
    end

    def closest_ancestor_key(key, keys_set)
      parts = key.split('.')
      (parts.length - 1).downto(1) do |index|
        prefix = parts[0, index].join('.')
        return prefix if keys_set.include?(prefix)
      end
      nil
    end
  end
end
