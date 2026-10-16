# frozen_string_literal: true

namespace :oligos do
  desc 'Fill blank oligo locus/category from oligo names (all sets, or the given ids). Dry run unless APPLY=1'
  task :backfill_locus_category, [:primer_set_id] => :environment do |_, args|
    apply = ENV['APPLY'] == '1'
    sets = PrimerSet.includes(:oligos).order(:id)
    sets = sets.where(id: [args[:primer_set_id], *args.extras]) if args[:primer_set_id]

    sets.each do |primer_set|
      blank = primer_set.oligos.select { |o| o.locus.blank? || o.category.blank? }
      next if blank.empty?

      blank.each(&:valid?) # fills locus/category from the name
      filled = blank.select { |o| o.locus_changed? || o.category_changed? }
      unparsed = blank - filled
      puts "#{primer_set.id} #{primer_set.name}: #{filled.size} of #{blank.size} blank oligos filled, " \
           "#{filled.filter_map(&:locus).uniq.size} loci"
      filled.each { |o| puts "  #{o.name} -> #{o.locus.inspect} #{o.category.inspect}" } if ENV['VERBOSE'] == '1'
      puts "  unparsed: #{unparsed.map(&:name).join(', ')}" if unparsed.any?
      filled.each(&:save!) if apply
    end
    puts apply ? 'Saved.' : 'Dry run; rerun with APPLY=1 to save.'
  end
end
