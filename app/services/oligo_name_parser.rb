# frozen_string_literal: true

# Infers an oligo's locus (amplicon or target) and category from its name, since FASTA uploads carry only names.
# Handles tiled amplicon schemes (varskip-0317-1_01_LEFT_ALT3, nCoV-2019_7_RIGHT_alt5, PM1_rsva_Amp1a_l_alt1)
# and qPCR/LAMP target-suffix names (2019-nCoV_N1-F, N2-Probe, E1-FIP). Fields it can't infer come back nil.
module OligoNameParser
  module_function

  ALT_MARKER = /[_-]alt\d*(?=[_-]|\z)/i
  TILED_SIDES = { 'LEFT' => 'Forward', 'RIGHT' => 'Reverse' }.freeze
  TILED_SHORT_SIDES = { 'L' => 'Forward', 'R' => 'Reverse' }.freeze
  # RSV panels mix taxa in one set; keep the taxon so RSV-A and RSV-B amplicons stay distinct
  TAXON_TOKEN = /\A(rsva|rsvb|panr)\z/i

  SUFFIX_CATEGORIES = {
    'F3' => 'F3', 'B3' => 'B3', 'FIP' => 'FIP', 'BIP' => 'BIP', 'LF' => 'LF', 'LB' => 'LB',
    'F' => 'Forward', 'FP' => 'Forward', 'FWD' => 'Forward', 'FORWARD' => 'Forward',
    'R' => 'Reverse', 'RP' => 'Reverse', 'REV' => 'Reverse', 'REVERSE' => 'Reverse',
    'P' => 'Probe', 'PR' => 'Probe', 'PROBE' => 'Probe'
  }.freeze

  def parse(name)
    base = name.to_s.split.first.to_s.gsub(ALT_MARKER, '')
    return { locus: nil, category: nil } if base.empty?

    parse_tiled(base) || parse_suffix(base) || { locus: nil, category: nil }
  end

  def parse_tiled(base)
    tokens = base.split('_')
    side_index = tokens.rindex { |t| TILED_SIDES.key?(t.upcase) }
    category = side_index && TILED_SIDES[tokens[side_index].upcase]
    if category.nil? && tokens.size > 2 && TILED_SHORT_SIDES.key?(tokens.last.upcase)
      side_index = tokens.size - 1
      category = TILED_SHORT_SIDES[tokens.last.upcase]
    end
    return unless category

    { locus: tiled_amplicon(tokens[0...side_index]), category: }
  end

  def tiled_amplicon(before_side)
    amplicon = before_side.last
    return if amplicon.blank? || before_side.size < 2

    taxon = before_side[-2]
    taxon&.match?(TAXON_TOKEN) ? "#{taxon}_#{amplicon}" : amplicon
  end

  # target then category at the end of the name: N1-F, N2-Probe, E1-FIP, 2019-nCoV_N1-P, N1_R
  def parse_suffix(base)
    match = base.match(/(?:\A|_)(?<locus>[A-Za-z0-9.]+)[-_](?<suffix>[A-Za-z0-9]+)\z/)
    category = match && SUFFIX_CATEGORIES[match[:suffix].upcase]
    return unless category

    { locus: match[:locus], category: }
  end
end
