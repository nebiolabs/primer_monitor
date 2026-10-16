# frozen_string_literal: true

require 'test_helper'

class OligoNameParserTest < ActiveSupport::TestCase
  {
    'varskip-0317-1_01_LEFT' => %w[01 Forward],
    'varskip-0317-1_01_LEFT_ALT3' => %w[01 Forward],
    'varskip-0317-1_12_RIGHT_Alt3' => %w[12 Reverse],
    'nCoV-2019_7_LEFT_alt0' => %w[7 Forward],
    'nCoV-2019_7_RIGHT_alt5' => %w[7 Reverse],
    'SARS-CoV-2_1_LEFT' => %w[1 Forward],
    'PM1_rsva_Amp1a_l' => %w[rsva_Amp1a Forward],
    'PM1_rsva_Amp1a_l_alt1' => %w[rsva_Amp1a Forward],
    'PM1_rsvb_Amp3b-alt1_l' => %w[rsvb_Amp3b Forward],
    'PM1_rsva_Amp7a_alt_r' => %w[rsva_Amp7a Reverse],
    'N1-F' => %w[N1 Forward],
    'N1-R' => %w[N1 Reverse],
    'N2-Probe' => %w[N2 Probe],
    '2019-nCoV_N1-P' => %w[N1 Probe],
    'N1_R' => %w[N1 Reverse],
    'E1-FIP' => %w[E1 FIP],
    'E1-F3' => %w[E1 F3],
    'N2-LB' => %w[N2 LB],
    'N1-F some description' => %w[N1 Forward]
  }.each do |name, (locus, category)|
    test "parses #{name}" do
      assert_equal({ locus:, category: }, OligoNameParser.parse(name))
    end
  end

  test 'returns nils when unsure' do
    %w[probe1 MS1 nCoV_X].each do |name|
      assert_equal({ locus: nil, category: nil }, OligoNameParser.parse(name), name)
    end
    assert_equal({ locus: nil, category: nil }, OligoNameParser.parse(nil))
  end
end
