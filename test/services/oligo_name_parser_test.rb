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

  {
    'SARS-CoV-2-1200_10_LEFT' => '10L',
    'SARS-CoV-2-1200_10_RIGHT' => '10R',
    'varskip-0317-1_01_LEFT_ALT3' => '01L_3',
    'SARS-CoV-2_10_RIGHT_alt1' => '10R_1',
    'nCoV-2019_7_LEFT_alt0' => '7L_0',
    'SARS-CoV-2_400_10_LEFT_0' => '10L',
    'SARS-CoV-2_400_10_LEFT_1' => '10L_1',
    'PM1_rsva_Amp1a_l' => '1aL',
    'PM1_rsva_Amp15a_l_alt1' => '15aL1',
    'PM1_panr_Amp11_r' => '11R',
    'PM1_rsva_Amp7a_alt_r' => '7aR_a',
    'varskip-0317-1_100_RIGHT_ALT12' => '100R1',
    'N1-F' => 'N1-F',
    '2019-nCoV_N1-P' => 'N1-P',
    'E1-FIP' => 'E1FIP',
    'covid19genome_200-29703_s15860_U_42R' => 'U_42R',
    'covid19genome_200-29703_s15860_U_100R' => 'U100R',
    'MS1' => 'MS1'
  }.each do |name, short|
    test "short name for #{name}" do
      assert_equal short, OligoNameParser.short_name(name)
    end
  end

  test 'short name is nil when the name gives nothing distinctive' do
    assert_nil OligoNameParser.short_name('probe12345')
    assert_nil OligoNameParser.short_name(nil)
  end
end
