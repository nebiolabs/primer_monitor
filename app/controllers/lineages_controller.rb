# frozen_string_literal: true

class LineagesController < ApplicationController
  def index
    authorize! :index, LineagesController
    @primer_sets = {
      # '1000+_20210101-20210606' => '1000+ 20210101-20210606',
      # '3_Prime_Ends' => '3 Prime Ends',
      'ARTIC_v3' => 'ARTIC v3',
      'ARTIC_v4' => 'ARTIC v4',
      'ARTIC_v4.1' => 'ARTIC v4.1',
      'ARTIC_v4.1_(alts_only)' => 'ARTIC v4.1 (alts only)',
      'ARTIC_v4.1_(spiked_alts)' => 'ARTIC v4.1 (spiked alts)',
      'ARTIC_v5.3.2' => 'ARTIC v5.3.2',
      'ARTICv4pre' => 'ARTICv4pre',
      'All_qPCR_Primers' => 'All qPCR Primers',
      'Charité' => 'Charité',
      'China_CDC' => 'China CDC',
      'HKU' => 'HKU',
      'Japan_NIID' => 'Japan NIID',
      'LabCorp_N' => 'LabCorp N',
      'Mason_April_22' => 'Mason April 22',
      'Mason_Orf1ab2,Orf1ab4,N2' => 'Mason Orf1ab2,Orf1ab4,N2',
      'Mason_iCare' => 'Mason iCare',
      'Midnight_1200' => 'Midnight 1200',
      'Midnight_V3' => 'Midnight V3',
      'NEB_LAMP_E2019' => 'NEB LAMP E2019',
      'NEB_Luna_qPCR_E3019_CDC' => 'NEB Luna qPCR E3019 CDC',
      'NEB_Primers' => 'NEB Primers',
      'OmniSARS2' => 'OmniSARS2',
      'RGN_qPCR' => 'RGN qPCR',
      'Resende' => 'Resende',
      'Swift' => 'Swift',
      'Test_Primers' => 'Test Primers',
      'Thailand_NIH' => 'Thailand NIH',
      'UNZA' => 'UNZA',
      'US_CDC' => 'US CDC',
      'US_CDC_Flu_SC2' => 'US CDC Flu SC2',
      'USydney' => 'USydney',
      'VSL_1b_proposed' => 'VSL 1b proposed',
      'VSS2b_proposed' => 'VSS2b proposed',
      'VarSkip_(unused)' => 'VarSkip (unused)',
      'VarSkip_1a' => 'VarSkip 1a',
      'VarSkip_1a_supplements' => 'VarSkip 1a supplements',
      'VarSkip_2_(pre)' => 'VarSkip 2 (pre)',
      'VarSkip_2a_(spiked_alts)' => 'VarSkip 2a (spiked alts)',
      'VarSkip_Long_1a' => 'VarSkip Long 1a'
      # 'varskip_long_unique_regions' => 'VarSkip Long unique regions',
      # 'varskip_short_unique_regions' => 'VarSkip Short unique regions'
    }

    @config = { "data_server": ENV['IGV_DATA_SERVER'] }

  end
end
