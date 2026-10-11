# frozen_string_literal: true

# "Multiplex PCR" covered the whole-genome tiling schemes (ARTIC, VarSkip, Midnight, ...), but many qPCR panels are
# multiplexed too. "Tiled amplicon" names what they share — overlapping amplicons covering the genome for
# sequencing — and links to the method paper instead of NEB's general PCR page.
class RenameMultiplexPcrToTiledAmplicon < ActiveRecord::Migration[8.1]
  OLD = { name: 'Multiplex PCR',
          url: 'https://www.neb.com/en-us/applications/dna-amplification-pcr-and-qpcr/pcr' }.freeze
  NEW = { name: 'Tiled amplicon', url: 'https://doi.org/10.1038/nprot.2017.066' }.freeze

  def up
    rename(OLD, NEW)
  end

  def down
    rename(NEW, OLD)
  end

  private

  def rename(from, to)
    execute <<~SQL.squish
      UPDATE amplification_methods SET name = #{quote(to[:name])}, description_url = #{quote(to[:url])}
      WHERE name = #{quote(from[:name])}
    SQL
  end
end
