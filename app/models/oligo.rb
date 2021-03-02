# frozen_string_literal: true

class Oligo < ApplicationRecord
  belongs_to :primer_set
  has_many :blast_hits
  validates :short_name, length: 1..5, allow_blank: true
  after_validation :align_oligo

  enum category: { F3: 'F3', FIP: 'FIP', LF: 'LF', LB: 'LB', BIP: 'BIP', B3: 'B3',
                   Forward: 'Forward', Probe: 'Probe', Reverse: 'Reverse' }


  PATH_TO_BOWTIE2 = <path>
  PATH_TO_BOWTIE2_INDEX = <path>

  def to_s
    "#{name}: #{sequence}"
  end

  def align_oligo
    #do aligning
    alignment = `echo #{self.sequence} | #{PATH_TO_BOWTIE2} --ma 3 --mp 2 --np 0 --xeq -N 1 -L 7 -r -k 2 --local --quiet --no-head -x #{PATH_TO_BOWTIE2_INDEX} -U /dev/stdin | cut -f 2,3,4,6`

    strand_flag = alignment.split("\t")[0]
    ref_contig = alignment.split("\t")[1]
    start_position = alignment.split("\t")[2]
    cigar = alignment.split("\t")[3]

    cigar_lengths = cigar.scan(/\d+/)
    cigar_types = cigar.scan(/[HSID=X]/)

    if ref_contig == "*"
      # didn't align, break
    end

    # Determine if the primer aligned to F or R strand
    if (strand_flag & 0x10) == 16
      strand = "R"
    else
      strand = "F"
    end

    if cigar =~ "H|S|I|D"
      # Has softclip or indels, warning maybe?
    end

    align_length = 0
    cigar_types.each_with_index { |element, index|
      if element =~ "=|X|D"
        align_length += cigar_lengths[index]
      elsif element =~ "I"
        align_length -+ cigar_lengths[index]
      else
        # Soft or hard clip
        next
      end
    }
    end_position = start_position + align_length

    self.ref_start = start_position
    self.ref_end = end_position
    self.strand = strand
  end
end
