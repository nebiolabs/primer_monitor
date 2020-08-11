class Oligo < ApplicationRecord
    belongs_to :primer_set
    has_many :blast_hits
end
