class Organism < ApplicationRecord
    has_many :blast_hits
    has_many :amplicons
end
