class PrimerSet < ApplicationRecord
    belongs_to :user
    has_many :oligos

end
