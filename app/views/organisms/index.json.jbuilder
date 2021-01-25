# frozen_string_literal: true

json.array! @organisms, partial: 'organisms/organism', as: :organism
