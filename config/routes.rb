# frozen_string_literal: true

Rails.application.routes.draw do
  root 'welcome#index'
  get 'welcome/index'
  resources :organisms
  resources :oligos
  resources :primer_sets
  resources :users
  resource :user_sessions, only: %i[new create destroy]
end
