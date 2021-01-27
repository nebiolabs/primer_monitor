# frozen_string_literal: true

Rails.application.routes.draw do
  root 'welcome#index'
  get 'help', to: 'help#show'
  get 'about', to: 'about#show'
  get 'history', to: 'history#show'
  resources :organisms
  resources :oligos
  resources :primer_sets
  resources :users
  resource :user_sessions, only: %i[new create destroy]
  resources :user_email_confirmations, only: :show
end
