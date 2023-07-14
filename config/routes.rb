# frozen_string_literal: true

Rails.application.routes.draw do
  devise_for :users, controllers: {
    sessions: 'users/sessions',
    registrations: 'users/registrations',
    confirmations: 'users/confirmations',
    unlocks: 'users/unlocks',
    passwords: 'users/passwords',
    omniauth_callbacks: 'users/omniauth_callbacks'
  }

  root 'welcome#index'
  get 'help', to: 'help#show'
  get 'about', to: 'about#show'
  get 'history', to: 'history#show'
  get 'lineage_variants', to: 'lineage_variants#index'
  resources :organisms, param: :name do
    resources :lineages, param: :name, constraints: { name: /[A-Z]+(\.\d+)*/ }
  end

  get 'lineages', to: 'lineages#index'

  resources :oligos
  resources :primer_sets
  resources :users
  resources :primer_set_subscriptions, only: [:create, :destroy]
end
