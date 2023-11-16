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
  get 'about', to: 'about#show'
  get 'history', to: 'history#show'
  resources :organisms, param: :name do
    get 'lineage_variants', to: 'lineage_variants#index'
    resources :lineages, param: :name, constraints: { name: /[A-Z]+(\.\d+)*/ }
    resources :primer_sets, only: [:index]
    resources :primer_variant_summaries, only: [:index]
  end

  get 'lineages', to: 'lineages#index'

  resources :oligos
  resources :users
  resources :primer_set_subscriptions, only: [:create, :destroy]
  resources :primer_sets, only: [:new, :show, :create, :destroy, :update]
end
