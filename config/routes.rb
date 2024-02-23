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

  # hardcoded legacy redirect
  get 'lineages', to: 'lineages#index'

  resources :organisms, param: :name do
    resources :lineage_variants, only: [:index]
    resources :lineages, param: :name, constraints: { name: /[A-z0-9.]+/ }
    resources :primer_sets, only: [:index]
    resource :primer_variant_summary, only: [:show]
  end

  resources :oligos
  resources :users
  resources :primer_set_subscriptions, only: [:create, :destroy]
  resources :primer_sets, only: [:new, :show, :create, :destroy, :update, :edit]
end
