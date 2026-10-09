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
  devise_scope :user do
    post 'users/password_reset', to: 'users/registrations#send_password_reset', as: :user_password_reset
  end

  root 'welcome#index'
  get 'about', to: 'about#show'
  get 'history', to: 'history#show'

  # hardcoded legacy redirect
  get 'lineages', to: 'lineages#index'

  resources :organisms, param: :slug do
    resources :lineage_variants, only: [:index] do
      get :variant_overlaps, on: :collection
    end
    resources :lineages, param: :name, constraints: { name: /[A-z0-9.]+/ }
    resources :primer_sets, only: [:index]
  end

  resources :oligos
  # accounts are created by signing up (Devise) or with Google/Microsoft, so admins only list, edit and delete
  resources :users, except: %i[new create]
  resources :primer_set_subscriptions, only: [:create, :destroy]
  resources :primer_sets, only: [:new, :show, :create, :destroy, :update, :edit]
end
