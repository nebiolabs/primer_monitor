// Configure your import map in config/importmap.rb. Read more: https://github.com/rails/importmap-rails
import "@hotwired/turbo-rails"

import * as ActiveStorage from "@rails/activestorage"
ActiveStorage.start()

import 'searchable_select';
import 'nested_fields';

import 'data_tables';

// Toggle the "is-active" class on both the navbar burger and the navbar menu (Bulma's mobile menu)
document.addEventListener('click', (event) => {
    if (!event.target.closest('.navbar-burger')) return;
    document.querySelectorAll('.navbar-burger, .navbar-menu').forEach(el => el.classList.toggle('is-active'));
});
