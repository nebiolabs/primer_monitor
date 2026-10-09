# Pin npm packages by running ./bin/importmap

# Turbo and the Rails libraries come from their gems, so they always match the installed versions
pin "@hotwired/turbo-rails", to: "turbo.min.js"
pin "@rails/ujs", to: "rails-ujs.esm.js"
pin "@rails/activestorage", to: "activestorage.esm.js"

pin "jquery" # @3.7.1
pin "select2" # @4.1.0
pin "@nathanvda/cocoon", to: "@nathanvda--cocoon.js", preload: false # @1.2.14
pin "igv" # @3.8.9
pin_all_from "app/javascript"
pin "datatables.net", to: "https://cdn.datatables.net/v/dt/dt-1.12.0/b-1.7.1/b-html5-1.7.1/b-print-1.7.1/r-2.3.0/sl-1.4.0/datatables.min.js"
