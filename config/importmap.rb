# Pin npm packages by running ./bin/importmap

pin "jquery" # @3.7.1
pin "select2" # @4.1.0
pin "@rails/activestorage", to: "@rails--activestorage.js" # @7.1.3
pin "@rails/ujs", to: "@rails--ujs.js" # @7.1.3
pin "@nathanvda/cocoon", to: "@nathanvda--cocoon.js", preload: false # @1.2.14
pin "@hotwired/turbo", to: "@hotwired--turbo.js" # @8.0.4
pin "@hotwired/turbo-rails", to: "@hotwired--turbo-rails.js" # @8.0.4
pin "@rails/actioncable/src", to: "@rails--actioncable--src.js" # @7.1.3
pin "@hotwired/turbo-rails", to: "turbo.min.js"
pin "igv" # @2.15.11
pin_all_from "app/javascript"
pin "datatables.net", to: "https://cdn.datatables.net/v/dt/dt-1.12.0/b-1.7.1/b-html5-1.7.1/b-print-1.7.1/r-2.3.0/sl-1.4.0/datatables.min.js"
