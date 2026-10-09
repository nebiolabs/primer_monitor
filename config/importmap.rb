# Pin npm packages by running ./bin/importmap

# Turbo and the Rails libraries come from their gems, so they always match the installed versions
pin "@hotwired/turbo-rails", to: "turbo.min.js"
pin "@rails/activestorage", to: "activestorage.esm.js"

pin "tom-select" # @2.6.2 (vendored tom-select.complete.min.js; defines window.TomSelect)
pin "igv" # @3.8.9
pin "jszip" # @3.10.1 (vendored dist/jszip.min.js; defines window.JSZip for the Excel button)
pin_all_from "app/javascript"
pin "datatables.net" # @3.1.3 (dependency-free; no jQuery)
pin "datatables.net-bm" # @3.1.3
pin "datatables.net-buttons" # @4.1.2 (includes the copy/csv/excel buttons)
pin "datatables.net-buttons-bm" # @4.1.2
pin "datatables.net-responsive" # @4.1.1
pin "datatables.net-responsive-bm" # @4.1.1
