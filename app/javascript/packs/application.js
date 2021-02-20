// This file is automatically compiled by Webpack, along with any other files
// present in this directory. You're encouraged to place your actual application logic in
// a relevant structure within app/javascript and only use these pack files to reference
// that code so it'll be compiled.

require("@rails/ujs").start()
require("turbolinks").start()
require("@rails/activestorage").start()
require("channels")
require("@nathanvda/cocoon")

let $ = require('jquery');

require( 'datatables.net-dt/css/jquery.dataTables.min.css' );
require( 'datatables.net-dt' );
require( 'datatables.net-buttons-dt' );
require( 'datatables.net-buttons/js/buttons.html5.js' );
require( 'datatables.net-buttons/js/buttons.print.js' );
require( 'datatables.net-responsive' );
require( 'datatables.net-select-dt' );

// Load datatables styles
import 'datatables.net-dt/css/jquery.dataTables.css'

$(document).on('turbolinks:load', () => {
    if($('[id^=DataTables_Table]').length == 0) {
        $('table').DataTable({
            dom: 'lfBrtip',
            buttons: [
                'copy', 'excel'
            ],
            responsive: true
        });
    }

    // Check for click events on the navbar burger icon
    $(".navbar-burger").click(function() {
        // Toggle the "is-active" class on both the "navbar-burger" and the "navbar-menu"
        $(".navbar-burger").toggleClass("is-active");
        $(".navbar-menu").toggleClass("is-active");

    });
})
