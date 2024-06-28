import 'init_jquery'

import select2 from "select2";
select2();

import Rails from '@rails/ujs';
Rails.start();

import "@hotwired/turbo-rails"

import * as ActiveStorage from "@rails/activestorage"
ActiveStorage.start()

let tables = []; //list of datatables to destroy on turbo:unload

import 'datatables.net';

import "@nathanvda/cocoon";

//old requires for reference until everything is tested
/*
require( 'datatables.net-dt/css/jquery.dataTables.min.css' );
require( 'datatables.net-dt' );
require( 'datatables.net-buttons-dt' );
require( 'datatables.net-buttons/js/buttons.html5.js' );
require( 'datatables.net-buttons/js/buttons.print.js' );
require( 'datatables.net-responsive' );
require( 'datatables.net-select-dt' );
require("channels");
*/

// Configure your import map in config/importmap.rb. Read more: https://github.com/rails/importmap-rails

$(document).on('turbo:load', () => {
    $('select.wide-select2').select2({width: '80%'});
    $('select.select2').select2();


    if($('[id^=DataTables_Table]').length == 0 || true) {
        let table = new DataTable('table', {
            dom: 'lfBrtip',
            buttons: [
                'copy', 'excel'
            ],
            responsive: true
        });
        tables.push(table);
        //doing an array so this can later be made to support multiple tables per page
    }

    // Check for click events on the navbar burger icon
    $(".navbar-burger").click(function() {
        // Toggle the "is-active" class on both the "navbar-burger" and the "navbar-menu"
        $(".navbar-burger").toggleClass("is-active");
        $(".navbar-menu").toggleClass("is-active");

    });
});

$(document).on('turbo:before-cache', () => {
    tables.forEach((table) => {
        table.destroy();
    });
    tables = [];
});
