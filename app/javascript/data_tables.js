import DataTable from 'datatables.net-bm';
import 'datatables.net-buttons-bm';
import 'datatables.net-responsive-bm';
import 'datatables.net-columncontrol-bm';
import 'jszip';

// Every table on a page gets sorting, search, paging, Copy/Excel export and responsive collapsing
// (DataTables 3, no jQuery). Tables marked data-column-filters also get a header dropdown per column for sorting and
// picking which values to show. Header cells with class filter-search get a typed search instead (number and date
// comparisons, or text for columns with too many distinct values to list); no-filter cells (download links, sequences,
// action buttons) only sort.
// Tables marked data-datatable="off" are left alone: the lineage variant table pairs each row with hidden sub-rows
// that sorting and paging would separate, and a Turbo restore brings it back already built.
DataTable.Buttons.jszip(window.JSZip);

let tables = [];

const columnFilterOptions = {
    columnControl: ['order', ['orderAsc', 'orderDesc', 'orderClear', 'searchList']],
    columnDefs: [
        { targets: '.no-filter', columnControl: ['order'] },
        { targets: '.filter-search', columnControl: ['order', ['orderAsc', 'orderDesc', 'orderClear', 'search']] }
    ],
    // ColumnControl draws the sort indicator and handles header clicks itself
    ordering: { indicators: false, handler: false }
};

document.addEventListener('turbo:load', () => {
    document.querySelectorAll('table:not([data-datatable="off"])').forEach(table => {
        if (DataTable.isDataTable(table)) return;
        table.classList.add('table', 'is-striped', 'is-hoverable', 'is-fullwidth'); // Bulma table styling
        tables.push(new DataTable(table, {
            layout: {
                topStart: ['pageLength', 'buttons'],
                topEnd: 'search',
                bottomStart: 'info',
                bottomEnd: 'paging'
            },
            // small, to match the page-length dropdown beside them
            buttons: [{ extend: 'copy', className: 'is-small' }, { extend: 'excel', className: 'is-small' }],
            responsive: true,
            ...('columnFilters' in table.dataset ? columnFilterOptions : {})
        }));
    });
});

// Restore the plain tables before Turbo snapshots the page
document.addEventListener('turbo:before-cache', () => {
    tables.forEach(table => table.destroy());
    tables = [];
});
