import "tom-select";

// Searchable dropdowns: <select class="searchable-select"> becomes a Tom Select widget. Tom Select keeps the
// original <select> in sync and dispatches its change events, so pages listen to the select itself.
const SELECTOR = 'select.searchable-select';

function enhance() {
    document.querySelectorAll(SELECTOR).forEach(select => {
        if (select.tomselect) return;
        new window.TomSelect(select, {
            plugins: select.multiple ? ['remove_button'] : [],
            maxOptions: 500,
            hidePlaceholder: true
        });
    });
}

// Restore the plain <select> before Turbo snapshots the page, or the cached copy shows a dead widget.
function teardown() {
    document.querySelectorAll(SELECTOR).forEach(select => select.tomselect?.destroy());
}

document.addEventListener('turbo:load', enhance);
document.addEventListener('turbo:before-cache', teardown);
