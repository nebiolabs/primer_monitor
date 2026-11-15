import { addNestedRow, removeNestedRow } from 'nested_fields';

// The oligo table on the primer set form. It sorts, filters and bulk-edits the rows in place. All rows stay in the
// form whatever the sort or filter, so every oligo is submitted (a paged table would drop the rows it isn't showing).
//
// Rows are .nested-fields <tr>s in tbody#samples; each editable cell has data-col naming its column. All listeners
// are delegated from the document, so they register once and survive Turbo navigation and failed-save re-renders.

const COLUMNS = ['name', 'short_name', 'locus', 'category', 'sequence'];
// the bases Oligo accepts (IUPAC codes and whitespace)
const INVALID_BASES = /[^ATUCGWSMKRYBDHVN\s]/i;
// numeric so amplicon 2 sorts before amplicon 10
const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });

let sortState = { key: null, dir: 1 };
let submitting = false;

const editor = () => document.querySelector('[data-oligo-editor]');
const tbody = () => document.getElementById('samples');
const rows = () => [...tbody().querySelectorAll(':scope > tr.nested-fields')];
const cell = (row, col) => row.querySelector(`[data-col="${col}"]`);
const isRemoved = row => row.classList.contains('is-removed');
const isShown = row => !row.hidden;
const sequenceLength = row => cell(row, 'sequence').value.replace(/\s/g, '').length;

function cellValue(row, key) {
    return key === 'length' ? sequenceLength(row) : cell(row, key).value;
}

function isEdited(element) {
    if (element.tagName === 'SELECT') return [...element.options].some(o => o.selected !== o.defaultSelected);
    return element.value !== element.defaultValue;
}

// Sets a cell as if the user had typed it. A category must name one of the select's options (any case), else the
// cell is left alone.
function setCell(row, col, value) {
    const element = cell(row, col);
    if (!element) return;
    if (col === 'category') {
        const option = [...element.options].find(o => o.value.toLowerCase() === value.toLowerCase());
        if (option) element.value = option.value;
    } else {
        element.value = col === 'short_name' ? value.slice(0, 5) : value;
    }
}

function refreshRow(row) {
    const sequence = cell(row, 'sequence');
    row.querySelector('[data-oligo-length]').textContent = sequence.value ? sequenceLength(row) : '';
    sequence.classList.toggle('is-danger', INVALID_BASES.test(sequence.value));
    sequence.title = INVALID_BASES.test(sequence.value) ? 'Contains characters that are not bases' : '';
    const persisted = row.dataset.persisted === 'true';
    row.classList.toggle('is-new', !persisted);
    row.classList.toggle('is-changed', persisted && COLUMNS.some(col => isEdited(cell(row, col))));
}

function refreshSummary() {
    const all = rows();
    const kept = all.filter(row => !isRemoved(row));

    // flag names used by more than one oligo
    const nameCounts = new Map();
    kept.forEach(row => {
        const name = cell(row, 'name').value.trim();
        if (name) nameCounts.set(name, (nameCounts.get(name) || 0) + 1);
    });
    all.forEach(row => {
        const name = cell(row, 'name');
        const duplicate = !isRemoved(row) && nameCounts.get(name.value.trim()) > 1;
        name.classList.toggle('is-warning', duplicate);
        name.title = duplicate ? 'Another oligo has this name' : '';
    });

    const count = cls => all.filter(row => row.classList.contains(cls)).length;
    const shown = all.filter(isShown).length;
    const parts = [`${kept.length} oligos`];
    if (shown < all.length) parts.push(`${shown} shown`);
    if (count('is-changed')) parts.push(`${count('is-changed')} changed`);
    if (count('is-new')) parts.push(`${count('is-new')} new`);
    if (count('is-removed')) parts.push(`${count('is-removed')} removed`);
    editor().querySelector('[data-oligo-summary]').textContent = parts.join(' · ');

    const selected = selectedRows().length;
    editor().querySelector('[data-oligo-selected-count]').textContent = `${selected} selected`;
    const selectAll = editor().querySelector('[data-select-all]');
    const shownRows = all.filter(isShown);
    selectAll.checked = shownRows.length > 0 && selected === shownRows.length;
    selectAll.indeterminate = selected > 0 && selected < shownRows.length;
}

function refreshAll() {
    rows().forEach(refreshRow);
    refreshSummary();
}

function selectedRows() {
    return rows().filter(row => isShown(row) && row.querySelector('[data-oligo-select]').checked);
}

function sortRows(key, dir) {
    sortState = { key, dir };
    const template = tbody().querySelector('template');
    rows()
        .sort((a, b) => {
            const va = cellValue(a, key);
            const vb = cellValue(b, key);
            return (key === 'length' ? va - vb : collator.compare(va, vb)) * dir;
        })
        .forEach(row => tbody().insertBefore(row, template));
    editor().querySelectorAll('th[data-sort-key]').forEach(th => {
        if (th.dataset.sortKey === key) th.setAttribute('aria-sort', dir > 0 ? 'ascending' : 'descending');
        else th.removeAttribute('aria-sort');
    });
}

function applyFilter() {
    const query = editor().querySelector('[data-oligo-filter]').value.trim().toLowerCase();
    rows().forEach(row => {
        row.hidden = !!query && !COLUMNS.some(col => cellValue(row, col).toLowerCase().includes(query));
        // bulk actions only reach rows the user can see
        if (row.hidden) row.querySelector('[data-oligo-select]').checked = false;
    });
    refreshSummary();
}

function bulkAction(action) {
    const selected = selectedRows();
    if (selected.length === 0) {
        alert('Check the rows to change first.');
        return;
    }
    selected.forEach(row => {
        if (action === 'category') setCell(row, 'category', editor().querySelector('[data-bulk-category]').value);
        else if (action === 'locus') setCell(row, 'locus', editor().querySelector('[data-bulk-locus]').value.trim());
        else if (action === 'short_name') setCell(row, 'short_name', '');
        else if (action === 'remove') removeNestedRow(row, true);
        if (row.isConnected) refreshRow(row);
    });
    refreshSummary();
}

function nextEditableRow(row) {
    let next = row.nextElementSibling;
    while (next && !(next.matches('tr.nested-fields') && isShown(next) && !isRemoved(next))) {
        next = next.nextElementSibling;
    }
    return next?.matches('tr.nested-fields') ? next : null;
}

// Several lines pasted into a cell fill that column down the following rows, adding rows as needed. Tab-separated
// values fill the cells to the right, in table order: name, short name, locus, category, sequence.
function pasteIntoRows(input, text) {
    const lines = text.replace(/\r/g, '').replace(/\n+$/, '').split('\n').map(line => line.split('\t'));
    const firstCol = COLUMNS.indexOf(input.dataset.col);
    let row = input.closest('tr');
    lines.forEach((values, i) => {
        if (i > 0) row = nextEditableRow(row) || addNestedRow(tbody());
        values.forEach((value, j) => {
            const col = COLUMNS[firstCol + j];
            if (col) setCell(row, col, value.trim());
        });
        refreshRow(row);
    });
    refreshSummary();
}

// FASTA records replace the sequence of the oligo with the same name, and the rest are added as new rows
function loadFasta(text, fileName) {
    if (text.trimStart()[0] !== '>') {
        alert(`'${fileName}' is not a valid FASTA file.`);
        return;
    }
    const byName = new Map(rows().filter(row => !isRemoved(row))
        .map(row => [cell(row, 'name').value.trim(), row]));
    let updated = 0;
    let added = 0;
    let invalid = 0;
    text.trimStart().slice(1).split(/\r?\n>/).forEach(record => {
        const [header, ...sequenceLines] = record.split(/\r?\n/);
        const name = header.trim();
        const sequence = sequenceLines.join('').replace(/\s/g, '');
        if (!name && !sequence) return;
        if (INVALID_BASES.test(sequence)) invalid++;
        let row = byName.get(name);
        if (row) {
            updated++;
        } else {
            row = addNestedRow(tbody());
            setCell(row, 'name', name);
            byName.set(name, row);
            added++;
        }
        setCell(row, 'sequence', sequence);
        refreshRow(row);
    });
    refreshSummary();
    const notes = [`${added} added`, `${updated} updated`];
    if (invalid) notes.push(`${invalid} with non-base characters (highlighted)`);
    showFileName(`${fileName}: ${notes.join(', ')}`);
}

function showFileName(text) {
    let label = document.getElementById('fasta_name');
    if (!label) {
        label = document.createElement('span');
        label.id = 'fasta_name';
        label.className = 'file-name';
        document.getElementById('fasta_upload_label').appendChild(label);
        document.getElementById('fasta_div').classList.add('has-name');
    }
    label.textContent = text;
}

function hasUnsavedOligoChanges() {
    return rows().some(row => ['is-changed', 'is-removed'].some(cls => row.classList.contains(cls)) ||
        (row.classList.contains('is-new') && COLUMNS.some(col => cell(row, col).value)));
}

function init() {
    const root = editor();
    if (!root || root.dataset.ready) return;
    root.dataset.ready = 'true';
    submitting = false;
    sortRows('name', 1);
    refreshAll();
}

document.addEventListener('input', event => {
    if (!editor()?.contains(event.target)) return;
    if (event.target.matches('[data-oligo-filter]')) {
        applyFilter();
        return;
    }
    const row = event.target.closest('tr.nested-fields');
    if (row) {
        refreshRow(row);
        refreshSummary();
    }
});

document.addEventListener('change', event => {
    const target = event.target;
    if (target.id === 'fasta_upload' && target.files.length) {
        const file = target.files[0];
        file.text().then(text => loadFasta(text, file.name));
        target.value = ''; // lets the same file be uploaded again
        return;
    }
    if (!editor()?.contains(target)) return;
    if (target.matches('[data-select-all]')) {
        rows().filter(isShown).forEach(row => { row.querySelector('[data-oligo-select]').checked = target.checked; });
    }
    const row = target.closest('tr.nested-fields');
    if (row) refreshRow(row);
    refreshSummary();
});

document.addEventListener('click', event => {
    if (!editor()?.contains(event.target)) return;
    const header = event.target.closest('th[data-sort-key]');
    if (header) {
        const key = header.dataset.sortKey;
        sortRows(key, sortState.key === key ? -sortState.dir : 1);
        return;
    }
    const bulk = event.target.closest('[data-bulk-action]');
    if (bulk) bulkAction(bulk.dataset.bulkAction);
});

document.addEventListener('paste', event => {
    const input = event.target;
    if (!input.matches?.('#samples input[data-col]')) return;
    const text = event.clipboardData.getData('text/plain');
    // a single value pastes normally
    if (!/[\n\t]/.test(text.replace(/[\r\n]+$/, ''))) return;
    event.preventDefault();
    pasteIntoRows(input, text);
});

// Enter moves down the column, like a spreadsheet, instead of submitting the form
document.addEventListener('keydown', event => {
    if (event.key !== 'Enter' || !event.target.matches?.('#samples input[data-col]')) return;
    event.preventDefault();
    const next = nextEditableRow(event.target.closest('tr'));
    cell(next || event.target.closest('tr'), event.target.dataset.col)?.focus();
});

document.addEventListener('nested:change', event => {
    if (!editor()?.contains(event.target)) return;
    const row = event.detail.row;
    if (row.isConnected) refreshRow(row);
    refreshSummary();
});

document.addEventListener('submit', event => {
    if (event.target.contains(editor())) submitting = true;
});

// a failed save re-renders the form, and the user may leave it after that
document.addEventListener('turbo:submit-end', event => {
    if (!event.detail.success) submitting = false;
});

window.addEventListener('beforeunload', event => {
    if (!editor() || submitting || !hasUnsavedOligoChanges()) return;
    event.preventDefault();
    event.returnValue = '';
});

document.addEventListener('turbo:before-visit', event => {
    if (!editor() || submitting || !hasUnsavedOligoChanges()) return;
    if (!confirm('Leave this page and lose your unsaved oligo changes?')) event.preventDefault();
});

document.addEventListener('turbo:load', init);
document.addEventListener('turbo:render', init);
init();
