// Add/remove rows of a Rails nested form (accepts_nested_attributes_for), replacing cocoon.
//
//   <div data-nested-fields>
//     <div class="nested-fields">... f.hidden_field :_destroy ... <button data-nested-remove></button></div>
//     <template data-nested-template> fields_for ... child_index: 'NEW_RECORD' </template>
//     <button type="button" data-nested-add>Add</button>
//   </div>
//
// An add button outside the container names it with a selector instead: data-nested-add="#samples". This is how
// a table body holds the rows while its add button sits below the table.
//
// Each added row gets a unique index in place of NEW_RECORD. Removing a saved row (data-persisted="true") marks it
// with the is-removed class and sets _destroy, and removing it again undoes that. Removing a row that was never saved
// deletes it. Either way the container gets a bubbling nested:change event.

let counter = 0;

export function addNestedRow(container) {
    const template = container.querySelector('template[data-nested-template]');
    const html = template.innerHTML.replace(/NEW_RECORD/g, `${Date.now()}${counter++}`);
    template.insertAdjacentHTML('beforebegin', html);
    const rows = container.querySelectorAll('.nested-fields');
    const row = rows[rows.length - 1];
    container.dispatchEvent(new CustomEvent('nested:change', { bubbles: true, detail: { row, action: 'add' } }));
    return row;
}

// Removes a row, or toggles the removal of a saved one. Pass removed to set the state instead of toggling it.
export function removeNestedRow(row, removed = undefined) {
    const container = row.parentElement;
    if (row.dataset.persisted === 'true') {
        const remove = removed ?? !row.classList.contains('is-removed');
        row.querySelector('input[name$="[_destroy]"]').value = remove ? '1' : '0';
        row.classList.toggle('is-removed', remove);
        const button = row.querySelector('[data-nested-remove]');
        if (button) button.title = remove ? 'Undo remove' : 'Remove';
    } else {
        row.remove();
    }
    container.dispatchEvent(new CustomEvent('nested:change', { bubbles: true, detail: { row, action: 'remove' } }));
}

document.addEventListener('click', (event) => {
    const add = event.target.closest('[data-nested-add]');
    if (add) {
        event.preventDefault();
        const target = add.dataset.nestedAdd;
        addNestedRow(target ? document.querySelector(target) : add.closest('[data-nested-fields]'));
        return;
    }

    const remove = event.target.closest('[data-nested-remove]');
    if (!remove) return;
    event.preventDefault();
    removeNestedRow(remove.closest('.nested-fields'));
});
