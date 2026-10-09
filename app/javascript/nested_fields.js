// Add/remove rows of a Rails nested form (accepts_nested_attributes_for), replacing cocoon.
//
//   <div data-nested-fields>
//     <div class="nested-fields">... f.hidden_field :_destroy ... <button data-nested-remove></button></div>
//     <template data-nested-template> fields_for ... child_index: 'NEW_RECORD' </template>
//     <button type="button" data-nested-add>Add</button>
//   </div>
//
// Each added row gets a unique index in place of NEW_RECORD. Removing a saved row (data-persisted="true") hides it
// and sets _destroy; removing a row that was never saved just deletes it.

let counter = 0;

export function addNestedRow(container) {
    const template = container.querySelector('template[data-nested-template]');
    const html = template.innerHTML.replace(/NEW_RECORD/g, `${Date.now()}${counter++}`);
    template.insertAdjacentHTML('beforebegin', html);
    const rows = container.querySelectorAll('.nested-fields');
    return rows[rows.length - 1];
}

document.addEventListener('click', (event) => {
    const add = event.target.closest('[data-nested-add]');
    if (add) {
        event.preventDefault();
        addNestedRow(add.closest('[data-nested-fields]'));
        return;
    }

    const remove = event.target.closest('[data-nested-remove]');
    if (!remove) return;
    event.preventDefault();
    const row = remove.closest('.nested-fields');
    if (row.dataset.persisted === 'true') {
        row.querySelector('input[name$="[_destroy]"]').value = '1';
        row.style.display = 'none';
    } else {
        row.remove();
    }
});
