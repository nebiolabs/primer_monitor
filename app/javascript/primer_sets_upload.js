import { addNestedRow } from 'nested_fields';

// Fills the primer set form's oligo rows from an uploaded FASTA file, one row per sequence.

document.addEventListener('change', (event) => {
    const upload = event.target.closest('#fasta_upload');
    if (!upload || upload.files.length === 0) return;

    let fastaName = document.getElementById('fasta_name');
    if (!fastaName) {
        fastaName = document.createElement('span');
        fastaName.id = 'fasta_name';
        fastaName.className = 'file-name';
        document.getElementById('fasta_upload_label').appendChild(fastaName);
        document.getElementById('fasta_div').classList.add('has-name');
    }
    fastaName.textContent = upload.files[0].name;
    upload.files[0].text().then(text => processFasta(text, upload.files[0].name));
});

function fillRow(row, longName, shortName, sequence) {
    row.querySelectorAll('input').forEach(input => {
        if (input.id.endsWith('short_name')) input.value = shortName;
        else if (input.id.endsWith('name')) input.value = longName;
        else if (input.id.endsWith('sequence')) input.value = sequence;
    });
}

function processFasta(fastaText, fileName) {
    if (fastaText[0] !== '>') {
        alert(`'${fileName}' is not a valid FASTA file!`);
        return;
    }
    const container = document.querySelector('#samples[data-nested-fields]');
    // matches everything except bases, ambiguity codes, - (for gaps), and whitespace
    const invalidSeqRE = /[^autcgnbvdhrykmws\- \t\r]/i;
    fastaText.slice(1).split('\n>').forEach(seq => {
        const lines = seq.split('\n');
        const longName = lines[0];
        const sequence = lines.slice(1).join('');
        const shortName = longName.split(' ')[0].slice(0, 5);
        if (invalidSeqRE.test(sequence)) {
            confirm(`Non-sequence character detected in file '${fileName}'.\n\nSequence:\n>${seq}\n\nContinue?`);
        }
        fillRow(addNestedRow(container), longName, shortName, sequence);
    });
}
