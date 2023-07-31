const $ = require("jquery");



$(document).ready(function(){

    const fasta_upload = $('#fasta_upload');

    fasta_upload.on('change', function(){
        if(fasta_upload[0].files.length > 0)
        {
            if($('#fasta_name').length === 0)
            {
                const fastaName = document.createElement("span");
                fastaName.id = 'fasta_name';
                fastaName.className = 'file-name';
                $('#fasta_upload_label')[0].appendChild(fastaName);
                $('#fasta_div').addClass('has-name');
            }
            $('#fasta_name')[0].innerHTML = fasta_upload[0].files[0].name;
            fasta_upload[0].files[0].text().then(processFasta);
        }
    });


});

function addPrimer(lastPrimer)
{
    let longName = lastPrimer[0];
    let shortName = lastPrimer[1];
    let sequence = lastPrimer[2];
    let newPrimerContainer = $('#samples > div:nth-last-of-type(2) > div.field-body')[0];
    for(const child of newPrimerContainer.children)
    {
        if(child.tagName.toUpperCase() === "INPUT")
        {
            if(child.id.endsWith("short_name"))
            {
                child.value = shortName;
            }
            else if(child.id.endsWith("name"))
            {
                child.value = longName;
            }
            else if(child.id.endsWith("sequence"))
            {
                child.value = sequence;
            }
        }
    }
}

function processFasta(fastaText)
{
    const fastaName = $('#fasta_name')[0].innerHTML;
    if(fastaText[0] !== ">")
    {
        alert("'"+fastaName+"' is not a valid FASTA file!");
        return;
    }
    fastaText.slice(1).split("\n>").forEach(function(seq){
        let seqDataRaw = seq.split("\n");
        let longName = seqDataRaw[0];
        let sequence = seqDataRaw.slice(1).join("");
        let shortName = longName.split(" ")[0].slice(0,5);
        const invalidSeqRE = /[^autcgnbvdhrykmws\- \t]/i; //matches everything except for all bases and ambiguity codes, - (for gaps), and whitespace
        if(sequence.match(invalidSeqRE) !== null)
        {
            confirm("Non-sequence character detected in file '"+fastaName+"'.\n\nSequence:\n>"+seq+"\n\nContinue?");
        }
        let createOligo = $('.add_fields')[0];
        createOligo.click.call(createOligo);
        addPrimer([longName, shortName, sequence]);
    })
}