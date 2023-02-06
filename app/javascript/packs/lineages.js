require("igv");
let $ = require('jquery');

let igvBrowser = null; //declaring this as a global for later
let tracks = [];
let primerSetsToNames = {};


function setupNameMapping()
{
    /*$('.primer_set_checkbox_label').each(function (index, element) {
       const id = element.getAttribute("for");
       let name = "";
       const children = element.children;

       for(let i=0; i<children.length; i++)
       {
          if(children[i].classList.contains("primer_display_name"))
          {
              name = children[i].innerHTML;
          }
       }
       primerSetsToNames[id] = name;
    });*/
    primerSetsToNames = JSON.parse($('#primer_set_json')[0].innerHTML);


}

function setSelectFormDisabled(state)
{
    $('.primer_set_checkbox').each(function (index, element) {
        element.disabled = state;
    });

    $('#select_all_primers')[0].disabled = state;

    $('#unselect_all_primers')[0].disabled = state;
}

function updatePrimerSets()
{
    if(igvBrowser != null) { //if it's been loaded
        let activeSets = [];
        $('.primer_set_checkbox').each(function (index, element) {
            if (element.checked) {
                activeSets.push(element.name);
            }
        });
        loadPrimerSets(activeSets, igvBrowser);
    }
}

async function loadPrimerSets(activePrimerSets, igvBrowser)
{
    setSelectFormDisabled(true); //lock the form while changes are made
    tracks.forEach(function(track){
       igvBrowser.removeTrack(track);
    });

    tracks = [];

    const primerSetPromises = [];

    activePrimerSets.forEach(function(primerSetKey){
        let primerSetData = primerSetsToNames[primerSetKey];
        const newTrack = {
            "name": primerSetData[1],
            "url": "http://localhost:8080/primer_sets/"+encodeURIComponent(primerSetData[0])+"/color.bed",
            "format": "bed",
            "displayMode": "EXPANDED",
            "autoHeight": true
        }
        primerSetPromises.push(igvBrowser.loadTrack(newTrack));
    });

    Promise.all(primerSetPromises).then(function(addedTracks){
        addedTracks.forEach(function(addedTrack){
            tracks.push(addedTrack);
        });
        setSelectFormDisabled(false); //unlock the form
    }).catch(function(){
        setSelectFormDisabled(false); //also unlock the form
    });


    //.then(function(addedTrack){
    //             tracks[primerSetKey] = addedTrack;
    //         });
    //await (primerSetsToDo == 0); //wait for loading to finish
    //alert("done");

}

function initBrowser() {

    const config =
        {
            reference: {
                "id": "NC_045512.2",
                "name": "NC_045512.2 (SARS-CoV-2)",
                "fastaURL": "http://localhost:8080/ref/NC_045512.2.fasta",
                "indexURL": "http://localhost:8080/ref/NC_045512.2.fasta.fai",
                tracks: [
                    {
                        "name": "Genes",
                        "type": "annotation",
                        "url": "http://localhost:8080/ref/NC_045512.2.gff3",
                        "format": "gff3",
                        "filterTypes": ['CDS', 'mature_protein_region_of_CDS', 'region', 'stem_loop', 'five_prime_UTR', 'three_prime_UTR'],
                        "displayMode": "EXPANDED",
                        "colorBy": "gbkey",
                        "colorTable": {
                            "Gene": "rgb(0,190,0)",
                        }
                    }
                ]
            }
        };

    const browser_div = document.getElementById("igv");
    igv.createBrowser(browser_div, config).then(function (theBrowser) {
        igvBrowser = theBrowser;
        setupNameMapping();
        loadPrimerSets(Object.keys(primerSetsToNames), igvBrowser);
    });

}

function setCheckboxes(state)
{
    $('.primer_set_checkbox').each(function (index, element) {
        element.checked = state;
    });
}


$(document).ready(function(){
    $('.primer_set_checkbox').on("click", function(){
        updatePrimerSets();
    });

    $('#select_all_primers').on("click", function(){
        setCheckboxes(true);
        updatePrimerSets();
    });

    $('#unselect_all_primers').on("click", function(){
        setCheckboxes(false);
        updatePrimerSets();
    });

    $('#primer_set_selection').on("submit", function(event){
        event.preventDefault(); //don't refresh
    });

    initBrowser();
})





