require("igv");
let $ = require('jquery');

let igvBrowser = null; //declaring this as a global for later
let tracks = [];
let primerSetsToNames = {};
let lineageSetsToNames = {};
let defaultCheckboxStatus = {};
let config = {};
let lineageSetNameMap = {};

function loadConfig()
{
    primerSetsToNames = JSON.parse($('#primer_set_json')[0].innerHTML);
    lineageSetsToNames = JSON.parse($('#lineage_set_json')[0].innerHTML);
    defaultCheckboxStatus = JSON.parse($('#default_tracks_json')[0].innerHTML);
    config = JSON.parse($('#config')[0].innerHTML);

    for (let lineageSetKey in lineageSetsToNames) {
        let lineageSet = lineageSetsToNames[lineageSetKey];
        lineageSetNameMap[lineageSet[0]] = lineageSet[1];
    }

}

function setSelectFormDisabled(state)
{
    $('.primer_set_checkbox').each(function (index, element) {
        element.disabled = state;
    });

    $('.lineage_set_radiobutton').each(function (index, element) {
        element.disabled = state;
    });

    $('#select_all_primers')[0].disabled = state;

    $('#unselect_all_primers')[0].disabled = state;

    $('#apply')[0].disabled = state;

    if(state)
    {
        $('#apply').addClass('is-loading');
    }
    else
    {
        $('#apply').removeClass('is-loading');
    }
}

function updatePrimerSets()
{
    const activeLineageGroup = lineageSetsToNames[$('input[name=lineage]:checked').val()][0];

    if(igvBrowser != null) { //if it's been loaded
        let activeSets = [];
        $('.primer_set_checkbox').each(function (index, element) {
            if (element.checked) {
                activeSets.push(element.name);
            }
        });
        loadPrimerSets(activeSets, igvBrowser, activeLineageGroup);
    }
}

function loadPrimerSets(activePrimerSets, igvBrowser, activeLineageGroup)
{
    //activeLineageGroup = "color"; //debug
    setSelectFormDisabled(true); //lock the form while changes are made
    tracks.forEach(function(track){
       igvBrowser.removeTrack(track);
    });

    tracks = [];

    const primerSetPromises = [];

    const variantsTrack = {
        "name": lineageSetNameMap[activeLineageGroup]+" Variants",
        "url": config['data_server']+"/"+config['organism_slug']+"/lineage_variants/"+encodeURIComponent(activeLineageGroup)+".bed",
        "format": "bed",
        "color": "#575757",
        "displayMode": "COLLAPSED",
        "autoHeight": true
    }
    igvBrowser.loadTrack(variantsTrack).then(function(addedTrack){
        tracks.push(addedTrack);

        activePrimerSets.forEach(function(primerSetKey){
            let primerSetData = primerSetsToNames[primerSetKey];
            const newTrack = {
                "name": primerSetData[1],
                "url": config['data_server']+"/"+config['organism_slug']+"/primer_sets_status/"+encodeURIComponent(primerSetData[0])+"/"+encodeURIComponent(activeLineageGroup)+".bed",
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

    });

}

function initBrowser() {

    const browserConfig =
        {
            reference: {
                "id": config['organism_slug'],
                "name": config['organism_name'],
                "fastaURL": config['data_server']+"/"+config['organism_slug']+"/ref/"+config['organism_slug']+".fasta",
                "indexURL": config['data_server']+"/"+config['organism_slug']+"/ref/"+config['organism_slug']+".fasta.fai",
                tracks: [
                    {
                        "name": "Genes",
                        "type": "annotation",
                        "url": config['data_server']+"/"+config['organism_slug']+"/ref/"+config['organism_slug']+".gff3",
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

    //remove duplicate browsers to prevent it from multiplying on back/forward
    $('.igv_div').children('.igv-container').remove();

    const browser_div = document.getElementById("igv");
    igv.createBrowser(browser_div, browserConfig).then(function (theBrowser) {
        igvBrowser = theBrowser;
        initCheckboxes();
        setRadiobuttons("XBB");
        updatePrimerSets();
    });

}

function setRadiobuttons(selected)
{
    $('.lineage_set_radiobutton').each(function (index, element) {
        element.checked = (lineageSetsToNames[element.value][0] === selected);
    });
}

function setCheckboxes(state)
{
    $('.primer_set_checkbox').each(function (index, element) {
        element.checked = state;
    });
}

function initCheckboxes()
{
    $('.primer_set_checkbox').each(function (index, element) {
        if(defaultCheckboxStatus[element.id])
        {
            element.checked = true;
        }
    });
}


$(document).ready(function(){
    $('#select_all_primers').on("click", function(){
        setCheckboxes(true);
    });

    $('#unselect_all_primers').on("click", function(){
        setCheckboxes(false);
    });

    $('#apply').on("click", function(){
        updatePrimerSets();
    });

    $('#primer_set_selection').on("submit", function(event){
        event.preventDefault(); //don't refresh
    });

    initCheckboxes();
    loadConfig();
    initBrowser();
})





