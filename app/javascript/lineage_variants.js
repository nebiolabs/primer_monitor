import 'init_jquery';
import "igv";

let igvBrowser = null;
let tracks = [];
let primerSetsToNames = {};
let lineageSetsToNames = {};
let config = {};

let activeLineageGroup = null;
let activeSets = [];

function loadConfig() {
    primerSetsToNames = JSON.parse($('#primer_set_json')[0].innerHTML);
    lineageSetsToNames = JSON.parse($('#lineage_set_json')[0].innerHTML);
    config = JSON.parse($('#config')[0].innerHTML);
}

function setSelectFormDisabled(state) {
    $('#apply').prop('disabled', state);
    if (state) {
        $('#apply').addClass('is-loading');
    } else {
        $('#apply').removeClass('is-loading');
    }
}

function updateLink() {
    let link_div_wrapper = $('#link_div_wrapper');
    if (link_div_wrapper.hasClass('invisible')) {
        let base_link = location.protocol + '//' + location.host + location.pathname;
        let full_link = base_link + "?primer_sets=" + activeSets.join(',') + ";lineage=" + activeLineageGroup;
        let link_element = $('#link')[0];
        link_element.innerHTML = full_link;
        link_element.href = full_link;
        link_div_wrapper.removeClass('invisible');
        $('#show_link')[0].innerHTML = 'Hide Link';
    } else {
        link_div_wrapper.addClass('invisible');
        $('#show_link')[0].innerHTML = 'Shareable Link';
    }
}

function updatePrimerSets() {
    $('#link_div_wrapper').addClass('invisible');
    $('#show_link')[0].innerHTML = 'Shareable Link';
    let link_element = $('#link')[0];
    link_element.innerHTML = "";
    link_element.href = "";

    activeLineageGroup = $('#lineage_select').val();

    if (igvBrowser != null) {
        activeSets = $('#primer_set_select').val() || [];
        loadPrimerSets(activeSets, igvBrowser, activeLineageGroup);
    }
}

function loadPrimerSets(activePrimerSets, igvBrowser, activeLineageGroup) {
    setSelectFormDisabled(true);
    tracks.forEach(function(track) {
        igvBrowser.removeTrack(track);
    });
    tracks = [];

    const primerSetPromises = [];

    const variantsTrack = {
        "name": (lineageSetsToNames[activeLineageGroup] || activeLineageGroup) + " Variants",
        "url": config['data_server'] + "/" + config['organism_slug'] + "/lineage_variants/" + encodeURIComponent(activeLineageGroup) + ".bed",
        "format": "bed",
        "color": "#575757",
        "displayMode": "COLLAPSED",
        "autoHeight": true
    };

    igvBrowser.loadTrack(variantsTrack).then(function(addedTrack) {
        tracks.push(addedTrack);

        activePrimerSets.forEach(function(primerKey) {
            const newTrack = {
                "name": primerSetsToNames[primerKey] || primerKey,
                "url": config['data_server'] + "/" + config['organism_slug'] + "/primer_sets_status/" + encodeURIComponent(primerKey) + "/" + encodeURIComponent(activeLineageGroup) + ".bed",
                "format": "bed",
                "displayMode": "EXPANDED",
                "autoHeight": true
            };
            primerSetPromises.push(igvBrowser.loadTrack(newTrack));
        });

        Promise.all(primerSetPromises).then(function(addedTracks) {
            addedTracks.forEach(function(addedTrack) {
                tracks.push(addedTrack);
            });
            setSelectFormDisabled(false);
        }).catch(function() {
            setSelectFormDisabled(false);
        });
    });
}

function initBrowser() {
    const browserConfig = {
        reference: {
            "id": config['organism_slug'],
            "name": config['organism_name'],
            "fastaURL": config['data_server'] + "/" + config['organism_slug'] + "/ref/" + config['organism_slug'] + ".fasta",
            "indexURL": config['data_server'] + "/" + config['organism_slug'] + "/ref/" + config['organism_slug'] + ".fasta.fai",
            tracks: [
                {
                    "name": "Genes",
                    "type": "annotation",
                    "url": config['data_server'] + "/" + config['organism_slug'] + "/ref/" + config['organism_slug'] + ".gff3",
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

    $('.igv_div').children('.igv-container').remove();

    const browser_div = document.getElementById("igv");
    igv.createBrowser(browser_div, browserConfig).then(function(theBrowser) {
        igvBrowser = theBrowser;
        updatePrimerSets();
    });
}

$(document).ready(function() {
    $('#apply').on("click", function() {
        updatePrimerSets();
    });

    $('#show_link').on("click", function() {
        updateLink();
    });

    $('#primer_set_selection').on("submit", function(event) {
        event.preventDefault();
    });

    loadConfig();
    initBrowser();
});
