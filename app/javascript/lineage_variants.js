import { createBrowser, removeBrowser } from 'igv_browser';
import { registerPageModule } from 'turbo_page_module';
import { loadVariantTable } from 'variant_table';

let igvBrowser = null;
let tracks = [];
let primerSetsToNames = {};
let lineageSetsToNames = {};
let config = {};

let activeLineageGroup = null;
let activeSets = [];

function loadConfig() {
    primerSetsToNames = JSON.parse(document.getElementById('primer_set_json').textContent);
    lineageSetsToNames = JSON.parse(document.getElementById('lineage_set_json').textContent);
    config = JSON.parse(document.getElementById('config').textContent);
}


function updateLink() {
    const wrapper = document.getElementById('link_div_wrapper');
    const linkElement = document.getElementById('link');
    if (wrapper.classList.contains('invisible')) {
        const baseLink = location.protocol + '//' + location.host + location.pathname;
        const fullLink = baseLink + "?primer_sets=" + activeSets.join(',') + ";lineage=" + activeLineageGroup;
        linkElement.textContent = fullLink;
        linkElement.href = fullLink;
        wrapper.classList.remove('invisible');
        document.getElementById('show_link').textContent = 'Hide Link';
    } else {
        wrapper.classList.add('invisible');
        document.getElementById('show_link').textContent = 'Shareable Link';
    }
}

function updatePrimerSets() {
    document.getElementById('link_div_wrapper').classList.add('invisible');
    document.getElementById('show_link').textContent = 'Shareable Link';
    const linkElement = document.getElementById('link');
    linkElement.textContent = "";
    linkElement.href = "";

    activeLineageGroup = document.getElementById('lineage_select').value;

    if (igvBrowser != null) {
        activeSets = [...document.getElementById('primer_set_select').selectedOptions].map(option => option.value);
        loadPrimerSets(activeSets, igvBrowser, activeLineageGroup);
        showVariantTable(activeLineageGroup, activeSets);
    }
}

function loadPrimerSets(activePrimerSets, igvBrowser, activeLineageGroup) {
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
        });
    });
}

function showVariantTable(lineage, primerSets) {
    const lineageLabel = lineageSetsToNames[lineage] || lineage;
    // primerSets contains URL-safe keys from tracks.json; send display names to match primer_sets.name in DB
    loadVariantTable(
        `${location.pathname}/variant_overlaps.json`,
        lineage,
        primerSets.map(ps => primerSetsToNames[ps] || ps),
        `Variants in ${lineageLabel} Overlapping Selected Primers`,
        config['reference_accession']
    );
}

function initBrowser() {
    const browserConfig = {
        reference: {
            "id": config['reference_accession'],
            "name": config['organism_name'] + " (" + config['reference_accession'] + ")",
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

    createBrowser(document.getElementById("igv"), browserConfig).then(function(theBrowser) {
        igvBrowser = theBrowser;
        document.getElementById('igv_loading').classList.add('invisible');
        activeLineageGroup = config['initial_lineage'];
        activeSets = config['initial_primer_sets'] || [];
        if (activeLineageGroup) {
            loadPrimerSets(activeSets, igvBrowser, activeLineageGroup);
            showVariantTable(activeLineageGroup, activeSets);
        }
    });
}

// Use document-level delegation so handlers are registered once and don't stack
// across Turbo navigations.
document.addEventListener('click', function(event) {
    if (event.target.closest('#show_link')) updateLink();
});

document.addEventListener('submit', function(event) {
    if (event.target.closest('#primer_set_selection')) event.preventDefault();
});

let updateTimer = null;
function debouncedUpdate() {
    clearTimeout(updateTimer);
    updateTimer = setTimeout(updatePrimerSets, 250);
}

document.addEventListener('change', function(event) {
    if (event.target.matches('#lineage_select, #primer_set_select')) debouncedUpdate();
});

registerPageModule(
    () => !!document.getElementById('lineage_select'),
    () => { loadConfig(); initBrowser(); },
    () => {
        removeBrowser();
        document.getElementById('igv_loading').classList.remove('invisible');
        igvBrowser = null;
        tracks = [];
        activeSets = [];
        activeLineageGroup = null;
    }
);
