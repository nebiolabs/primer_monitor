import { createBrowser, removeBrowser } from 'igv_browser';
import { registerPageModule } from 'turbo_page_module';
import { loadVariantTable } from 'variant_table';

let config = {};
let igvBrowser = null;
let lineageTracks = [];
let lineageRequest = 0;

function primerSetTrack() {
    return {
        "name": config['primer_set_display_name'],
        "url": config['data_server']+"/"+config['organism_slug']+"/primer_sets_bed/"+config['primer_set_name']+".bed",
        "format": "bed",
        "displayMode": "EXPANDED",
        "autoHeight": true
    };
}

// The lineage's variants, then this primer set colored by whether those variants likely impact each primer
function lineageTrackConfigs(lineage) {
    const base = config['data_server'] + "/" + config['organism_slug'];
    return [
        {
            "name": (config['lineage_names'][lineage] || lineage) + " Variants",
            "url": base + "/lineage_variants/" + encodeURIComponent(lineage) + ".bed",
            "format": "bed",
            "color": "#575757",
            "displayMode": "COLLAPSED",
            "autoHeight": true
        },
        {
            "name": config['primer_set_display_name'],
            "url": base + "/primer_sets_status/" + encodeURIComponent(config['primer_set_name']) + "/" + encodeURIComponent(lineage) + ".bed",
            "format": "bed",
            "displayMode": "EXPANDED",
            "autoHeight": true
        }
    ];
}

async function showLineage(lineage) {
    const request = ++lineageRequest;
    loadVariantTable(
        config['variant_overlaps_path'],
        lineage,
        [config['primer_set_display_name']],
        `Variants in ${config['lineage_names'][lineage] || lineage} Overlapping This Primer Set`,
        config['reference_accession']
    );
    lineageTracks.forEach(track => igvBrowser.removeTrack(track));
    lineageTracks = [];
    // load in order so the variants stay above the primers
    for (const trackConfig of lineageTrackConfigs(lineage)) {
        const track = await igvBrowser.loadTrack(trackConfig);
        // a newer selection arrived while this track loaded
        if (request !== lineageRequest) {
            igvBrowser?.removeTrack(track);
            return;
        }
        lineageTracks.push(track);
    }
}

function initBrowser() {

    const browserConfig =
        {
            reference: {
                "id": config['reference_accession'],
                "name": config['organism_name']+" ("+config['reference_accession']+")",
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
                    },
                    // without lineage data, show just where the primers bind
                    ...(config['initial_lineage'] ? [] : [primerSetTrack()])
                ]
            }
        };

    createBrowser(document.getElementById("igv"), browserConfig).then(function(theBrowser) {
        igvBrowser = theBrowser;
        if (config['initial_lineage']) showLineage(config['initial_lineage']);
    });

}

// Document-level delegation so the handlers are registered once and don't stack across Turbo navigations
document.addEventListener('change', function(event) {
    if (!event.target.matches('#primer_set_lineage_select') || !igvBrowser) return;
    const lineage = event.target.value;
    // keep the selection in the URL so it can be shared or bookmarked
    const url = new URL(location.href);
    url.searchParams.set('lineage', lineage);
    history.replaceState(history.state, '', url);
    showLineage(lineage);
});

document.addEventListener('submit', function(event) {
    if (event.target.closest('#primer_set_lineage_selection')) event.preventDefault();
});

registerPageModule(
    () => !!document.querySelector('#config[data-igv-page="primer_set"]'),
    () => {
        config = JSON.parse(document.getElementById('config').innerHTML);
        if ('primer_set_name' in config) initBrowser();
    },
    () => {
        removeBrowser();
        igvBrowser = null;
        lineageTracks = [];
        lineageRequest++;
    }
);
