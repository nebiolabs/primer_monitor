import igv from "igv";

// The page's one igv.js browser, created and removed through igv's API so Turbo navigation never leaves a stale
// or duplicate browser behind. igv.createBrowser is async, so removal waits for a browser that is still loading.
let pending = null;

export function createBrowser(div, config) {
    removeBrowser();
    suppressTrackNamePopovers(div);
    pending = igv.createBrowser(div, config);
    return pending;
}

// igv.js opens a popover with each track's URL when its name is clicked, which we don't want. The click is stopped
// while it is still passing down the page (capture phase), so igv's own handler on the track name never sees it. The
// listener is added once per container, however often the page is visited.
function suppressTrackNamePopovers(div) {
    if (div.dataset.trackPopoversSuppressed) return;
    div.dataset.trackPopoversSuppressed = 'true';
    div.addEventListener('click', event => {
        if (event.target.closest('.igv-track-label')) event.stopPropagation();
    }, true);
}

export function removeBrowser() {
    if (!pending) return;
    const loading = pending;
    pending = null;
    loading.then(browser => igv.removeBrowser(browser)).catch(() => {});
}

// The promise of the page's browser, or null when there is none
export function currentBrowser() {
    return pending;
}

// The organism's reference genome with its gene annotations, for the diagrams on the organism's pages
export function referenceConfig(config) {
    const reference = config['reference_accession'];
    const slug = config['organism_slug'];
    const dataServer = config['data_server'];
    return {
        "id": reference,
        "name": config['organism_name'] + " (" + reference + ")",
        "fastaURL": dataServer + "/" + slug + "/ref/" + slug + ".fasta",
        "indexURL": dataServer + "/" + slug + "/ref/" + slug + ".fasta.fai",
        tracks: [
            {
                "name": "Genes",
                "type": "annotation",
                "url": dataServer + "/" + slug + "/ref/" + slug + ".gff3",
                "format": "gff3",
                "description": "",
                "filterTypes": ['CDS', 'mature_protein_region_of_CDS', 'region', 'stem_loop', 'five_prime_UTR', 'three_prime_UTR'],
                "displayMode": "EXPANDED",
                "colorBy": "gbkey",
                "colorTable": {
                    "Gene": "rgb(0,190,0)",
                }
            }
        ]
    };
}
