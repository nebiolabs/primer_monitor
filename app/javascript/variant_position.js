import { createBrowser, removeBrowser, referenceConfig } from 'igv_browser';
import { drawPointMap, linkTableRows } from 'point_map';
import { registerPageModule } from 'turbo_page_module';

let mapHandle = null;

// Diagram of the reference around the position: genes, then the primers and the variants overlapping it
function initDiagram(config) {
    const browserConfig = {
        reference: referenceConfig(config),
        locus: config['locus'],
        tracks: [
            // one track per primer, most impacted first, coloured by impact
            ...(config['primer_tracks'] || []).map(primer => ({
                "name": primer.name,
                "url": primer.url,
                "format": "bed",
                "color": primer.color,
                "displayMode": "SQUISHED",
                "height": 17,
                "autoHeight": false
            })),
            // one track per allele, most frequent first; the name carries its frequency
            ...(config['variant_tracks'] || []).map(variant => ({
                "name": variant.name,
                "url": variant.url,
                "format": "bed",
                "color": variant.color || "#d62828",
                "displayMode": "SQUISHED",
                "height": 12,
                "autoHeight": false
            }))
        ]
    };

    createBrowser(document.getElementById('variant_igv'), browserConfig);
}

function initMap(element) {
    const locations = JSON.parse(element.dataset.locations);
    mapHandle = drawPointMap(element, locations);
    linkTableRows(document, mapHandle);
}

function cleanup() {
    removeBrowser();
    if (mapHandle) mapHandle.map.remove();
    mapHandle = null;
}

registerPageModule(
    () => !!document.getElementById('variant_position_config'),
    () => {
        initDiagram(JSON.parse(document.getElementById('variant_position_config').textContent));
        const mapElement = document.getElementById('variant_map');
        if (mapElement) initMap(mapElement);
    },
    cleanup
);
