import { createBrowser, removeBrowser, referenceConfig, trackDescription } from 'igv_browser';
import { drawPointMap, linkTableRows } from 'point_map';
import { registerPageModule } from 'turbo_page_module';

let mapHandle = null;

// Diagram of the reference with one track per primer set, each coloured by the impact of this lineage's variants
function initDiagram(config) {
    createBrowser(document.getElementById('lineage_igv'), {
        reference: referenceConfig(config),
        tracks: (config['primer_set_tracks'] || []).map(primerSet => ({
            "name": primerSet.name,
            "url": primerSet.url,
            "format": "bed",
            "color": primerSet.color,
            "displayMode": "SQUISHED",
            "height": 17,
            "autoHeight": false,
            ...trackDescription(primerSet.description)
        }))
    });
}

function initMap(element) {
    mapHandle = drawPointMap(element, JSON.parse(element.dataset.points));
    linkTableRows(document, mapHandle);
}

function cleanup() {
    removeBrowser();
    if (mapHandle) mapHandle.map.remove();
    mapHandle = null;
}

registerPageModule(
    () => !!document.getElementById('lineage_config'),
    () => {
        initDiagram(JSON.parse(document.getElementById('lineage_config').textContent));
        const mapElement = document.getElementById('lineage_map');
        if (mapElement) initMap(mapElement);
    },
    cleanup
);
