import 'leaflet'; // defines window.L
import { createBrowser, removeBrowser } from 'igv_browser';
import { registerPageModule } from 'turbo_page_module';

const TILE_URL = 'https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png';

let leafletMap = null;

// Diagram of the reference around the position: genes, then the variant and the primers overlapping it
function initDiagram(config) {
    const reference = config['reference_accession'];
    const slug = config['organism_slug'];
    const dataServer = config['data_server'];
    const browserConfig = {
        reference: {
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
                    "filterTypes": ['CDS', 'mature_protein_region_of_CDS', 'region', 'stem_loop', 'five_prime_UTR', 'three_prime_UTR'],
                    "displayMode": "EXPANDED",
                    "colorBy": "gbkey",
                    "colorTable": {
                        "Gene": "rgb(0,190,0)",
                    }
                }
            ]
        },
        locus: config['locus'],
        tracks: [
            {
                "name": "Primers overlapping position",
                "url": config['primer_track_url'],
                "format": "bed",
                "displayMode": "SQUISHED",
                "autoHeight": true
            },
            // one track per allele, most frequent first; the name carries its frequency
            ...(config['variant_tracks'] || []).map(variant => ({
                "name": variant.name,
                "url": variant.url,
                "format": "bed",
                "color": "#d62828",
                "displayMode": "SQUISHED",
                "height": 12,
                "autoHeight": false
            }))
        ]
    };

    createBrowser(document.getElementById('variant_igv'), browserConfig);
}

// Circles at each location's coordinates: bigger for more sequences, redder for a larger share of them
function initMap(element) {
    const locations = JSON.parse(element.dataset.locations)
        .filter(location => location.latitude != null && location.longitude != null);

    leafletMap = window.L.map(element, { scrollWheelZoom: false }).setView([20, 0], 2);
    window.L.tileLayer(TILE_URL, {
        maxZoom: 8,
        attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
    }).addTo(leafletMap);

    const points = locations.map(location => [location.latitude, location.longitude]);
    const showAll = () => {
        if (points.length > 0) leafletMap.fitBounds(points, { padding: [30, 30], maxZoom: 5 });
    };

    locations.forEach(location => {
        const point = [location.latitude, location.longitude];
        const hue = Math.round((1 - location.fraction) * 60); // yellow (all sequences) to red (none)
        const marker = window.L.circleMarker(point, {
            radius: Math.min(5 + 3 * Math.sqrt(location.variant_count), 30),
            color: '#333',
            weight: 1,
            fillColor: `hsl(${hue}, 90%, 50%)`,
            fillOpacity: 0.8,
            bubblingMouseEvents: false // so clicking a circle does not also count as a background click
        }).addTo(leafletMap);
        marker.on('click', () => leafletMap.flyTo(point, Math.max(leafletMap.getZoom(), 5)));

        const popup = document.createElement('div');
        const name = document.createElement('strong');
        name.textContent = location.name;
        const detail = document.createElement('div');
        detail.textContent = `${location.variant_count} of ${location.sequence_count} sequences ` +
            `(${(location.fraction * 100).toFixed(1)}%)`;
        popup.append(name, detail);
        marker.bindPopup(popup);
    });

    leafletMap.on('click', showAll);
    showAll();
}

function cleanup() {
    removeBrowser();
    if (leafletMap) leafletMap.remove();
    leafletMap = null;
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
