import 'leaflet'; // defines window.L

const TILE_URL = 'https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png';

// Circles at each point: bigger for more sequences, redder for a larger share of the location's sequences. Each
// location is {name, latitude, longitude, variant_count, sequence_count, fraction}; locations without coordinates are
// skipped. Returns the map and its markers by name.
export function drawPointMap(element, locations) {
    const placed = locations.filter(location => location.latitude != null && location.longitude != null);
    const map = window.L.map(element, { scrollWheelZoom: false }).setView([20, 0], 2);
    window.L.tileLayer(TILE_URL, {
        maxZoom: 8,
        attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
    }).addTo(map);

    const markersByName = new Map();
    const points = placed.map(location => [location.latitude, location.longitude]);
    const showAll = () => {
        if (points.length > 0) map.fitBounds(points, { padding: [30, 30], maxZoom: 5 });
    };
    const focus = marker => {
        map.flyTo(marker.getLatLng(), Math.max(map.getZoom(), 5));
        marker.openPopup();
    };

    placed.forEach(location => {
        const point = [location.latitude, location.longitude];
        const hue = Math.round((1 - location.fraction) * 60); // yellow (all sequences) to red (none)
        const marker = window.L.circleMarker(point, {
            radius: Math.min(5 + 3 * Math.sqrt(location.variant_count), 30),
            color: '#333',
            weight: 1,
            fillColor: `hsl(${hue}, 90%, 50%)`,
            fillOpacity: 0.8,
            bubblingMouseEvents: false // so clicking a circle does not also count as a background click
        }).addTo(map);
        marker.on('click', () => focus(marker));
        markersByName.set(location.name, marker);

        const popup = document.createElement('div');
        const name = document.createElement('strong');
        name.textContent = location.name;
        const detail = document.createElement('div');
        detail.textContent = `${location.variant_count} of ${location.sequence_count} sequences ` +
            `(${(location.fraction * 100).toFixed(1)}%)`;
        popup.append(name, detail);
        marker.bindPopup(popup);
    });

    map.on('click', showAll);
    showAll();
    return { map, markersByName, focus };
}

// Table rows with a data-map-point attribute zoom the map to that point when clicked
export function linkTableRows(element, { map, markersByName, focus }) {
    element.querySelectorAll('tr[data-map-point]').forEach(row => {
        const marker = markersByName.get(row.dataset.mapPoint);
        if (!marker) return;
        row.addEventListener('click', () => {
            map.getContainer().scrollIntoView({ behavior: 'smooth', block: 'center' });
            focus(marker);
        });
    });
}
