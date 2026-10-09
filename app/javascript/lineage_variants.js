import { createBrowser, removeBrowser } from 'igv_browser';
import { registerPageModule } from 'turbo_page_module';

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
        loadVariantTable(activeLineageGroup, activeSets);
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
            loadVariantTable(activeLineageGroup, activeSets);
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

function formatVariant(v) {
    const ref = v.ref;
    const alt = v.variant;
    if (v.variant_type === 'X') {
        return ref ? `${escapeHtml(ref)}>${escapeHtml(alt)}` : `>${escapeHtml(alt)}`;
    }
    if (v.variant_type === 'I') {
        return ref ? `${escapeHtml(ref)}>${escapeHtml(ref)}${escapeHtml(alt)}` : `>${escapeHtml(alt)}`;
    }
    if (v.variant_type === 'D') {
        return ref ? `\u0394${escapeHtml(ref)}` : `\u0394${escapeHtml(alt)}`;
    }
    return escapeHtml(alt);
}

function renderOligoSvg(sequence, oligoStart, oligoEnd, strand, variantStart, variantEnd) {
    const CELL_W = 12;
    const CELL_H = 20;
    const LABEL_W = 24;
    const ARROW_W = 8;
    const isPlus = strand !== '-';
    const n = sequence.length;
    const seqW = n * CELL_W;
    const totalW = LABEL_W + seqW + ARROW_W + LABEL_W;
    const totalH = CELL_H + 4;

    const fivePrimeX  = isPlus ? 0              : LABEL_W + seqW + ARROW_W;
    const threePrimeX = isPlus ? LABEL_W + seqW + ARROW_W : 0;

    const arrowX = isPlus ? LABEL_W + seqW : LABEL_W;
    const arrowPoints = isPlus
        ? `${arrowX},0 ${arrowX + ARROW_W},${CELL_H / 2} ${arrowX},${CELL_H}`
        : `${arrowX + ARROW_W},0 ${arrowX},${CELL_H / 2} ${arrowX + ARROW_W},${CELL_H}`;

    const cells = sequence.split('').map((base, i) => {
        const gPos = isPlus ? oligoStart + i : oligoEnd - 1 - i;
        const isVariant = gPos >= variantStart && gPos < variantEnd;
        const x = LABEL_W + i * CELL_W;
        const fill = isVariant ? '#cc0000' : 'rgb(0,0,200)';
        return `
            <rect x="${x}" y="0" width="${CELL_W}" height="${CELL_H}" fill="${fill}"/>
            <text x="${x + CELL_W / 2}" y="${CELL_H - 5}"
                  fill="white" font-size="10" text-anchor="middle"
                  font-family="monospace">${escapeHtml(base)}</text>`;
    }).join('');

    return `<svg xmlns="http://www.w3.org/2000/svg"
                 width="${totalW}" height="${totalH}"
                 style="display:block; cursor:pointer"
                 class="oligo-svg"
                 data-start="${oligoStart}" data-end="${oligoEnd}">
        <title>Show in IGV</title>
        <text x="${fivePrimeX + 2}" y="${CELL_H - 5}"
              fill="#666" font-size="10" font-family="monospace">5'</text>
        <text x="${threePrimeX + 2}" y="${CELL_H - 5}"
              fill="#666" font-size="10" font-family="monospace">3'</text>
        <polygon points="${arrowPoints}" fill="rgb(0,0,200)"/>
        ${cells}
    </svg>`;
}

function escapeHtml(str) {
    return String(str ?? '').replace(/[&<>"']/g, c =>
        ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])
    );
}

function formatSeen(seen) {
    if (!seen || !seen.date) return '—';
    const parts = [seen.date, seen.lineage].filter(Boolean);
    const line = parts.map(escapeHtml).join(' · ');
    const loc = seen.location ? `<br><small class="has-text-grey">${escapeHtml(seen.location)}</small>` : '';
    return line + loc;
}

function buildVariantTable(variants) {
    if (variants.length === 0) {
        return '<p class="has-text-grey">No variants overlap the selected primers for this lineage.</p>';
    }

    const rows = variants.map(v => {
        const oligoSubRows = v.oligos.map(o => `
            <tr class="variant-oligo-row" style="display:none">
                <td colspan="3" style="padding-left:2rem; background:#f9f9f9">
                    <strong>${escapeHtml(o.name)}</strong>
                    <span class="has-text-grey"> — ${escapeHtml(o.primer_set)}</span>
                    <div style="overflow-x:auto; margin-top:0.4rem">
                        ${renderOligoSvg(o.sequence, o.oligo_start, o.oligo_end, o.strand, v.ref_start, v.ref_end)}
                    </div>
                </td>
                <td style="background:#f9f9f9; vertical-align:top">${formatSeen(v.first_seen)}</td>
                <td style="background:#f9f9f9; vertical-align:top">${formatSeen(v.last_seen)}</td>
            </tr>
        `).join('');

        return `
            <tr class="variant-row is-clickable" data-expanded="false">
                <td>
                    <button class="variant-expand-btn button is-small is-white" aria-label="expand">▶</button>
                    ${v.ref_start}–${v.ref_end}
                </td>
                <td>${formatVariant(v)}</td>
                <td>${v.frequency_pct.toFixed(1)}%</td>
                <td>${formatSeen(v.first_seen)}</td>
                <td>${formatSeen(v.last_seen)}</td>
            </tr>
            ${oligoSubRows}
        `;
    }).join('');

    return `
        <table class="table is-fullwidth is-hoverable is-narrow">
            <thead>
                <tr>
                    <th>Position</th>
                    <th>Variant</th>
                    <th>Frequency</th>
                    <th>First Seen</th>
                    <th>Last Seen</th>
                </tr>
            </thead>
            <tbody>${rows}</tbody>
        </table>
    `;
}

function loadVariantTable(lineage, primerSets) {
    const section = document.getElementById('variant_table_section');
    const loading = document.getElementById('variant_table_loading');
    const content = document.getElementById('variant_table_content');
    if (!section) return;

    if (!lineage || primerSets.length === 0) {
        section.style.display = 'none';
        return;
    }

    const lineageLabel = lineageSetsToNames[lineage] || lineage;
    const heading = document.getElementById('variant_table_heading');
    if (heading) heading.textContent = `Variants in ${lineageLabel} Overlapping Selected Primers`;

    section.style.display = 'block';
    loading.style.display = 'block';
    content.innerHTML = '';

    const params = new URLSearchParams();
    params.set('lineage', lineage);
    // primerSets contains URL-safe keys from tracks.json; send display names to match primer_sets.name in DB
    primerSets.forEach(ps => params.append('primer_sets[]', primerSetsToNames[ps] || ps));

    const url = `${location.pathname}/variant_overlaps.json?${params}`;

    fetch(url)
        .then(r => r.json())
        .then(data => {
            loading.style.display = 'none';
            content.innerHTML = buildVariantTable(data.variants || []);
        })
        .catch(() => {
            loading.style.display = 'none';
            content.innerHTML = '<p class="has-text-danger">Error loading variant data.</p>';
        });
}

document.addEventListener('click', function(event) {
    const oligo = event.target.closest('.oligo-svg');
    if (!oligo || !igvBrowser || !config['reference_accession']) return;
    const start = Math.max(0, parseInt(oligo.dataset.start) - 50);
    const end = parseInt(oligo.dataset.end) + 50;
    igvBrowser.search(`${config['reference_accession']}:${start}-${end}`)
        .then(() => document.getElementById('igv').scrollIntoView({ behavior: 'smooth', block: 'start' }));
});

document.addEventListener('click', function(event) {
    const row = event.target.closest('.variant-row');
    if (!row) return;
    const expanded = row.dataset.expanded === 'true';
    row.dataset.expanded = String(!expanded);
    row.querySelector('.variant-expand-btn').textContent = expanded ? '▶' : '▼';
    for (let next = row.nextElementSibling; next?.classList.contains('variant-oligo-row'); next = next.nextElementSibling) {
        next.style.display = expanded ? 'none' : '';
    }
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
