import 'init_jquery';
import "igv";
import { registerPageModule } from 'turbo_page_module';

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
        $('#igv_loading').addClass('invisible');
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
$(document).on('click', '#show_link', function() {
    updateLink();
});

$(document).on('submit', '#primer_set_selection', function(event) {
    event.preventDefault();
});

let updateTimer = null;
function debouncedUpdate() {
    clearTimeout(updateTimer);
    updateTimer = setTimeout(updatePrimerSets, 250);
}

$(document).on('change', '#lineage_select', debouncedUpdate);
$(document).on('change', '#primer_set_select', debouncedUpdate);

const VARIANT_TYPE_LABELS = { X: 'SNP', D: 'deletion', I: 'insertion' };

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
    const mid = CELL_H / 2;
    const arrowPoints = isPlus
        ? `${arrowX},${mid - 4} ${arrowX + ARROW_W},${mid} ${arrowX},${mid + 4}`
        : `${arrowX + ARROW_W},${mid - 4} ${arrowX},${mid} ${arrowX + ARROW_W},${mid + 4}`;

    const cells = sequence.split('').map((base, i) => {
        const gPos = isPlus ? oligoStart + i : oligoEnd - 1 - i;
        const isVariant = gPos >= variantStart && gPos < variantEnd;
        const x = LABEL_W + i * CELL_W;
        const fill = isVariant ? '#cc0000' : '#444444';
        return `
            <rect x="${x}" y="0" width="${CELL_W}" height="${CELL_H}" fill="${fill}"/>
            <text x="${x + CELL_W / 2}" y="${CELL_H - 5}"
                  fill="white" font-size="10" text-anchor="middle"
                  font-family="monospace">${escapeHtml(base)}</text>`;
    }).join('');

    return `<svg xmlns="http://www.w3.org/2000/svg"
                 width="${totalW}" height="${totalH}"
                 style="display:block">
        <text x="${fivePrimeX + 2}" y="${CELL_H - 5}"
              fill="#666" font-size="10" font-family="monospace">5'</text>
        <text x="${threePrimeX + 2}" y="${CELL_H - 5}"
              fill="#666" font-size="10" font-family="monospace">3'</text>
        <polygon points="${arrowPoints}" fill="#888"/>
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
    const parts = [seen.date, seen.lineage, seen.location].filter(Boolean);
    return parts.map(escapeHtml).join(' · ');
}

function buildVariantTable(variants) {
    if (variants.length === 0) {
        return '<p class="has-text-grey">No variants overlap the selected primers for this lineage.</p>';
    }

    const rows = variants.map(v => {
        const oligoSubRows = v.oligos.map(o => `
            <tr class="variant-oligo-row" style="display:none">
                <td colspan="6" style="padding-left:2rem; background:#f9f9f9">
                    <strong>${escapeHtml(o.name)}</strong>
                    <span class="has-text-grey"> — ${escapeHtml(o.primer_set)}</span><br>
                    <div style="overflow-x:auto; margin-top:0.4rem">
                        ${renderOligoSvg(o.sequence, o.oligo_start, o.oligo_end, o.strand, v.ref_start, v.ref_end)}
                    </div>
                </td>
            </tr>
        `).join('');

        return `
            <tr class="variant-row is-clickable" data-expanded="false">
                <td>
                    <button class="variant-expand-btn button is-small is-white" aria-label="expand">▶</button>
                    ${v.ref_start}–${v.ref_end}
                </td>
                <td>${escapeHtml(VARIANT_TYPE_LABELS[v.variant_type] || v.variant_type)}</td>
                <td><code>${escapeHtml(v.variant)}</code></td>
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
                    <th>Type</th>
                    <th>Change</th>
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

$(document).on('click', '.variant-expand-btn', function() {
    const row = $(this).closest('tr');
    const expanded = row.data('expanded') === true;
    row.data('expanded', !expanded);
    $(this).text(expanded ? '▶' : '▼');
    row.nextUntil('tr:not(.variant-oligo-row)').toggle(!expanded);
});

registerPageModule(
    () => !!document.getElementById('lineage_select'),
    () => { loadConfig(); initBrowser(); },
    () => {
        $('.igv_div').children('.igv-container').remove();
        $('#igv_loading').removeClass('invisible');
        igvBrowser = null;
        tracks = [];
        activeSets = [];
        activeLineageGroup = null;
    }
);
