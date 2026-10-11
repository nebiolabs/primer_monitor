import { currentBrowser } from 'igv_browser';

// The table of variants overlapping primers, shared by the Primer Status and primer set pages. Each page provides
// a #variant_table_section with #variant_table_heading, #variant_table_loading and #variant_table_content.

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
        <table class="table is-fullwidth is-hoverable is-narrow" data-datatable="off">
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

// Fetches variant overlaps from the lineage_variants#variant_overlaps endpoint and renders them. primerSetNames are
// primer_sets.name values (display names), not the URL-safe keys from tracks.json.
export function loadVariantTable(overlapsPath, lineage, primerSetNames, headingText, referenceAccession) {
    const section = document.getElementById('variant_table_section');
    const loading = document.getElementById('variant_table_loading');
    const content = document.getElementById('variant_table_content');
    if (!section) return;

    if (!lineage || primerSetNames.length === 0) {
        section.style.display = 'none';
        return;
    }

    const heading = document.getElementById('variant_table_heading');
    if (heading) heading.textContent = headingText;

    section.style.display = 'block';
    loading.style.display = 'block';
    content.innerHTML = '';
    content.dataset.referenceAccession = referenceAccession || '';

    const params = new URLSearchParams();
    params.set('lineage', lineage);
    primerSetNames.forEach(name => params.append('primer_sets[]', name));

    fetch(`${overlapsPath}?${params}`)
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

// Registered once when this module loads, so they don't stack across Turbo navigations or the pages importing it.
document.addEventListener('click', function(event) {
    const oligo = event.target.closest('.oligo-svg');
    const browser = currentBrowser();
    const reference = document.getElementById('variant_table_content')?.dataset.referenceAccession;
    if (!oligo || !browser || !reference) return;
    const start = Math.max(0, parseInt(oligo.dataset.start) - 50);
    const end = parseInt(oligo.dataset.end) + 50;
    browser.then(igvBrowser => igvBrowser.search(`${reference}:${start}-${end}`))
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
