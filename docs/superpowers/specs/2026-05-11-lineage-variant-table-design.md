# Lineage Variant Table & Primer Visualization

**Date:** 2026-05-11
**Page:** `/organisms/:slug/lineage_variants`

## Overview

Add two linked features to the lineage variants page:

1. A table of variants (filtered to those overlapping primers in the selected primer sets, for the selected lineage) showing position, type, frequency, and global first/last seen dates with location and lineage context.
2. An expandable per-variant primer visualization: each overlapping oligo rendered as an SVG sequence diagram with variant positions highlighted in red, 5′/3′ end labels, and a direction arrow.

Both update when the lineage or primer set selection changes, using the same debounce pattern as the existing IGV browser.

---

## Section 1: Database changes

### Migration 1 — `organisms.pct_cutoff`

Add `pct_cutoff float NOT NULL DEFAULT 1.0` to the organisms table. The pipeline currently receives this as a CLI parameter; moving it to the DB means both the pipeline and the Rails app read from one place. The hard floor of 10 occurrences remains a hardcoded constant (`MINIMUM_OCCURRENCE_COUNT = 10`) in both the shell script and any Ruby code.

### Migration 2 — `lineage_variant_primer_overlaps`

New table populated by the pipeline, read by the Rails endpoint.

```sql
CREATE TABLE lineage_variant_primer_overlaps (
  id                bigserial PRIMARY KEY,
  organism_id       bigint NOT NULL REFERENCES organisms(id),
  lineage_group_key varchar NOT NULL,
  ref_start         integer NOT NULL,
  ref_end           integer NOT NULL,
  variant_type      varchar NOT NULL,
  variant           varchar NOT NULL,
  frequency_pct     float NOT NULL,
  oligo_id          bigint NOT NULL REFERENCES oligos(id),
  created_at        timestamp NOT NULL DEFAULT now()
);

CREATE INDEX ON lineage_variant_primer_overlaps (organism_id, lineage_group_key);

ALTER TABLE lineage_variant_primer_overlaps
  ADD CONSTRAINT lineage_variant_primer_overlaps_unique
  UNIQUE (organism_id, lineage_group_key, ref_start, variant_type, variant, oligo_id);
```

First/last seen (date, location, lineage) is **not** stored here — it changes with every sequence ingestion and is computed at request time via a batch query on existing tables.

---

## Section 2: Pipeline changes

### `pct_cutoff` from DB

`recompute_affected_primers.sh` already reads `variant_bed_lookback_days` from the organisms table via psql. Add the same pattern for `pct_cutoff`. Remove the `--pct_cutoff` CLI parameter from `update_visualization.nf` and `update_visualization_data.sh`.

### New script: `lib/visualization/store_variant_overlaps.sh`

Called from `recompute_affected_primers.sh` once per lineage group, immediately after `count_variants.sh` generates `lineage_variants/{lineage}.bed`.

**Arguments:** organism slug, lineage group key, path to the variants BED file.

**Uses:** the same `$DB_HOST`, `$DB_NAME`, `$DB_USER_RO` env vars already exported by the calling script. Connects with write-capable credentials (`$DB_USER`) for the INSERT.

**Steps:**
1. COPY the variants BED (`chrom, ref_start, ref_end, variant_name, frequency_pct`) into a postgres temp table.
2. Parse `variant_name` (format: `variant_type/variant`, e.g. `X/A`, `D/3-`) into separate columns.
3. Resolve `organism_taxon_id` from the `chrom` column (which is `organism_taxa.reference_accession`).
4. JOIN `oligo_alignment_positions` (filtering by `organism_taxon_id` and overlap condition: `NOT (oap.ref_start >= variant.ref_end OR oap.ref_end <= variant.ref_start)`) → `oligos` → `primer_sets` to find overlapping oligos.
5. DELETE all existing rows for `(organism_id, lineage_group_key)` before inserting, so variants that dropped below threshold are removed.
6. `INSERT INTO lineage_variant_primer_overlaps ... ON CONFLICT DO NOTHING` — the unique constraint acts as a bug/concurrency guard.

**Credentials:** uses `$DB_USER` (write-capable) for the DELETE and INSERT; `$DB_USER_RO` is insufficient. `$DB_USER` must be present in the pipeline `.env` (it is already used by migrations).

---

## Section 3: Rails endpoint

### Route

```ruby
resources :organisms, param: :slug do
  resources :lineage_variants, only: [:index] do
    get :variant_overlaps, on: :collection
  end
end
```

URL: `GET /organisms/:organism_slug/lineage_variants/variant_overlaps.json?lineage=X&primer_sets[]=Y&primer_sets[]=Z`

### Controller action: `LineageVariantsController#variant_overlaps`

1. Look up organism by slug; return 404 if missing.
2. Validate `lineage` and `primer_sets` params present; return 400 otherwise.
3. Apply existing `authorize! :index, LineageVariantsController` guard.
4. **Overlap query:** SELECT from `lineage_variant_primer_overlaps` joined with `oligos`, `oligo_alignment_positions`, and `primer_sets`, filtered to `organism_id`, `lineage_group_key`, and the requested primer set names. Collect overlapping oligos (name, sequence, oligo_start, oligo_end, strand, primer_set_name) grouped by variant.
5. **First/last seen query:** Single batch query for all distinct `(ref_start, variant_type, variant)` in the result. Uses two CTEs with `DISTINCT ON (ref_start, variant_type, variant) ORDER BY ref_start, variant_type, variant, COALESCE(date_collected, date_submitted) ASC, fasta_records.id ASC` (first seen) and `... DESC, fasta_records.id DESC` (last seen) to break ties deterministically. Joined to `variant_sites → fasta_records → lineage_calls → lineages` and `detailed_geo_locations → detailed_geo_location_aliases` to retrieve date, lineage name, and location (region + division).
6. Respond with JSON only; merge first/last seen into the variant objects.

### Response shape

```json
{
  "variants": [
    {
      "ref_start": 22995,
      "ref_end": 22996,
      "variant_type": "X",
      "variant": "A",
      "frequency_pct": 45.2,
      "first_seen": {
        "date": "2021-11-15",
        "lineage": "BA.1",
        "location": "South Africa / Gauteng"
      },
      "last_seen": {
        "date": "2024-01-10",
        "lineage": "JN.1",
        "location": "USA / New York"
      },
      "oligos": [
        {
          "name": "nCoV_400_F",
          "sequence": "ATCGATCG...",
          "oligo_start": 22990,
          "oligo_end": 23010,
          "strand": "+",
          "primer_set": "NEB SARS-CoV-2"
        }
      ]
    }
  ]
}
```

---

## Section 4: Frontend

### Placement

New `<div id="variant_table_section">` in `lineage_variants/index.html.erb`, inserted between the IGV div and the `<hr>` before the scoring explanation. Hidden initially; shown once data loads.

### JS: `lineage_variants.js`

New function `loadVariantTable(lineage, primerSets)` called from `updatePrimerSets()` alongside the existing `loadPrimerSets()` call — same debounce, same trigger conditions. Fetches `variant_overlaps.json`, then re-renders the table section.

On lineage/primer-set change: show a loading indicator in the table section while the fetch is in flight.

### Table

One row per variant:

| Position | Type | Change | Frequency | First Seen | Last Seen |
|---|---|---|---|---|---|
| 22995–22996 | SNP | A | 45.2% | 2021-11-15 · BA.1 · South Africa/Gauteng | 2024-01-10 · JN.1 · USA/New York |

Frontend maps raw `variant_type` values for display: `X` → `SNP`, `D` → `deletion`, `I` → `insertion`.

Each row has an expand toggle. Expanding reveals one sub-row per overlapping oligo containing the oligo name, primer set name, and the SVG sequence visualization.

### SVG sequence visualization

Generated in JS from the oligo JSON data. No server rendering.

**Layout:**
- Fixed-width cells (12 px per base), monospace character per cell.
- Cell fill: dark gray (`#444444`) with white text. Cells where the oligo position falls within `[ref_start, ref_end)` of the variant: red fill (`#cc0000`) with white text.
- `5′` label to the left of the first cell; `3′` label to the right of the last cell — positions determined by `strand` (`+` → 5′ left, 3′ right; `-` → 5′ right, 3′ left).
- Small arrowhead on the 3′ end indicating direction.
- Long sequences scroll horizontally within the container rather than wrapping.

**Overlap calculation:** `overlap = oligo position is within [variant.ref_start, variant.ref_end)`, where oligo position = `oligo_start + character_index` for `+` strand; `oligo_end - 1 - character_index` for `-` strand.
