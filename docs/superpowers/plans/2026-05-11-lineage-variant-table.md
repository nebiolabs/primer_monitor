# Lineage Variant Table & Primer Visualization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a variant table with first/last seen dates and an SVG primer visualization to the lineage variants page, backed by a pipeline-populated DB table.

**Architecture:** The Nextflow pipeline populates `lineage_variant_primer_overlaps` (one row per variant × oligo per lineage group) after each run. A new `variant_overlaps.json` endpoint reads from that table and enriches with first/last seen via a batch query on existing tables. JS fetches and renders the table + SVG visualizations, updating on lineage/primer-set change.

**Spec:** `docs/superpowers/specs/2026-05-11-lineage-variant-table-design.md`

**Tech Stack:** Rails 7 / Minitest / PostgreSQL / vanilla JS + jQuery / inline SVG

---

## File Map

| Action | Path |
|---|---|
| Create | `db/migrate/20260511100000_add_pct_cutoff_to_organisms.rb` |
| Create | `db/migrate/20260511100001_create_lineage_variant_primer_overlaps.rb` |
| Create | `app/models/lineage_variant_primer_overlap.rb` |
| Create | `test/models/lineage_variant_primer_overlap_test.rb` |
| Create | `test/fixtures/lineage_variant_primer_overlaps.yml` |
| Create | `test/controllers/lineage_variants_controller_test.rb` |
| Create | `lib/visualization/store_variant_overlaps.sh` |
| Modify | `test/fixtures/organisms.yml` — add `pct_cutoff` |
| Modify | `config/routes.rb` — add `variant_overlaps` collection route |
| Modify | `app/controllers/lineage_variants_controller.rb` — add `variant_overlaps` action |
| Modify | `app/views/lineage_variants/index.html.erb` — add `#variant_table_section` |
| Modify | `app/javascript/lineage_variants.js` — add fetch + table + SVG rendering |
| Modify | `lib/visualization/recompute_affected_primers.sh` — read `pct_cutoff` from DB; call `store_variant_overlaps.sh` |
| Modify | `lib/visualization/update_visualization_data.sh` — remove `pct_cutoff` CLI param |
| Modify | `lib/update_visualization.nf` — remove `pct_cutoff` param |

---

## Task 1: Add pct_cutoff to organisms

**Files:**
- Create: `db/migrate/20260511100000_add_pct_cutoff_to_organisms.rb`
- Modify: `test/fixtures/organisms.yml`

- [ ] **Step 1: Write the failing test**

Add to `test/models/lineage_variant_primer_overlap_test.rb` (create the file):

```ruby
# frozen_string_literal: true

require 'test_helper'

class LineageVariantPrimerOverlapTest < ActiveSupport::TestCase
  test 'organism has pct_cutoff' do
    assert_respond_to organisms(:sars_cov2), :pct_cutoff
    assert_equal 1.0, organisms(:sars_cov2).pct_cutoff
  end
end
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
bin/rails test test/models/lineage_variant_primer_overlap_test.rb
```

Expected: FAIL — `undefined method 'pct_cutoff'`

- [ ] **Step 3: Create migration**

```ruby
# db/migrate/20260511100000_add_pct_cutoff_to_organisms.rb
# frozen_string_literal: true

class AddPctCutoffToOrganisms < ActiveRecord::Migration[7.0]
  def change
    add_column :organisms, :pct_cutoff, :float, null: false, default: 1.0
  end
end
```

- [ ] **Step 4: Run migration**

```bash
bin/rails db:migrate
```

- [ ] **Step 5: Update organisms fixture**

In `test/fixtures/organisms.yml`, add `pct_cutoff: 1.0` to the `sars_cov2` entry:

```yaml
sars_cov2:
  name: "SARS-CoV-2"
  slug: "sars-cov-2"
  public: true
  pct_cutoff: 1.0
```

- [ ] **Step 6: Run test to confirm it passes**

```bash
bin/rails test test/models/lineage_variant_primer_overlap_test.rb
```

Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add db/migrate/20260511100000_add_pct_cutoff_to_organisms.rb \
        db/structure.sql \
        test/fixtures/organisms.yml \
        test/models/lineage_variant_primer_overlap_test.rb
git commit -m "Add pct_cutoff column to organisms table"
```

---

## Task 2: Create lineage_variant_primer_overlaps table and model

**Files:**
- Create: `db/migrate/20260511100001_create_lineage_variant_primer_overlaps.rb`
- Create: `app/models/lineage_variant_primer_overlap.rb`
- Create: `test/fixtures/lineage_variant_primer_overlaps.yml`

- [ ] **Step 1: Add model validity test to existing test file**

Append to `test/models/lineage_variant_primer_overlap_test.rb`:

```ruby
  test 'valid overlap can be created' do
    overlap = lineage_variant_primer_overlaps(:one)
    assert overlap.valid?
    assert_equal 'XBB', overlap.lineage_group_key
    assert_equal 100, overlap.ref_start
    assert_equal 'X', overlap.variant_type
  end
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
bin/rails test test/models/lineage_variant_primer_overlap_test.rb
```

Expected: FAIL — table or fixture not found

- [ ] **Step 3: Create migration**

```ruby
# db/migrate/20260511100001_create_lineage_variant_primer_overlaps.rb
# frozen_string_literal: true

class CreateLineageVariantPrimerOverlaps < ActiveRecord::Migration[7.0]
  def change
    create_table :lineage_variant_primer_overlaps do |t|
      t.references :organism, null: false, foreign_key: true
      t.string :lineage_group_key, null: false
      t.integer :ref_start, null: false
      t.integer :ref_end, null: false
      t.string :variant_type, null: false
      t.string :variant, null: false
      t.float :frequency_pct, null: false
      t.references :oligo, null: false, foreign_key: true
      t.datetime :created_at, null: false, default: -> { 'NOW()' }
    end

    add_index :lineage_variant_primer_overlaps,
              %i[organism_id lineage_group_key]

    add_index :lineage_variant_primer_overlaps,
              %i[organism_id lineage_group_key ref_start variant_type variant oligo_id],
              unique: true,
              name: 'lineage_variant_primer_overlaps_unique'
  end
end
```

- [ ] **Step 4: Run migration**

```bash
bin/rails db:migrate
```

- [ ] **Step 5: Create model**

```ruby
# app/models/lineage_variant_primer_overlap.rb
# frozen_string_literal: true

class LineageVariantPrimerOverlap < ApplicationRecord
  belongs_to :organism
  belongs_to :oligo
end
```

- [ ] **Step 6: Create fixture**

```yaml
# test/fixtures/lineage_variant_primer_overlaps.yml
one:
  organism: sars_cov2
  lineage_group_key: "XBB"
  ref_start: 100
  ref_end: 101
  variant_type: "X"
  variant: "T"
  frequency_pct: 15.5
  oligo: one
```

- [ ] **Step 7: Run test to confirm it passes**

```bash
bin/rails test test/models/lineage_variant_primer_overlap_test.rb
```

Expected: both tests PASS

- [ ] **Step 8: Commit**

```bash
git add db/migrate/20260511100001_create_lineage_variant_primer_overlaps.rb \
        db/structure.sql \
        app/models/lineage_variant_primer_overlap.rb \
        test/models/lineage_variant_primer_overlap_test.rb \
        test/fixtures/lineage_variant_primer_overlaps.yml
git commit -m "Add lineage_variant_primer_overlaps table and model"
```

---

## Task 3: Add route and controller action skeleton

**Files:**
- Modify: `config/routes.rb`
- Create: `test/controllers/lineage_variants_controller_test.rb`
- Modify: `app/controllers/lineage_variants_controller.rb`

- [ ] **Step 1: Write failing controller test**

```ruby
# test/controllers/lineage_variants_controller_test.rb
# frozen_string_literal: true

require 'test_helper'

class LineageVariantsControllerTest < ActionDispatch::IntegrationTest
  setup do
    @organism = organisms(:sars_cov2)
    sign_in(users(:admin_user))
  end

  test 'variant_overlaps returns 400 without lineage param' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :bad_request
  end

  test 'variant_overlaps returns 400 without primer_sets param' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      format: :json
    )
    assert_response :bad_request
  end

  test 'variant_overlaps returns 404 for unknown organism' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: 'no-such-organism',
      lineage: 'XBB',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :not_found
  end

  test 'variant_overlaps returns JSON for valid params' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :success
    data = JSON.parse(response.body)
    assert data.key?('variants')
  end
end
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
bin/rails test test/controllers/lineage_variants_controller_test.rb
```

Expected: FAIL — route not found

- [ ] **Step 3: Add route**

In `config/routes.rb`, change:

```ruby
resources :lineage_variants, only: [:index]
```

to:

```ruby
resources :lineage_variants, only: [:index] do
  get :variant_overlaps, on: :collection
end
```

- [ ] **Step 4: Add skeleton action to controller**

Append to the public section of `app/controllers/lineage_variants_controller.rb` (before the `private` keyword):

```ruby
  def variant_overlaps
    @organism = Organism.find_by(slug: params[:organism_slug])
    return head :not_found unless @organism

    authorize! :index, LineageVariantsController

    lineage_param = params[:lineage].presence
    primer_set_params = Array(params[:primer_sets]).reject(&:blank?)

    return head :bad_request if lineage_param.nil? || primer_set_params.empty?

    render json: { variants: [] }
  end
```

- [ ] **Step 5: Run tests to confirm they pass**

```bash
bin/rails test test/controllers/lineage_variants_controller_test.rb
```

Expected: all 4 tests PASS

- [ ] **Step 6: Commit**

```bash
git add config/routes.rb \
        app/controllers/lineage_variants_controller.rb \
        test/controllers/lineage_variants_controller_test.rb
git commit -m "Add variant_overlaps route and controller skeleton"
```

---

## Task 4: Implement overlap query

**Files:**
- Modify: `app/controllers/lineage_variants_controller.rb`
- Modify: `test/controllers/lineage_variants_controller_test.rb`

The fixture `lineage_variant_primer_overlaps(:one)` links organism `sars_cov2`, lineage_group_key `"XBB"`, oligo `one` (which belongs to primer_set `Charité`).

- [ ] **Step 1: Add failing test for overlap data**

Append to `test/controllers/lineage_variants_controller_test.rb`:

```ruby
  test 'variant_overlaps returns matching variants with oligos' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :success
    data = JSON.parse(response.body)
    assert_equal 1, data['variants'].length

    v = data['variants'].first
    assert_equal 100, v['ref_start']
    assert_equal 101, v['ref_end']
    assert_equal 'X', v['variant_type']
    assert_equal 'T', v['variant']
    assert_in_delta 15.5, v['frequency_pct'], 0.01

    assert_equal 1, v['oligos'].length
    o = v['oligos'].first
    assert_equal 'probe1', o['name']
    assert_equal 'GCAATTTATATACATATA', o['sequence']
    assert_equal 'Charité', o['primer_set']
  end

  test 'variant_overlaps excludes variants for other primer sets' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      primer_sets: ['CDC'],
      format: :json
    )
    assert_response :success
    data = JSON.parse(response.body)
    assert_equal 0, data['variants'].length
  end

  test 'variant_overlaps excludes variants for other lineage groups' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'JN.1',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :success
    data = JSON.parse(response.body)
    assert_equal 0, data['variants'].length
  end
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
bin/rails test test/controllers/lineage_variants_controller_test.rb
```

Expected: the 3 new tests FAIL — variants array is always empty

- [ ] **Step 3: Implement the overlap query**

Replace the `variant_overlaps` action in `app/controllers/lineage_variants_controller.rb`:

```ruby
  def variant_overlaps
    @organism = Organism.find_by(slug: params[:organism_slug])
    return head :not_found unless @organism

    authorize! :index, LineageVariantsController

    lineage_param = params[:lineage].presence
    primer_set_params = Array(params[:primer_sets]).reject(&:blank?)

    return head :bad_request if lineage_param.nil? || primer_set_params.empty?

    rows = LineageVariantPrimerOverlap
      .joins(oligo: [:primer_set, { oligo_alignment_positions: :organism_taxon }])
      .where(
        lineage_variant_primer_overlaps: {
          organism_id: @organism.id,
          lineage_group_key: lineage_param
        },
        primer_sets: { name: primer_set_params },
        organism_taxa: { organism_id: @organism.id }
      )
      .select(
        'lineage_variant_primer_overlaps.ref_start',
        'lineage_variant_primer_overlaps.ref_end',
        'lineage_variant_primer_overlaps.variant_type',
        'lineage_variant_primer_overlaps.variant',
        'lineage_variant_primer_overlaps.frequency_pct',
        'oligos.name AS oligo_name',
        'oligos.sequence AS oligo_sequence',
        'oligos.strand AS oligo_strand',
        'oligo_alignment_positions.ref_start AS oligo_start',
        'oligo_alignment_positions.ref_end AS oligo_end',
        'primer_sets.name AS primer_set_name'
      )

    variants_map = {}
    rows.each do |row|
      key = [row.ref_start, row.ref_end, row.variant_type, row.variant]
      variants_map[key] ||= {
        ref_start: row.ref_start,
        ref_end: row.ref_end,
        variant_type: row.variant_type,
        variant: row.variant,
        frequency_pct: row.frequency_pct,
        oligos: []
      }
      variants_map[key][:oligos] << {
        name: row.oligo_name,
        sequence: row.oligo_sequence,
        oligo_start: row.oligo_start,
        oligo_end: row.oligo_end,
        strand: row.oligo_strand,
        primer_set: row.primer_set_name
      }
    end

    render json: { variants: variants_map.values }
  end
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
bin/rails test test/controllers/lineage_variants_controller_test.rb
```

Expected: all tests PASS

- [ ] **Step 5: Commit**

```bash
git add app/controllers/lineage_variants_controller.rb \
        test/controllers/lineage_variants_controller_test.rb
git commit -m "Implement variant_overlaps overlap query"
```

---

## Task 5: Add first/last seen to response

**Files:**
- Modify: `app/controllers/lineage_variants_controller.rb`
- Modify: `test/controllers/lineage_variants_controller_test.rb`

The existing `fasta_records` fixtures have no `date_collected`/`date_submitted`, so the first/last seen values will be `nil` in tests — that's fine; we're verifying the keys exist and the query doesn't crash.

- [ ] **Step 1: Add test asserting first_seen/last_seen keys are present**

Append to `test/controllers/lineage_variants_controller_test.rb`:

```ruby
  test 'variant_overlaps includes first_seen and last_seen keys' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :success
    v = JSON.parse(response.body)['variants'].first
    assert v.key?('first_seen'), 'missing first_seen'
    assert v.key?('last_seen'),  'missing last_seen'
    assert v['first_seen'].key?('date')
    assert v['first_seen'].key?('lineage')
    assert v['first_seen'].key?('location')
  end
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
bin/rails test test/controllers/lineage_variants_controller_test.rb
```

Expected: FAIL — `first_seen` key missing

- [ ] **Step 3: Implement first/last seen query**

Add a private method to `app/controllers/lineage_variants_controller.rb`:

```ruby
  def fetch_first_last_seen(variants)
    return {} if variants.empty?

    conditions = variants.map do |v|
      ApplicationRecord.sanitize_sql_array(
        ['(vs.ref_start = ? AND vs.variant_type = ? AND vs.variant = ?)',
         v[:ref_start], v[:variant_type], v[:variant]]
      )
    end.join(' OR ')

    sql = <<~SQL
      WITH observations AS (
        SELECT
          vs.ref_start,
          vs.variant_type,
          vs.variant,
          COALESCE(fr.date_collected, fr.date_submitted) AS obs_date,
          l.name AS lineage_name,
          COALESCE(dga.region, '') ||
            CASE WHEN dga.division IS NOT NULL THEN ' / ' || dga.division ELSE '' END AS location,
          ROW_NUMBER() OVER (
            PARTITION BY vs.ref_start, vs.variant_type, vs.variant
            ORDER BY COALESCE(fr.date_collected, fr.date_submitted) ASC NULLS LAST, fr.id ASC
          ) AS rn_first,
          ROW_NUMBER() OVER (
            PARTITION BY vs.ref_start, vs.variant_type, vs.variant
            ORDER BY COALESCE(fr.date_collected, fr.date_submitted) DESC NULLS LAST, fr.id DESC
          ) AS rn_last
        FROM variant_sites vs
        JOIN fasta_records fr ON fr.id = vs.fasta_record_id
        JOIN lineage_calls lc ON lc.id = fr.lineage_call_id
        JOIN lineages l ON l.id = lc.lineage_id
        JOIN detailed_geo_locations dgl ON dgl.id = fr.detailed_geo_location_id
        JOIN detailed_geo_location_aliases dga ON dga.id = dgl.detailed_geo_location_alias_id
        WHERE #{conditions}
      )
      SELECT
        ref_start, variant_type, variant,
        MAX(CASE WHEN rn_first = 1 THEN obs_date END)    AS first_date,
        MAX(CASE WHEN rn_first = 1 THEN lineage_name END) AS first_lineage,
        MAX(CASE WHEN rn_first = 1 THEN location END)    AS first_location,
        MAX(CASE WHEN rn_last  = 1 THEN obs_date END)    AS last_date,
        MAX(CASE WHEN rn_last  = 1 THEN lineage_name END) AS last_lineage,
        MAX(CASE WHEN rn_last  = 1 THEN location END)    AS last_location
      FROM observations
      WHERE rn_first = 1 OR rn_last = 1
      GROUP BY ref_start, variant_type, variant
    SQL

    ApplicationRecord.connection.select_all(sql).each_with_object({}) do |row, h|
      key = [row['ref_start'].to_i, row['variant_type'], row['variant']]
      h[key] = {
        first_seen: { date: row['first_date'], lineage: row['first_lineage'], location: row['first_location'] },
        last_seen:  { date: row['last_date'],  lineage: row['last_lineage'],  location: row['last_location'] }
      }
    end
  end
```

- [ ] **Step 4: Call it in variant_overlaps and merge results**

In the `variant_overlaps` action, replace `render json: { variants: variants_map.values }` with:

```ruby
    variants = variants_map.values
    first_last = fetch_first_last_seen(variants)

    variants.each do |v|
      seen = first_last.fetch([v[:ref_start], v[:variant_type], v[:variant]],
                              { first_seen: { date: nil, lineage: nil, location: nil },
                                last_seen:  { date: nil, lineage: nil, location: nil } })
      v.merge!(seen)
    end

    render json: { variants: }
```

- [ ] **Step 5: Run tests to confirm all pass**

```bash
bin/rails test test/controllers/lineage_variants_controller_test.rb
```

Expected: all tests PASS

- [ ] **Step 6: Commit**

```bash
git add app/controllers/lineage_variants_controller.rb \
        test/controllers/lineage_variants_controller_test.rb
git commit -m "Add first/last seen to variant_overlaps response"
```

---

## Task 6: ERB view update

**Files:**
- Modify: `app/views/lineage_variants/index.html.erb`

- [ ] **Step 1: Add variant table section to the view**

In `app/views/lineage_variants/index.html.erb`, find the line:

```erb
<hr/>
```

(the `<hr/>` before the "Calculating variant impact" block) and insert before it:

```erb
<div id="variant_table_section" class="block" style="display:none">
  <h2 class="title is-5">Variants Overlapping Selected Primers</h2>
  <div id="variant_table_loading" class="has-text-centered" style="display:none; padding: 1rem;">
    <progress class="progress is-primary" style="max-width: 16rem; margin: 0 auto;" max="100"></progress>
  </div>
  <div id="variant_table_content"></div>
</div>
```

- [ ] **Step 2: Verify the page still loads**

Start the server and visit `/organisms/sars-cov-2/lineage_variants`. The new `div` is hidden — page should look identical to before. Check browser console for JS errors.

```bash
bin/dev
```

- [ ] **Step 3: Commit**

```bash
git add app/views/lineage_variants/index.html.erb
git commit -m "Add variant_table_section placeholder to lineage_variants view"
```

---

## Task 7: JS fetch and table rendering

**Files:**
- Modify: `app/javascript/lineage_variants.js`

- [ ] **Step 1: Add helper functions, stub renderOligoSvg, and loadVariantTable**

Append to `app/javascript/lineage_variants.js` before the `registerPageModule` call:

`renderOligoSvg` is stubbed here so `buildVariantTable` doesn't throw; Task 8 replaces it with the real implementation.

```javascript
const VARIANT_TYPE_LABELS = { X: 'SNP', D: 'deletion', I: 'insertion' };

// Stub — replaced by real implementation in Task 8
function renderOligoSvg(_sequence, _oligoStart, _oligoEnd, _strand, _variantStart, _variantEnd) {
    return '';
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

    section.style.display = 'block';
    loading.style.display = 'block';
    content.innerHTML = '';

    const params = new URLSearchParams();
    params.set('lineage', lineage);
    primerSets.forEach(ps => params.append('primer_sets[]', ps));

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
```

- [ ] **Step 2: Wire loadVariantTable into updatePrimerSets**

Find the `updatePrimerSets` function and add a call to `loadVariantTable` after the existing `loadPrimerSets` call:

```javascript
function updatePrimerSets() {
    // ... existing code ...
    activeLineageGroup = $('#lineage_select').val();

    if (igvBrowser != null) {
        activeSets = $('#primer_set_select').val() || [];
        loadPrimerSets(activeSets, igvBrowser, activeLineageGroup);
        loadVariantTable(activeLineageGroup, activeSets);  // ← add this line
    }
}
```

- [ ] **Step 3: Wire expand toggles with event delegation**

Add to the `$(document).on(...)` section at the bottom of the file (before `registerPageModule`):

```javascript
$(document).on('click', '.variant-expand-btn', function() {
    const row = $(this).closest('tr');
    const expanded = row.data('expanded') === true;
    row.data('expanded', !expanded);
    $(this).text(expanded ? '▶' : '▼');
    row.nextUntil('tr:not(.variant-oligo-row)').toggle(!expanded);
});
```

- [ ] **Step 4: Call loadVariantTable on initial page load**

In the `initBrowser` function's `.then` callback, after `loadPrimerSets(activeSets, igvBrowser, activeLineageGroup)`, add:

```javascript
        if (activeLineageGroup) {
            loadPrimerSets(activeSets, igvBrowser, activeLineageGroup);
            loadVariantTable(activeLineageGroup, activeSets);  // ← add this line
        }
```

- [ ] **Step 5: Verify in browser**

Start the server, visit `/organisms/sars-cov-2/lineage_variants`, select a lineage and primer set. The variant table should appear below the IGV browser. Open browser DevTools Network tab and confirm a request to `variant_overlaps.json` is made when selections change.

```bash
bin/dev
```

- [ ] **Step 6: Commit**

```bash
git add app/javascript/lineage_variants.js
git commit -m "Add variant table fetch and rendering to lineage_variants page"
```

---

## Task 8: SVG sequence visualization

**Files:**
- Modify: `app/javascript/lineage_variants.js`

- [ ] **Step 1: Replace the renderOligoSvg stub with the real implementation**

Find the stub added in Task 7:

```javascript
// Stub — replaced by real implementation in Task 8
function renderOligoSvg(_sequence, _oligoStart, _oligoEnd, _strand, _variantStart, _variantEnd) {
    return '';
}
```

Replace it with the full implementation:

```javascript
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

    // 5' is left for + strand, right for - strand
    const fivePrimeX  = isPlus ? 0              : LABEL_W + seqW + ARROW_W;
    const threePrimeX = isPlus ? LABEL_W + seqW + ARROW_W : 0;

    // Arrowhead polygon points (direction of 3' end)
    const arrowX = isPlus ? LABEL_W + seqW : LABEL_W;
    const mid = CELL_H / 2;
    const arrowPoints = isPlus
        ? `${arrowX},${mid - 4} ${arrowX + ARROW_W},${mid} ${arrowX},${mid + 4}`
        : `${arrowX + ARROW_W},${mid - 4} ${arrowX},${mid} ${arrowX + ARROW_W},${mid + 4}`;

    const cells = sequence.split('').map((base, i) => {
        // Genomic position of this character
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
```

- [ ] **Step 2: Verify SVG in browser**

With the server running, select a lineage and primer set that has overlapping variants. Expand a variant row. You should see the SVG sequence visualization: dark gray cells for non-variant positions, red cells for variant positions, `5'` and `3'` labels, and an arrow at the 3′ end.

Test both `+` and `-` strand oligos if available — confirm the `5'` label is on the correct side.

Quick console test (paste into browser DevTools):

```javascript
document.getElementById('variant_table_content').innerHTML +=
  renderOligoSvg('ATCGATCGATCG', 100, 112, '+', 104, 106);
```

Expected: 12 gray cells with cells at index 4 and 5 (positions 104–105) in red, `5'` on left, `3'` + arrow on right.

- [ ] **Step 3: Commit**

```bash
git add app/javascript/lineage_variants.js
git commit -m "Add SVG oligo sequence visualization with variant highlighting"
```

---

## Task 9: Pipeline script — store_variant_overlaps.sh

**Files:**
- Create: `lib/visualization/store_variant_overlaps.sh`

- [ ] **Step 1: Create the script**

```bash
#!/usr/bin/env bash

# Stores variant-primer overlaps for one lineage group into the DB.
#
# Usage: store_variant_overlaps.sh <organism_slug> <lineage_group_key> <variants_bed_path>
#
# Requires: DB_HOST, DB_NAME, DB_USER (write-capable) env vars exported by caller.
# variants_bed_path: the lineage_variants/{lineage}.bed file produced by count_variants.sh
#   Columns: chrom, ref_start, ref_end, variant_name (type/allele), frequency_pct

set -e

organism_slug="$1"
lineage_group_key="$2"
variants_bed="$3"

if [[ -z "$organism_slug" || -z "$lineage_group_key" || -z "$variants_bed" ]]; then
  echo "Usage: store_variant_overlaps.sh <organism_slug> <lineage_group_key> <variants_bed_path>" >&2
  exit 1
fi

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" \
  -v "organism_slug=$organism_slug" \
  -v "lineage_key=$lineage_group_key" \
  -v "variants_bed=$variants_bed" \
  -f - <<'PSQL'

CREATE TEMP TABLE tmp_variants (
  chrom        varchar,
  ref_start    integer,
  ref_end      integer,
  variant_name varchar,
  frequency_pct float
) ON COMMIT DROP;

\COPY tmp_variants FROM :'variants_bed' WITH (FORMAT csv, DELIMITER E'\t', HEADER false)

-- Remove stale rows from previous pipeline runs for this organism + lineage group
DELETE FROM lineage_variant_primer_overlaps
WHERE organism_id = (SELECT id FROM organisms WHERE slug = :'organism_slug')
  AND lineage_group_key = :'lineage_key';

-- Insert overlaps: variant positions x oligos aligned to the same reference
INSERT INTO lineage_variant_primer_overlaps
  (organism_id, lineage_group_key, ref_start, ref_end,
   variant_type, variant, frequency_pct, oligo_id)
SELECT
  org.id,
  :'lineage_key',
  tv.ref_start,
  tv.ref_end,
  split_part(tv.variant_name, '/', 1),
  split_part(tv.variant_name, '/', 2),
  tv.frequency_pct,
  o.id
FROM tmp_variants tv
JOIN organism_taxa ot  ON ot.reference_accession = tv.chrom
JOIN organisms org     ON org.id = ot.organism_id AND org.slug = :'organism_slug'
JOIN oligo_alignment_positions oap
  ON  oap.organism_taxon_id = ot.id
  AND NOT (oap.ref_start >= tv.ref_end OR oap.ref_end <= tv.ref_start)
JOIN oligos o          ON o.id = oap.oligo_id
JOIN primer_sets ps    ON ps.id = o.primer_set_id AND ps.organism_id = org.id
ON CONFLICT DO NOTHING;

PSQL
```

- [ ] **Step 2: Make it executable**

```bash
chmod +x lib/visualization/store_variant_overlaps.sh
```

- [ ] **Step 3: Smoke-test manually (if DB is accessible)**

```bash
# Requires: DB_HOST, DB_NAME, DB_USER set, and a real variants BED file available
# Example:
export DB_HOST=localhost DB_NAME=primer_monitor_dev DB_USER=primer_monitor
lib/visualization/store_variant_overlaps.sh sars-cov-2 XBB /path/to/XBB.bed
```

Expected: psql executes without error; rows appear in `lineage_variant_primer_overlaps`.

- [ ] **Step 4: Commit**

```bash
git add lib/visualization/store_variant_overlaps.sh
git commit -m "Add store_variant_overlaps.sh pipeline script"
```

---

## Task 10: Update pipeline to use DB pct_cutoff and call store script

**Files:**
- Modify: `lib/visualization/recompute_affected_primers.sh`
- Modify: `lib/visualization/update_visualization_data.sh`
- Modify: `lib/update_visualization.nf`

- [ ] **Step 1: Update recompute_affected_primers.sh**

Read `pct_cutoff` from the DB instead of accepting it as a CLI argument, and call `store_variant_overlaps.sh` after each `count_variants.sh` run.

In `lib/visualization/recompute_affected_primers.sh`:

Change the argument list comment and variable assignment at the top:

```bash
# Usage: recompute_affected_primers.sh <cutoff date> <output path> <score cutoff> <primer sets list> <organism slug> <cpus> [lineage sets...]
# (pct_cutoff is now read from the organisms DB table, not passed as a CLI arg)

cutoff_date="$1";
output_path="$2";
score_cutoff="$3";   # was previously arg 4; shift everything left by 1
primer_sets_list_path="$4";
organism_slug="$5";
threads="$6"
```

After the `lookback` block (which already reads `variant_bed_lookback_days`), add:

```bash
min_pct=$(psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -v "organism_slug=$organism_slug" \
  <<< "SELECT pct_cutoff FROM organisms WHERE slug=:'organism_slug';" --csv -t)
```

Update the `shift 7` to `shift 6` (one fewer positional arg).

In the per-lineage-set loop, after the `count_variants.sh` call and before `process_primer_sets.sh`, add:

```bash
  lineage_variants_bed="$output_path/lineage_variants/$lineage_set_name.bed"
  "$(dirname "$0")/store_variant_overlaps.sh" "$organism_slug" "$lineage_set_name" "$lineage_variants_bed"
```

- [ ] **Step 2: Update update_visualization_data.sh**

Remove `pct_cutoff` as a positional argument and stop passing it to `recompute_affected_primers.sh`.

In `lib/visualization/update_visualization_data.sh`:

Remove the line:
```bash
pct_cutoff="$3"
```

Change the script parameter comment from:
```
<primer_monitor_path> <organism_slug> <pct_cutoff> <score_cutoff> <cpus>
```
to:
```
<primer_monitor_path> <organism_slug> <score_cutoff> <cpus>
```

Update the positional assignments:
```bash
primer_monitor_path="$1"
organism_dirname="$2"
score_cutoff="$3"   # was $4
cpus="$4"           # was $5
```

Update the `recompute_affected_primers.sh` call (remove the `"$pct_cutoff"` argument):
```bash
| xargs "$primer_monitor_path/lib/visualization/recompute_affected_primers.sh" - "./$organism_dirname" "$score_cutoff" \
"$primer_sets_data_file" "$organism_dirname" "$cpus"
```

- [ ] **Step 3: Update update_visualization.nf**

In `lib/update_visualization.nf`:

Remove:
```
assert params.pct_cutoff != null : "--pct_cutoff must be specified"
pct_cutoff = params.pct_cutoff
```

Update the shell command in `update_visualization_data` process (remove `!{pct_cutoff}` from the call):
```groovy
    !{primer_monitor_path}/lib/visualization/update_visualization_data.sh -o !{override_path} !{primer_monitor_path} \
    !{organism} !{score_cutoff} !{task.cpus}
```

- [ ] **Step 4: Verify pipeline changes are consistent**

Check there are no remaining references to `pct_cutoff` as a CLI parameter:

```bash
grep -n "pct_cutoff" \
  lib/visualization/recompute_affected_primers.sh \
  lib/visualization/update_visualization_data.sh \
  lib/update_visualization.nf
```

Expected: references only in the `psql` DB query line in `recompute_affected_primers.sh` and in the `pct_cutoff` Nextflow parameter check that was removed. Any remaining hits should be the DB-query usage, not CLI param passing.

- [ ] **Step 5: Run the full test suite**

```bash
bin/rails test
```

Expected: all tests PASS

- [ ] **Step 6: Commit**

```bash
git add lib/visualization/recompute_affected_primers.sh \
        lib/visualization/update_visualization_data.sh \
        lib/update_visualization.nf
git commit -m "Read pct_cutoff from DB and populate lineage_variant_primer_overlaps in pipeline"
```
