# Amplicon Dropout by Lineage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** For tiled amplicon primer sets, report per lineage the fraction of sequences in which each amplicon still
amplifies — i.e. at least one forward and at least one reverse primer of that amplicon bind — instead of flagging
every primer that overlaps a variant.

**Why:** Today a primer is "affected" if any variant falls inside its reference footprint
(`variant_overlaps`, `lineage_variant_primer_overlaps`). The primer's own sequence is never consulted. ALT primers
are usually designed to carry the variant base, so they sit exactly on the variant they fix and get flagged too: the
remedy shows up as a second problem, and a set patched with ALTs looks worse than the unpatched set.

**Architecture:** Oligos get an amplicon and a side (reusing `oligos.locus` and `oligos.category`
Forward/Reverse), parsed from names with a manual override. Primer alignment records the strand. A new step in
the visualization pipeline walks the per-sequence variants already extracted for each lineage group, reconstructs
each sequence's bases under every primer, compares them to the primer's own sequence, and counts per amplicon how
many sequences have ≥1 working forward and ≥1 working reverse primer. Results go to a new
`lineage_amplicon_coverage` table read by the primer set page.

**Tech Stack:** Rails 7 / Minitest / PostgreSQL / Nextflow / bash + Python (pipeline) / vanilla JS

---

## Current state (verified 2026-10-11 against staging)

- Tiled amplicon sets (NEBNext 105, VarSkip 2d/2e 102/103, Swift 99, RSV 101/104) have `locus`, `category`, `strand`
  all NULL. Amplicon and side exist only in the name, in several conventions:
  - `varskip-0317-1_01_LEFT`, `…_01_LEFT_ALT3`, `…_01_LEFT_Alt3`
  - `nCoV-2019_7_LEFT_alt0`, `nCoV-2019_7_RIGHT_alt5` (ARTICv3, set 16 — `category`/`strand` already filled)
  - `PM1_rsva_Amp1a_l`, `PM1_rsva_Amp1a_l_alt1`, `PM1_rsvb_Amp3b-alt1_l`, `PM1_rsva_Amp7a_alt_r`
  - RSV sets mix rsva and rsvb amplicons (`Amp3a` vs `Amp3b`) in one set; these are different taxa.
- `oligo_alignment_positions` stores `ref_start`/`ref_end` per taxon but not the strand the primer aligned on
  (bowtie2 in `lib/update_primers.sh`), so a reverse primer can't be compared to the reference without guessing.
- `lib/visualization/extract_all_variants.sh` already emits one row per (sequence, variant) with lineage;
  `count_variants.sh` filters it per lineage group via `filter_variants.awk`. Sequences with no variants at all
  don't appear in that output, so per-group sequence totals must come from the lineage membership, not the
  variant rows.
- `initial_score` treats a variant within 5 bases of the 3′ end as critical; reuse that as the default.

## Definitions (confirmed 2026-10-11)

- **Primer works on a sequence** if, comparing the sequence's bases at the primer site to the primer sequence
  (reverse-complemented for minus-strand primers):
  - no mismatch, insertion or deletion in the 3′-most `three_prime_window` bases (default 5), and
  - total mismatches ≤ `max_mismatches` (default 2); any indel inside the primer counts as failing.
  - IUPAC degenerate bases in the primer match any compatible base.
  - `N` / ambiguous calls in the sequence: treat as unknown, not mismatch; if any of the 3′ window is unknown the
    sequence is excluded from that amplicon's denominator.
- **Amplicon works on a sequence** if ≥1 Forward oligo and ≥1 Reverse oligo of that amplicon work.
- **Dropout flag:** amplicon works in < `dropout_threshold` (default 90%) of a lineage group's sequences.

Store the three parameters on `organisms` next to `pct_cutoff` so they're tunable without a deploy.

---

## File Map

| Action | Path |
|---|---|
| Create | `app/models/concerns/oligo_name_parser.rb` (or `app/services/oligo_name_parser.rb`) |
| Create | `test/services/oligo_name_parser_test.rb` |
| Create | `lib/tasks/backfill_oligo_amplicons.rake` |
| Create | `db/migrate/…_add_strand_to_oligo_alignment_positions.rb` |
| Create | `db/migrate/…_add_amplicon_params_to_organisms.rb` |
| Create | `db/migrate/…_create_lineage_amplicon_coverage.rb` |
| Create | `app/models/lineage_amplicon_coverage.rb` + test + fixtures |
| Create | `lib/visualization/amplicon_coverage.py` + `test/lib/amplicon_coverage_test.py` (pytest) |
| Create | `lib/visualization/store_amplicon_coverage.sh` |
| Modify | `app/models/oligo.rb` — fill locus/category from parser on save when blank |
| Modify | `app/views/primer_sets/_form` / oligo fields — show and allow editing amplicon + side |
| Modify | `lib/update_primers.sh` — capture SAM strand into `oligo_alignment_positions.strand` |
| Modify | `lib/visualization/update_visualization_data.sh` — run coverage step per lineage group |
| Modify | `lib/visualization/extract_all_variants.sh` — also emit lineage membership (all sequences) |
| Modify | `app/controllers/primer_sets_controller.rb` — `amplicon_coverage.json` |
| Modify | `app/views/primer_sets/show.html.erb`, `app/javascript/primer_sets_show.js` — amplicon table |
| Modify | `db/structure.sql` (regenerated by migrations) |

---

## Task 1: Amplicon and side for every oligo

- [ ] Write `OligoNameParser.parse(name) → { amplicon:, side:, alt: }` with tests covering every convention
      listed above (pull the full distinct name list from staging into the test as fixtures). Side tokens:
      `LEFT|_l$|F` → Forward, `RIGHT|_r$|R` → Reverse. Strip ALT markers (`_ALT\d*`, `_Alt\d*`, `_alt\d*`,
      `-alt\d*`, `_alt_`) before extracting the amplicon. Keep the taxon token (`rsva`/`rsvb`/`panr`) as part of
      the amplicon key so RSV-A and RSV-B amplicons don't merge. Return nil fields when unsure rather than guess.
- [ ] `Oligo` `before_validation`: when `locus` / `category` blank and amplification method is Tiled amplicon,
      fill from the parser.
- [ ] Rake task `oligos:backfill_amplicons[primer_set_id]` (default: all tiled amplicon sets), dry-run by default,
      prints a per-set summary: amplicons found, oligos unparsed, amplicons missing a side. Run on staging,
      review, then `APPLY=1`.
- [ ] Primer set edit form: show `locus` (labelled "Amplicon") and `category` for each oligo; add a set-level
      warning listing amplicons with no Forward or no Reverse.

## Task 2: Record alignment strand

- [ ] Migration: `oligo_alignment_positions.strand char(1)` (`+`/`-`), nullable.
- [ ] `lib/update_primers.sh`: derive strand from SAM flag 16 in the existing awk stage and insert it.
- [ ] Re-run alignment for all primer sets on staging (`update_primers.sh` with all ids) to populate it.
- [ ] Cross-check: for parsed Forward oligos strand should be `+`, Reverse `-`; report mismatches (they indicate
      a naming or alignment problem worth surfacing on the set page).

## Task 3: Organism parameters and result table

- [ ] Migration on `organisms`: `primer_three_prime_window int default 5`, `primer_max_mismatches int default 2`,
      `amplicon_dropout_pct float default 90`.
- [ ] Migration: `lineage_amplicon_coverage` (UNLOGGED, like `lineage_variant_primer_overlaps`):
      `organism_id, lineage_group_key, primer_set_id, amplicon, sequences int, working int,
      working_pct float, failing_forward int, failing_reverse int, top_failing_oligo_ids bigint[], created_at`.
      Unique on `(organism_id, lineage_group_key, primer_set_id, amplicon)`.
- [ ] Model + fixtures + basic test.

## Task 4: Coverage computation (`lib/visualization/amplicon_coverage.py`)

Inputs: the lineage group's filtered variant rows (same columns as `extract_all_variants.sh`), the group's full
sequence list, reference FASTA, and a primer table exported from the DB
(`oligo_id, primer_set_id, amplicon, side, sequence, taxon accession, ref_start, ref_end, strand`).

- [ ] For each primer, take the reference slice `[ref_start, ref_end)` and the primer sequence in reference
      orientation (revcomp if strand `-`). Precompute baseline mismatches primer-vs-reference — an ALT primer will
      mismatch the reference, which is expected.
- [ ] Index primers by interval; stream variant rows grouped by sequence. For each sequence, apply only the
      variants that intersect a primer to that primer's reference slice (SNPs substitute; deletions/insertions
      inside the primer → indel; deletion spanning the whole primer → primer can't bind) and score the primer
      with the rules above. Primers untouched by any variant keep their baseline result.
- [ ] Per sequence and amplicon: works = any Forward works and any Reverse works. Accumulate counts per
      (primer_set, amplicon); sequences with no variant rows contribute the baseline result.
- [ ] Output TSV for `store_amplicon_coverage.sh`.
- [ ] Tests (pytest, small synthetic reference): plain primer broken by 3′ SNP; ALT primer that carries the SNP
      rescues the amplicon; SNP outside 3′ window under `max_mismatches` still works; deletion inside primer
      fails; reverse primer 3′ end is at `ref_start` (not `ref_end`); degenerate base matches; sequence with no
      variants uses baseline; two ALTs with correlated variants (the case per-lineage frequencies get wrong).
- [ ] Performance check on the sars-cov-2 "all" group on staging: one pass over the variant stream, memory
      bounded by primers × amplicons, target well under the current `count_variants.sh` time.

## Task 5: Pipeline wiring

- [ ] `extract_all_variants.sh` (or a sibling query): emit sequence → lineage membership for all sequences in the
      date window, so groups include variant-free sequences.
- [ ] `update_visualization_data.sh`: after `count_variants.sh` for each lineage group, run
      `amplicon_coverage.py` for tiled amplicon primer sets; collect outputs; `store_amplicon_coverage.sh` replaces
      the organism's rows in one transaction (same pattern as `store_variant_overlaps.sh`).
- [ ] Python env: add to the existing pixi/conda env used by `update_visualization.nf`; no new deps beyond stdlib
      if possible.
- [ ] Fail loudly (non-zero exit) so `primer-monitor.sh` reports `FAILED`.

## Task 6: Primer set page

- [ ] `GET /primer_sets/:id/amplicon_coverage.json?lineage_group=…` → rows for that set and group.
- [ ] On the primer set page, for tiled amplicon sets, an amplicon table for the selected lineage: amplicon,
      working %, sequences, which oligo(s) carry it, which fail; rows under the organism's dropout threshold
      highlighted. Sort by working % ascending.
- [ ] Keep the existing per-oligo variant overlap view, but label ALT-rescued overlaps ("covered by …_ALT2")
      instead of presenting them as problems.
- [ ] System test with fixtures: an amplicon rescued by an ALT shows as working.

## Task 7: Notifications (follow-up, separate PR)

- [ ] Decide whether `identify_primers_for_notifications` should key on amplicon dropout for tiled amplicon sets
      rather than per-oligo overlap. Out of scope for the first release; note the gap on the page.

---

## Decisions (2026-10-11)

1. "Works" rule and defaults confirmed: 3′ window 5, ≤2 mismatches, indels fail, 90% dropout threshold.
2. Amplicon/side are inferred from oligo names, with a manual override in the primer set editor; upload does not
   require explicit columns.
3. The feature applies to tiled amplicon sets (the method formerly named "Multiplex PCR", renamed 2026-10-11).
   ARTICv3 (set 16) was misfiled as qPCR and is reclassified as
   Tiled amplicon (staging and production, 2026-10-11). Midnight (15), Varskip Short 1a (17) and Varskip Long (18),
   also tiling sets filed as qPCR, were reclassified in production only; staging is left as is.
