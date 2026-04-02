#!/usr/bin/env python3
"""
Replace compound ID columns (sample_name / SampleName / nvp) with inchi_key
in all relevant files under data/ and output/dr8.tsv.

Files where the NVP column is the ONLY compound identifier: rename + map values.
Files that already have an inchi_key column: drop the NVP column.
"""
import csv
import gzip
import io
import os

REPO = os.path.dirname(os.path.abspath(__file__))
MAP_FILE = os.path.join(REPO, 'ignore', 'nvp_inchi_map.tsv')

# Load mapping: compound_id -> inchi_key
mapping = {}
with open(MAP_FILE) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        mapping[row['nvp']] = row['inchi_key']

print(f"Loaded {len(mapping)} mappings")


def process_file(path, delimiter, old_col, new_col, already_has_inchi, is_gz=False):
    """
    Read file, transform compound ID column, write back in place.

    old_col: column name to replace/drop
    new_col: replacement column name (only used when not already_has_inchi)
    already_has_inchi: if True, just drop old_col; if False, rename + map
    """
    # Read
    open_fn = gzip.open if is_gz else open
    mode_r = 'rt'
    with open_fn(path, mode_r, encoding='utf-8-sig') as f:
        content = f.read()

    reader = csv.DictReader(io.StringIO(content), delimiter=delimiter)
    fieldnames = list(reader.fieldnames)
    rows = list(reader)

    if old_col not in fieldnames:
        print(f"  SKIP: '{old_col}' not found in {os.path.basename(path)}")
        return

    unmapped = set()

    if already_has_inchi:
        # Drop old_col
        new_fieldnames = [f for f in fieldnames if f != old_col]
        new_rows = [{k: v for k, v in row.items() if k != old_col} for row in rows]
    else:
        # Rename old_col -> new_col in place, mapping values
        idx = fieldnames.index(old_col)
        new_fieldnames = fieldnames[:idx] + [new_col] + fieldnames[idx+1:]
        new_rows = []
        for row in rows:
            val = row.pop(old_col)
            mapped = mapping.get(val, '') if val else ''
            if val and not mapped:
                unmapped.add(val)
            row[new_col] = mapped
            new_rows.append(row)

    if unmapped:
        print(f"  WARNING: {len(unmapped)} values not in mapping for {os.path.basename(path)}: "
              f"{sorted(unmapped)[:5]}{'...' if len(unmapped) > 5 else ''}")

    # Write back
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=new_fieldnames, delimiter=delimiter,
                            lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(new_rows)
    result = out.getvalue()

    mode_w = 'wt'
    with open_fn(path, mode_w, encoding='utf-8') as f:
        f.write(result)

    action = "dropped" if already_has_inchi else f"renamed '{old_col}' -> '{new_col}'"
    print(f"  OK: {os.path.basename(path)} ({action}, {len(new_rows)} rows)")


DATA = os.path.join(REPO, 'data')
OUT  = os.path.join(REPO, 'output')

# --- Files where sample_name must be renamed + mapped ---
for fname in [
    '8pt_singlecell_gfp.csv',
    '8pt_singlecell_gfp-gpi.csv',
    '8pt_well_gfp.csv',
    '8pt_well_gfp-gpi.csv',
]:
    process_file(os.path.join(DATA, fname),
                 delimiter=',', old_col='sample_name', new_col='inchi_key',
                 already_has_inchi=False)

# --- Gzipped files ---
for fname in [
    'gfp_well-level_correctedScoresOutput_data.csv.gz',
    'gfp-gpi_well-level_correctedScoresOutput_data.csv.gz',
]:
    process_file(os.path.join(DATA, fname),
                 delimiter=',', old_col='sample_name', new_col='inchi_key',
                 already_has_inchi=False, is_gz=True)

# --- Files that already have inchi_key: just drop the old column ---
process_file(os.path.join(DATA, 'hit_table_gfp.csv'),
             delimiter=',', old_col='SampleName', new_col='inchi_key',
             already_has_inchi=True)

process_file(os.path.join(DATA, 'hit_table_gfp-gpi.csv'),
             delimiter=',', old_col='SampleName', new_col='inchi_key',
             already_has_inchi=True)

process_file(os.path.join(DATA, 'moa_annotations.csv'),
             delimiter=',', old_col='nvp', new_col='inchi_key',
             already_has_inchi=True)

# --- hits_for_confirmation: SampleName only, rename + map ---
process_file(os.path.join(DATA, 'hits_for_confirmation.csv'),
             delimiter=',', old_col='SampleName', new_col='inchi_key',
             already_has_inchi=False)

# --- output/dr8.tsv: referenced by src/figure_2_and_s1.R ---
process_file(os.path.join(OUT, 'dr8.tsv'),
             delimiter='\t', old_col='nvp', new_col='inchi_key',
             already_has_inchi=False)

# --- output/ files ---

# moa.tsv: already has inchi_key, drop nvp
process_file(os.path.join(OUT, 'moa.tsv'),
             delimiter='\t', old_col='nvp', new_col='inchi_key',
             already_has_inchi=True)

# picks.tsv: only nvp column, map to inchi_key
process_file(os.path.join(OUT, 'picks.tsv'),
             delimiter='\t', old_col='nvp', new_col='inchi_key',
             already_has_inchi=False)

# primary.tsv: has nvp + inchi (full InChI string, different from inchi_key), map nvp -> inchi_key
process_file(os.path.join(OUT, 'primary.tsv'),
             delimiter='\t', old_col='nvp', new_col='inchi_key',
             already_has_inchi=False)

# erics_picks.tsv: compound column only
process_file(os.path.join(OUT, 'erics_picks.tsv'),
             delimiter='\t', old_col='compound', new_col='inchi_key',
             already_has_inchi=False)

# erics_picks_v2.tsv: compound column + other cols
process_file(os.path.join(OUT, 'erics_picks_v2.tsv'),
             delimiter='\t', old_col='compound', new_col='inchi_key',
             already_has_inchi=False)

print("\nDone.")
