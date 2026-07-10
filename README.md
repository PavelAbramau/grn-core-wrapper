# grn-core-wrapper

A minimal, self-contained Snakemake pipeline containing **only** the 6
GRN-core rules from scenicplus's official Snakemake pipeline:

- `tf_to_gene`
- `region_to_gene`
- `eGRN_direct`
- `eGRN_extended`
- `AUCell_direct`
- `AUCell_extended`

Every rule's shell command is copied verbatim from the officially installed
scenicplus package's own Snakefile (see "Provenance" below) -- nothing about
*how* these 6 steps run has been changed. What's different is that this
repo replaces every upstream input those 6 rules would normally depend on
with small, frozen, committed fixture files (`data/fixtures/`), so this
pipeline never needs network access, a real motif database, or any of the
slower/flakier upstream steps to run. It exists to let CI (and anyone else)
verify, cheaply and deterministically, that the GRN-core computation itself
still behaves as expected -- see `tests/compare_to_control.py`.

## Provenance: pinned scenicplus version

The fixtures and control comparison files in this repo were generated with
this exact scenicplus build (also recorded in `config/config.yaml`):

| | |
|---|---|
| pip version (`pip show scenicplus`) | `1.0a2` |
| git commit | `840dab85de3044846234c157e327b6ea0290abfb` |
| git describe | `v1.0a2-2-g840dab8` (2 commits past the `v1.0a2` tag) |
| commit date | 2026-01-16 14:46:10 +0100 |

scenicplus is not published on PyPI; it was pip-installed with `--no-deps`
from a local checkout of the above commit. See `../environment-fixes/` (in
the parent project) for the full, verified install sequence, and
`requirements-lock.txt` in this repo for the exact `pip freeze` snapshot of
the environment this was validated against.

A different scenicplus commit could change `tf_to_gene`/`region_to_gene`/
`eGRN`/`AUCell` behavior (argument names, defaults, internal algorithms) and
invalidate the comparison in `tests/compare_to_control.py` -- if you bump
the pinned version, re-run the full Phase-1 pipeline
(`scenicplus-validation-harness/`) to regenerate fresh fixtures and control
files before trusting this repo's tests again.

## What this repo excludes

This Snakefile deliberately does **not** contain 7 of the 13 rules in
scenicplus's official pipeline. Their outputs are supplied as static files
under `data/fixtures/` instead. If you need to actually run one of these
steps (e.g. against your own real data), here's what produces it upstream:

| Excluded rule | What produces this output upstream |
|---|---|
| `download_genome_annotations` | `scenicplus prepare_data download_genome_annotations` (`scenicplus.data_wrangling.gene_search_space.download_gene_annotation_and_chromsizes`) -- biomart + NCBI eutils; writes `genome_annotation.tsv` + `chromsizes.tsv`. |
| `prepare_GEX_ACC_multiome` | `scenicplus prepare_data prepare_GEX_ACC` (`scenicplus.data_wrangling.adata_cistopic_wrangling.process_multiome_data`) -- merges a cisTopic object + GEX AnnData into `combined_GEX_ACC_mudata` (our `data/fixtures/combined_GEX_ACC_mudata.h5mu`). |
| `get_search_space` | `scenicplus prepare_data search_spance` [sic -- upstream's own CLI subcommand name] (`scenicplus.data_wrangling.gene_search_space.get_search_space`) -- writes `search_space.tsv` (our `data/fixtures/search_space.tsv`). |
| `motif_enrichment_cistarget` | `scenicplus grn_inference motif_enrichment_cistarget` (`pycistarget.motif_enrichment_cistarget.cisTarget`) -- writes `ctx_results.hdf5`/`.html`. |
| `motif_enrichment_dem` | `scenicplus grn_inference motif_enrichment_dem` (`pycistarget.motif_enrichment_dem.DEM`) -- writes `dem_results.hdf5`/`.html`. |
| `prepare_menr` | `scenicplus prepare_data prepare_menr` (`scenicplus.data_wrangling.cistarget_wrangling.get_and_merge_cistromes`) -- merges the two motif-enrichment results above into `cistromes_direct.h5ad`/`cistromes_extended.h5ad`/`tf_names.txt` (our `data/fixtures/` copies of these). |
| `scplus_mudata` | `scenicplus grn_inference create_scplus_mudata` (`scenicplus.scenicplus_mudata.ScenicPlusMuData`) -- combines this pipeline's own `AUCell_direct`/`AUCell_extended`/`eRegulon_direct`/`eRegulons_extended` outputs with the GEX/ACC MuData into one final `scplusmdata.h5mu`. Not run here since it consumes this repo's own outputs, not an upstream dependency -- see `scenicplus-validation-harness/` for the full pipeline including this step.

For the full, unmodified 13-rule pipeline (including all of the above), see
`../scenicplus-validation-harness/` in the parent project -- that's also
where `data/fixtures/` and `results_control/` in this repo were copied from.

## Self-containment

- `data/fixtures/` and `results_control/` are **real file copies**, not
  symlinks -- verified with `find . -type l` (see below). This repo does
  not depend on `../scenicplus-validation-harness/` or `../results_control/`
  existing at runtime; those are only where the fixtures were originally
  produced.
- `requirements-lock.txt` + `config/config.yaml`'s pinned scenicplus commit
  fully document the environment this was validated against.
- `environment/` holds its own copies of `environment.yml` and
  `requirements.txt` (not a reference to `../environment-fixes/`), so the
  install sequence works even if this directory is ever extracted into its
  own standalone repo.
- CI (`.github/workflows/ci.yml`) runs entirely off what's committed in this
  repo -- it never invokes the excluded rules above, never hits the
  network for genome annotations or motif databases (scenicplus itself is
  installed pinned to the exact commit above, directly from GitHub), and
  never needs the parent `scenicplus-validation-harness/` project to exist.

```
$ find . -type l
# (no output -- nothing in this repo is a symlink)
```

## Running

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --directory results --cores 2
python tests/compare_to_control.py
```

## Layout

```
grn-core-wrapper/
  workflow/Snakefile          # the 6 rules, copied verbatim
  config/config.yaml          # points input_data at data/fixtures/, pins scenicplus version
  data/fixtures/               # frozen inputs (see table above)
  environment/                 # own copy of environment.yml + requirements.txt
  results/                     # this repo's own run output (gitignored)
  results_control/             # frozen control comparison files (committed, small)
  tests/compare_to_control.py  # exact-match test against results_control/
  requirements-lock.txt        # pip freeze snapshot
  .github/workflows/ci.yml
```
