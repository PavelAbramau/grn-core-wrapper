#!/usr/bin/env python
"""Compare the pruned 6-rule GRN-core Snakefile's outputs against the
control: the same eGRN/AUCell outputs produced by the full, unmodified
official scenicplus Snakemake pipeline (Phase 1) on the same toy dataset.

Both runs consume the exact same frozen fixtures (grn-core-wrapper/data/
fixtures/, copied from scenicplus-validation-harness/run/) for
tf_to_gene/region_to_gene/eGRN/AUCell -- so if the pruned Snakefile's copied
rule definitions and config are faithful to the original, re-running just
those 6 steps should reproduce the Phase-1 run's own eGRN_direct/extended
and AUCell_direct/extended outputs.

Prints a full diff on any mismatch. Does not adjust tolerances to force a
pass -- if something doesn't match, that's the answer.
"""
import sys
from pathlib import Path

import mudata
import numpy as np
import pandas as pd

GRN_CORE_ROOT = Path(__file__).resolve().parent.parent

RESULTS = GRN_CORE_ROOT / "results"
# Committed, real-copy control files living inside this repo (see
# README.md "Self-containment") -- NOT the parent project's
# ../results_control/, so this repo (and CI) never depends on anything
# outside itself.
CONTROL = GRN_CORE_ROOT / "results_control"

EREGULON_PAIRS = [
    ("direct", RESULTS / "eRegulon_direct.tsv", CONTROL / "snakemake_09a_eregulons_direct.tsv"),
    ("extended", RESULTS / "eRegulons_extended.tsv", CONTROL / "snakemake_09b_eregulons_extended.tsv"),
]
AUCELL_PAIRS = [
    ("direct", RESULTS / "AUCell_direct.h5mu", CONTROL / "snakemake_10a_auc_direct.h5mu"),
    ("extended", RESULTS / "AUCell_extended.h5mu", CONTROL / "snakemake_10b_auc_extended.h5mu"),
]

AUC_ATOL = 1e-6

failures: list[str] = []


def fail(msg: str) -> None:
    failures.append(msg)
    print(f"[FAIL] {msg}")


def ok(msg: str) -> None:
    print(f"[OK]   {msg}")


def compare_eregulons(label: str, test_path: Path, control_path: Path) -> None:
    print(f"\n=== eRegulons ({label}) ===")
    print(f"  test:    {test_path}")
    print(f"  control: {control_path}")
    test_df = pd.read_csv(test_path, sep="\t")
    control_df = pd.read_csv(control_path, sep="\t")

    test_tfs = set(test_df["TF"].unique())
    control_tfs = set(control_df["TF"].unique())
    if test_tfs != control_tfs:
        fail(
            f"eRegulons ({label}): TF sets differ.\n"
            f"    only in test:    {sorted(test_tfs - control_tfs)}\n"
            f"    only in control: {sorted(control_tfs - test_tfs)}"
        )
    else:
        ok(f"eRegulons ({label}): same {len(test_tfs)} TFs: {sorted(test_tfs)}")

    common_tfs = test_tfs & control_tfs
    test_targets = test_df.groupby("TF")["Gene"].apply(lambda s: set(s)).to_dict()
    control_targets = control_df.groupby("TF")["Gene"].apply(lambda s: set(s)).to_dict()
    for tf in sorted(common_tfs):
        t_targets = test_targets.get(tf, set())
        c_targets = control_targets.get(tf, set())
        if len(t_targets) != len(c_targets):
            fail(
                f"eRegulons ({label}): TF {tf} target gene COUNT differs: "
                f"test={len(t_targets)} control={len(c_targets)}"
            )
        elif t_targets != c_targets:
            fail(
                f"eRegulons ({label}): TF {tf} target gene SET differs (same count, different genes).\n"
                f"    only in test:    {sorted(t_targets - c_targets)}\n"
                f"    only in control: {sorted(c_targets - t_targets)}"
            )
        else:
            ok(f"eRegulons ({label}): TF {tf} -- {len(t_targets)} target genes, identical")

    if len(test_df) != len(control_df):
        fail(f"eRegulons ({label}): row count differs: test={len(test_df)} control={len(control_df)}")
    else:
        ok(f"eRegulons ({label}): same row count ({len(test_df)})")


def compare_aucell(label: str, test_path: Path, control_path: Path) -> None:
    print(f"\n=== AUCell ({label}) ===")
    print(f"  test:    {test_path}")
    print(f"  control: {control_path}")
    test_mdata = mudata.read(str(test_path))
    control_mdata = mudata.read(str(control_path))

    for modality in ["Gene_based", "Region_based"]:
        if modality not in test_mdata.mod or modality not in control_mdata.mod:
            fail(f"AUCell ({label}): modality {modality!r} missing from test or control MuData")
            continue

        t_ad = test_mdata[modality]
        c_ad = control_mdata[modality]

        t_cells, c_cells = list(t_ad.obs_names), list(c_ad.obs_names)
        if set(t_cells) != set(c_cells):
            fail(f"AUCell ({label}/{modality}): cell barcode sets differ ({len(t_cells)} vs {len(c_cells)})")
            continue

        t_eregs, c_eregs = list(t_ad.var_names), list(c_ad.var_names)
        if set(t_eregs) != set(c_eregs):
            fail(
                f"AUCell ({label}/{modality}): eRegulon column sets differ.\n"
                f"    only in test:    {sorted(set(t_eregs) - set(c_eregs))}\n"
                f"    only in control: {sorted(set(c_eregs) - set(t_eregs))}"
            )
            continue

        # Align both to the same cell/eRegulon order before comparing values.
        c_ad_aligned = c_ad[t_cells, t_eregs]
        t_x = np.asarray(t_ad.X.todense() if hasattr(t_ad.X, "todense") else t_ad.X)
        c_x = np.asarray(c_ad_aligned.X.todense() if hasattr(c_ad_aligned.X, "todense") else c_ad_aligned.X)

        diff = np.abs(t_x - c_x)
        max_diff = float(diff.max()) if diff.size else 0.0
        n_mismatch = int((diff > AUC_ATOL).sum())

        if n_mismatch > 0:
            mismatch_idx = np.argwhere(diff > AUC_ATOL)
            sample = mismatch_idx[:10]
            sample_lines = [
                f"      cell={t_cells[r]!r} eRegulon={t_eregs[c]!r} test={t_x[r, c]:.8f} control={c_x[r, c]:.8f} diff={diff[r, c]:.2e}"
                for r, c in sample
            ]
            fail(
                f"AUCell ({label}/{modality}): {n_mismatch}/{diff.size} values exceed atol={AUC_ATOL} "
                f"(max diff={max_diff:.2e}). First {len(sample)} mismatches:\n" + "\n".join(sample_lines)
            )
        else:
            ok(f"AUCell ({label}/{modality}): all {diff.size} values match within atol={AUC_ATOL} (max diff={max_diff:.2e})")


def main() -> int:
    for path in [p for _, a, b in EREGULON_PAIRS for p in (a, b)] + [p for _, a, b in AUCELL_PAIRS for p in (a, b)]:
        if not path.exists():
            print(f"[FAIL] missing file: {path}")
            return 1

    for label, test_path, control_path in EREGULON_PAIRS:
        compare_eregulons(label, test_path, control_path)
    for label, test_path, control_path in AUCELL_PAIRS:
        compare_aucell(label, test_path, control_path)

    print("\n" + "=" * 70)
    if failures:
        print(f"RESULT: FAIL -- {len(failures)} mismatch(es) found (see [FAIL] lines above)")
        return 1
    print("RESULT: PASS -- pruned 6-rule Snakefile reproduces the Phase-1 control run exactly")
    return 0


if __name__ == "__main__":
    sys.exit(main())
