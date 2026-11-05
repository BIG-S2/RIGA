#!/usr/bin/env python
"""benchmark_harness_model_ablation.py (2026-10-13, review point 2) -- progressive-model ablation on the 1,000 permuted-null
phenotypes of the big null (same Y* files, so every arm is paired with the published strat_project / pooled_project results).
Thin wrapper around benchmark_harness.py: it registers new arms and extends build_frames; the harness is not modified.

Identical participants (S-rows = every subject of every project x sequencing center cell with both labels):
  abl_addcoll_S    Age + Sex + PC1-20 + center + Study_cov  (single-label or < 50 projects collapsed to Other, as pooled_project)
  abl_addfull_S    Age + Sex + PC1-20 + center + Study      (every project its own level)
  strat_project    Age + Sex + PC1-20 + cell                (the stratified design; existing arm)
  abl_int_cls_S    strat_project + latent class + posterior + entropy
All informative participants of an additive model (rows whose project and sequencing center both hold both labels):
  pooled_project   the published pooled comparator (existing arm, all rows)
  (an all-informative-rows additive arm was dropped: its rows differ from the S-rows by 10 subjects)
Study-stratified meta-analysis:
  meta_project     per-project logistic null on that project's mixed-cell rows (Age + Sex + PC1-10, PC1-4 below 200
                   subjects, cell dummies when the project spans more than one sequencing center), fixed-effect score
                   meta-analysis U = sum U_k, V = sum V_k; projects with < 50 subjects or < 5 of either label are dropped.
Usage: as benchmark_harness.py, with --arms drawn from the list above."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import qr as scipy_qr
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_harness as P  # noqa: E402

H = P.H
NEW_UNCOND = ["abl_addcoll_S", "abl_addfull_S", "abl_int_cls_S"]
NEW_CLR = ["meta_project"]
for a in NEW_UNCOND:
    P.UNCOND_ARMS.append(a)
    P.ARM_ORDER.append(a)
P.CLR_ARMS.append("meta_project")
P.ARM_ORDER.append("meta_project")
P.ARM_LABEL.update(abl_addcoll_S="A1", abl_addfull_S="A2", abl_int_cls_S="A4", meta_project="MA")


class MetaProject:
    """Fixed-effect score meta-analysis over per-project logistic nulls (same interface as P.CLR)."""

    def __init__(self, pool: pd.DataFrame, y: np.ndarray, min_n: int = 50, min_label: int = 5):
        base = P.base_frame(pool, y)
        cell = (base["Study"].astype(str) + "|" + base["center"].astype(str)).to_numpy()
        mixed = np.zeros(len(base), dtype=bool)
        mixed[P.mixed_rows(cell, y)] = True
        pcs10 = [c for c in H.PCS10 if c in base.columns]
        self.parts, self.dropped = [], 0
        for proj, idx in base.groupby("Study").indices.items():
            idx = np.asarray(idx)[mixed[np.asarray(idx)]]
            yk = y[idx]
            n1 = int(yk.sum()); n0 = int(idx.size - n1)
            if idx.size < min_n or n1 < min_label or n0 < min_label:
                self.dropped += int(idx.size)
                continue
            d = base.iloc[idx].reset_index(drop=True)
            d["cell_k"] = cell[idx]
            covs = ["Age", "Sex"] + (pcs10 if idx.size >= 200 else pcs10[:4])
            if d["cell_k"].nunique() > 1:
                covs = covs + ["cell_k"]
            nf = H.null_fit_local(d, covs)
            self.parts.append((idx, nf))
        self.rows = np.concatenate([p[0] for p in self.parts]).astype(np.int64)
        self.codes = np.concatenate([np.full(p[0].size, k) for k, p in enumerate(self.parts)]).astype(np.int64)
        offs = np.cumsum([0] + [p[0].size for p in self.parts])
        self.slices = [slice(int(offs[i]), int(offs[i + 1])) for i in range(len(self.parts))]
        self.n_rows = int(self.rows.size)
        self.n_cases = int(y[self.rows].sum())
        self.n_strata = len(self.parts)
        n1, n0 = self.n_cases, self.n_rows - self.n_cases
        self.n_effective = float(4.0 * n1 * n0 / max(self.n_rows, 1))
        self.info_weight = float(n1 * n0 / max(self.n_rows, 1))

    def score(self, G: np.ndarray):
        U = np.zeros(G.shape[0]); V = np.zeros(G.shape[0])
        for (_, nf), sl in zip(self.parts, self.slices):
            Gk = G[:, sl].astype(np.float64)
            U += Gk @ nf["e"]
            Gt = Gk * nf["sqrt_w"][None, :]
            GQ = Gt @ nf["Q"]
            V += np.sum(Gt * Gt, axis=1) - np.sum(GQ * GQ, axis=1)
        V = np.where(V > 1e-12, V, np.nan)
        stat = U * U / V
        return chi2.sf(stat, df=1), stat, V


_orig_build_frames = P.build_frames


def full_rank_design(d: pd.DataFrame, numeric: list, factors: list):
    """Explicit treatment dummies for the factors, then drop dummy columns aliased with earlier columns (QR with pivoting),
    so that uncollapsed project + sequencing center (projects nested in one center are aliased with that center) is fit
    as a full-rank logistic model. Returns the frame with dummy columns added and the covariate list."""
    d = d.copy(); dums = []
    for f in factors:
        D = pd.get_dummies(d[f].astype(str), prefix=f"lvl_{f}", drop_first=True, dtype=float)
        d = pd.concat([d, D], axis=1); dums += list(D.columns)
    cols = numeric + dums
    M = np.column_stack([np.ones(len(d))] + [d[c].to_numpy(float) for c in cols])
    _, R, piv = scipy_qr(M, pivoting=True, mode="economic")
    rank = int(np.sum(np.abs(np.diag(R)) > 1e-8 * np.abs(R[0, 0])))
    keep_idx = set(int(i) for i in piv[:rank])
    kept = [c for j, c in enumerate(cols, start=1) if j in keep_idx or c in numeric]
    P.log(f"  full_rank_design: {len(dums)} dummies, {len(cols) + 1 - rank} aliased dropped")
    return d, kept


def build_frames(pool, y, matched, py_sets, arms):
    base_arms = [a for a in arms if a not in NEW_UNCOND + NEW_CLR]
    frames, clr = _orig_build_frames(pool, y, matched, py_sets, base_arms)
    base = P.base_frame(pool, y)
    pcs20 = [c for c in H.PCS20 if c in base.columns]
    cell = (base["Study"].astype(str) + "|" + base["center"].astype(str)).to_numpy()
    s_rows = P.mixed_rows(cell, y)
    if any(a in arms for a in ("abl_addcoll_S", "abl_addfull_S", "abl_int_cls_S")):
        d = base.iloc[s_rows].reset_index(drop=True)
        if "abl_addcoll_S" in arms:
            dc = d.copy(); dc["Study_cov"] = H.collapse_levels(dc["Study"], dc["Y"].to_numpy())
            frames["abl_addcoll_S"] = (dc, s_rows, ["Age", "Sex"] + pcs20 + ["center", "Study_cov"], None)
        if "abl_addfull_S" in arms:
            df_, covs_ = full_rank_design(d, ["Age", "Sex"] + pcs20, ["center", "Study_strat"])
            frames["abl_addfull_S"] = (df_, s_rows, covs_, None)
        if "abl_int_cls_S" in arms:
            di = d.copy(); di["cell_SC_strat"] = cell[s_rows]
            frames["abl_int_cls_S"] = (di, s_rows, ["Age", "Sex"] + pcs20 + ["cell_SC_strat", "lpa_cluster", "cluster_prob", "entropy_norm"], None)
    if "meta_project" in arms:
        clr["meta_project"] = MetaProject(pool, y)
        P.log(f"  [meta_project] projects={clr['meta_project'].n_strata} rows={clr['meta_project'].n_rows} dropped={clr['meta_project'].dropped}")
    return frames, clr


P.build_frames = build_frames

if __name__ == "__main__":
    sys.exit(P.main())
