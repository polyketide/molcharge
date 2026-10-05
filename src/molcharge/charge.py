#!/usr/bin/env python3
"""molcharge.charge — net charge Z of a small molecule versus pH, from a SMILES. Local, offline.

Engine: Dimorphite-DL 2.0.2 rule table (41 ionisable moieties, mean pKa ± σ; Apache-2.0) → independent-site
Henderson–Hasselbalch sum, with a ±precision·σ envelope. Every output block names the engine, its version and the
rule-table hash: a predicted pKa is a function of (molecule, method), never a property of the molecule alone.

    Z(pH) = q0 + Σ_base 1/(1+10^(pH−pKa)) − Σ_acid 1/(1+10^(pKa−pH))

q0 = permanent formal charge (sulfonium, quaternary N, ...). Invariants used as self-test gates: Z is monotone
non-increasing in pH, Z(pH→−∞) = q0 + N_base and Z(pH→+∞) = q0 − N_acid. Monotonicity and the limits cannot detect
a double-counted site, so site counting has its own gates.

Three corrections over calling Dimorphite-DL directly (each has a self-test gate):
  * one site per site: Dimorphite returns one record per SMARTS hit, so a secondary/tertiary amine N was counted
    once per C–N bond (triethylamine Z@7.4 +2.56 instead of +0.85) and H3PO4 got six terms for two ionisations;
  * no silent truncation: Dimorphite's default cap is 50 hits per pattern, which cut molecules with more than
    ~25 aliphatic N to 25 base sites; the cap is raised to 1000 and reaching it is reported;
  * table order: Phosphonate_ester is matched right after Phosphonate. In file order the amide/aniline/amine
    patterns took the N of a P(=O)(OR)(NR'R'')–OH first and the P–OH was lost (agrocin 84 Z@7.4 −0.995 → −1.996).

Usage:
  molcharge --smiles "NCC(=O)O" [--ph 0:14:0.1] [--at 7.4,8.0] [--precision 1.0] [--out STEM] [--json]
  molcharge --selftest
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys

LN10 = math.log(10.0)
ENGINES = ("dimorphite",)
VERSION = "0.1.0"


# ----------------------------------------------------------------------------------------------
# identity block — canonical SMILES, InChIKey, formula, MW, undefined-stereocentre warning (RDKit)
# ----------------------------------------------------------------------------------------------
def identity_block(smiles: str) -> dict:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors, inchi
    except Exception as e:  # pragma: no cover
        return {"ok": False, "error": f"rdkit unavailable: {e}"}
    m = Chem.MolFromSmiles(smiles or "")
    if m is None:
        return {"ok": False, "valid": False, "input": smiles, "note": "invalid SMILES — flagged, not coerced"}
    out = {"ok": True, "valid": True, "canonical_smiles": Chem.MolToSmiles(m),
           "formula": rdMolDescriptors.CalcMolFormula(m), "mol_weight": round(Descriptors.MolWt(m), 2),
           "n_heavy_atoms": m.GetNumHeavyAtoms(), "identity_source": "rdkit"}
    try:
        key = inchi.MolToInchiKey(m) if inchi.INCHI_AVAILABLE else None
        out["inchikey"] = key
        out["inchikey_skeleton"] = key.split("-")[0] if key else None
    except Exception:
        out["inchikey"] = None
    try:
        unspec = [i for i, t in Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
                  if t == "?"]
        out["n_stereocentres_undefined"] = len(unspec)
        if unspec:
            out["stereo_warning"] = f"{len(unspec)} undefined stereocentre(s) at atom idx {unspec}"
    except Exception:
        pass
    return out


# ----------------------------------------------------------------------------------------------
# grid / helpers
# ----------------------------------------------------------------------------------------------
def parse_grid(spec: str) -> tuple[float, float, float]:
    parts = [p.strip() for p in str(spec).split(":")]
    if len(parts) == 2:
        parts.append("0.1")
    if len(parts) != 3:
        raise ValueError(f"--ph must be lo:hi[:step], got {spec!r}")
    lo, hi, step = (float(p) for p in parts)
    if not (step > 0):
        raise ValueError(f"step must be > 0, got {step}")
    if not (hi > lo):
        raise ValueError(f"need hi > lo, got {lo}:{hi}")
    if (hi - lo) / step > 100000:
        raise ValueError("grid too fine (>100000 points)")
    return lo, hi, step


def make_grid(lo: float, hi: float, step: float) -> list[float]:
    n = int(round((hi - lo) / step))
    pts = [round(lo + i * step, 10) for i in range(n + 1)]
    if pts[-1] > hi + 1e-9:
        pts = pts[:-1]
    return pts


def _frac(x: float) -> float:
    """1/(1+10^x) with overflow safety (x may be ±1000 for Dimorphite's nitro sentinel)."""
    if x > 300:
        return 0.0
    if x < -300:
        return 1.0
    return 1.0 / (1.0 + 10.0 ** x)


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


# ----------------------------------------------------------------------------------------------
# Engine A — Dimorphite-DL rule table → site list → HH sum
# ----------------------------------------------------------------------------------------------
def dimorphite_info() -> dict:
    try:
        import dimorphite_dl
        from importlib.resources import files
        table = (files("dimorphite_dl.smarts") / "site_substructures.smarts").read_text(encoding="utf-8")
        import rdkit
        return {"available": True, "version": getattr(dimorphite_dl, "__version__", "?"),
                "table_sha256": sha256_text(table), "rdkit_version": rdkit.__version__,
                "n_table_lines": sum(1 for ln in table.splitlines() if ln.strip() and not ln.startswith("#"))}
    except Exception as e:
        return {"available": False, "error": f"{type(e).__name__}: {e}",
                "hint": "pip install dimorphite-dl==2.0.2"}


MAX_HITS_PER_PATTERN = 1000   # Dimorphite's own upper bound. Its default (50) is per PATTERN and counts hits, not sites:
                              # an internal sec. amine N is hit once per C–N bond, so >25 aliphatic N were silently cut
                              # to 25 (measured 2026-10-05 on linear oligo-ethyleneimines). Reaching the cap is reported.


# Dimorphite matches its table in FILE order and protects every matched atom. Phosphonate_ester (P(=O)(OR)(N/C)–OH,
# file position 38) sits after the amide/aniline/amine patterns, so when the N on P also carries a carbon those
# patterns take the N first and the whole P–OH match is discarded (agrocin 84 read Z@7.4 −0.995; −1.996 with this
# move, chemical judgement ≈ −2). Measured 2026-10-05. The table itself is untouched; only the iteration order moves.
TABLE_REORDER = (("Phosphonate_ester", "Phosphonate"),)   # (pattern, placed right after)


def _site_detector(**kw):
    from dimorphite_dl.protonate.detect import ProtonationSiteDetector

    class _Reordered(ProtonationSiteDetector):
        def _iterate_available_substructures(self):
            subs = list(super()._iterate_available_substructures())
            for moved, anchor in TABLE_REORDER:
                names = [x.name for x in subs]
                if moved not in names or anchor not in names:   # never fall back to the old order silently
                    raise RuntimeError(f"Dimorphite table has no {moved!r} or {anchor!r}; TABLE_REORDER needs review")
                item = subs.pop(names.index(moved))
                subs.insert([x.name for x in subs].index(anchor) + 1, item)
            yield from subs
    return _Reordered(**kw)


def detect_sites_dimorphite(smiles: str) -> dict:
    """Dimorphite-DL 2.x site detection. Returns the neutral reference state, q0, and one entry per DISTINCT site.

    kind rule mirrors dimorphite_dl.protonate.change.set_protonation_charge: a NITROGEN site whose moiety name
    does not start with '*' is protonated as +1 / deprotonated as 0 (a base); every other site (O, S, and the
    '*'-marked acidic N–H moieties) is protonated as 0 / deprotonated as −1 (an acid).

    Dedup (fixed 2026-10-05). find_sites() returns one record per SMARTS HIT, and several hits can name the same
    site. Counting each hit as a site adds spurious Henderson–Hasselbalch terms — measured before the fix:
    diethylamine Z@7.4 +1.70 and triethylamine +2.56 (ethylamine +0.85), phosphoric acid six terms for two
    ionizations. Two levels, because neither alone is enough:
      (1) the same pattern on the same set of HEAVY atoms is one group, collapsing to its first hit, when
          (a) the pattern carries ONE pKa and the hits name the same site atom — they differ only in which H
              was matched (an NH2 amide or sulfonamide); or
          (b) the pattern carries SEVERAL pKa (Phosphate, Phosphonate) — the hits are role permutations of one
              group (the three OH of phosphoric acid). HH depends only on the multiset of pKa, so which
              equivalent O carries which value does not matter.
          A single-pKa pattern whose hits name DIFFERENT site atoms is not merged even on one heavy-atom set:
          in 3- and 4-membered rings two distinct N–H sites can share one (diaziridinone O=C1NN1).
      (2) the same (site atom, kind, pKa) is one site: a secondary/tertiary amine N is hit once per C–N bond,
          each hit with a different carbon partner (so level 1 does not merge them) but the same site.
    The number of Henderson–Hasselbalch TERMS removed (not hits) is reported as duplicate_terms_collapsed.
    Measured on the 10,724 unique molecules of the local IUPAC high-confidence + SAMPL6 sets: 0 repeated
    (atom, kind, pKa), 0 atoms with two pKa of one kind, 0 site atoms lost."""
    from rdkit import Chem
    from dimorphite_dl.mol import MoleculeRecord
    rec = MoleculeRecord(smiles)
    rec, sites = _site_detector(max_sites_per_molecule=MAX_HITS_PER_PATTERN).find_sites(rec)
    per_pattern: dict[str, int] = {}
    for s in sites:
        per_pattern[s.name] = per_pattern.get(s.name, 0) + len(s.pkas)
    truncated = sorted(n for n, c in per_pattern.items() if c >= MAX_HITS_PER_PATTERN - 1)   # Dimorphite stops at cap−1
    mol = rec.mol                                   # prepared: neutralised, explicit Hs, heavy-atom indices preserved
    if mol is None:
        raise RuntimeError("Dimorphite-DL could not prepare the molecule")
    q0 = int(sum(a.GetFormalCharge() for a in mol.GetAtoms()))
    fixed = [{"atom_idx": a.GetIdx(), "element": a.GetSymbol(), "charge": a.GetFormalCharge()}
             for a in mol.GetAtoms() if a.GetFormalCharge() != 0]
    entries, seen_groups, seen_sites, collapsed = [], set(), set(), 0
    for s in sites:
        heavy = frozenset(i for i in s.idxs_match if mol.GetAtomWithIdx(i).GetAtomicNum() > 1)
        named = frozenset(s.idxs_match[pk.idx_site] for pk in s.pkas) if len(s.pkas) == 1 else None
        group = (s.name, heavy, named)              # named=None: multi-pKa pattern, role permutations merge
        if group in seen_groups:                    # level 1: H choice or role permutation of the same group
            collapsed += len(s.pkas)
            continue
        seen_groups.add(group)
        for pk in s.pkas:
            idx = s.idxs_match[pk.idx_site]
            atom = mol.GetAtomWithIdx(idx)
            is_base = atom.GetAtomicNum() == 7 and not s.name.startswith("*")
            site = (idx, "base" if is_base else "acid", round(pk.mean, 3))
            if site in seen_sites:                  # level 2: another hit naming the same site
                collapsed += 1
                continue
            seen_sites.add(site)
            entries.append({
                "atom_idx": idx, "element": atom.GetSymbol(), "moiety": s.name.lstrip("*"),
                "kind": "base" if is_base else "acid",
                "pka": round(pk.mean, 3), "sigma": round(pk.stdev, 3),
                "q_protonated": 1 if is_base else 0, "q_deprotonated": 0 if is_base else -1,
                "smarts": s.smarts,
            })
    entries.sort(key=lambda e: (e["atom_idx"], e["pka"]))
    ref = Chem.MolToSmiles(Chem.RemoveHs(mol))
    return {"reference_state_smiles": ref, "q0": q0, "fixed_charge_atoms": fixed, "sites": entries,
            "n_acid": sum(1 for e in entries if e["kind"] == "acid"),
            "n_base": sum(1 for e in entries if e["kind"] == "base"),
            "duplicate_terms_collapsed": collapsed, "truncated_patterns": truncated}


def hh_charge(sites: list[dict], q0: int, ph: float, shift: float = 0.0) -> float:
    z = float(q0)
    for e in sites:
        pka = e["pka"] + shift
        if e["kind"] == "base":
            z += _frac(ph - pka)          # fraction protonated (+1)
        else:
            z -= _frac(pka - ph)          # fraction deprotonated (−1)
    return z


def find_pI(zf, lo: float, hi: float, tol: float = 1e-6) -> dict:
    """Unique zero-crossing of a monotone non-increasing Z(pH) on [lo, hi] by bisection."""
    zlo, zhi = zf(lo), zf(hi)
    if zlo < 0 and zhi < 0:
        return {"pI": None, "note": "Z < 0 over the whole range (net negative everywhere; no isoelectric point in range)"}
    if zlo > 0 and zhi > 0:
        return {"pI": None, "note": "Z > 0 over the whole range (net positive everywhere; no isoelectric point in range)"}
    if abs(zlo) < 1e-12:
        return {"pI": lo, "note": "Z = 0 at the lower bound"}
    if abs(zhi) < 1e-12:
        return {"pI": hi, "note": "Z = 0 at the upper bound"}
    a, b = lo, hi
    for _ in range(200):
        m = 0.5 * (a + b)
        zm = zf(m)
        if abs(zm) < 1e-12 or (b - a) < tol:
            return {"pI": round(m, 4), "note": "bisection on the monotone curve"}
        if zm > 0:
            a = m
        else:
            b = m
    return {"pI": round(0.5 * (a + b), 4), "note": "bisection (iteration cap)"}


def engine_dimorphite(smiles: str, grid: list[float], at: list[float], precision: float) -> dict:
    info = dimorphite_info()
    if not info.get("available"):
        return {"ok": False, "engine": "dimorphite", "error": "engine not installed", **info}
    try:
        det = detect_sites_dimorphite(smiles)
    except Exception as e:
        return {"ok": False, "engine": "dimorphite", "error": f"site detection failed: {type(e).__name__}: {e}",
                "note": "flag-not-fabricate: no curve is emitted for a molecule the rule table cannot prepare"}
    sites, q0 = det["sites"], det["q0"]
    zf = lambda ph, sh=0.0: hh_charge(sites, q0, ph, sh)  # noqa: E731
    # envelope: every term is monotone in its pKa, so shifting all pKa by ±precision·σ bounds Z
    z = [zf(p) for p in grid]
    z_hi = [hh_charge([dict(e, pka=e["pka"] + precision * e["sigma"]) for e in sites], q0, p) for p in grid]
    z_lo = [hh_charge([dict(e, pka=e["pka"] - precision * e["sigma"]) for e in sites], q0, p) for p in grid]
    curve = [{"ph": p, "z": round(a, 6), "z_lo": round(b, 6), "z_hi": round(c, 6)} for p, a, b, c in zip(grid, z, z_lo, z_hi)]
    pi = find_pI(zf, grid[0], grid[-1])
    z_at = [{"ph": p, "z": round(zf(p), 6),
             "z_lo": round(hh_charge([dict(e, pka=e["pka"] - precision * e["sigma"]) for e in sites], q0, p), 6),
             "z_hi": round(hh_charge([dict(e, pka=e["pka"] + precision * e["sigma"]) for e in sites], q0, p), 6)} for p in at]
    return {
        "ok": True, "engine": "dimorphite", "model": "independent-site Henderson–Hasselbalch sum over Dimorphite-DL rule-table pKa",
        "version": info["version"], "table_sha256": info["table_sha256"], "rdkit_version": info["rdkit_version"],
        "table_order": "file order, except " + "; ".join(f"{m} moved after {a}" for m, a in TABLE_REORDER),
        "precision_sigma": precision,
        "reference_state_smiles": det["reference_state_smiles"], "q0": q0, "fixed_charge_atoms": det["fixed_charge_atoms"],
        "sites": sites, "n_acid": det["n_acid"], "n_base": det["n_base"],
        "duplicate_terms_collapsed": det["duplicate_terms_collapsed"], "truncated_patterns": det["truncated_patterns"],
        "limits": {"z_low_ph": q0 + det["n_base"], "z_high_ph": q0 - det["n_acid"]},
        "curve": curve, "pI": pi["pI"], "pI_note": pi["note"], "z_at": z_at,
        "caveats": [
            "rule-table pKa (moiety mean ± σ from Dimorphite-DL's training set), not a per-molecule prediction; σ is 1–3 units for many moieties — read z_lo/z_hi",
            "independent-site approximation: coupled/adjacent sites (zwitterions, polyphosphates) deviate; run engine=all and read the A/B disagreement (neither engine is authoritative)",
            "every aromatic N without H is counted as a base (pKa 4.35 ± 2.07) and every aromatic N–H as an acid (7.17 ± 2.95) — imidazole/adenine-type rings are over-counted",
            "25 °C, ionic strength 0; tautomers not enumerated",
        ] + ([f"site list TRUNCATED: pattern(s) {det['truncated_patterns']} reached Dimorphite's {MAX_HITS_PER_PATTERN}-hit cap; "
              "the charge is incomplete"] if det["truncated_patterns"] else []),
    }


# ----------------------------------------------------------------------------------------------
# outputs
# ----------------------------------------------------------------------------------------------
def write_csv(path: str, result: dict) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["engine", "ph", "z", "z_lo", "z_hi"])
        for name, blk in result["engines"].items():
            if not blk.get("ok"):
                continue
            for row in blk["curve"]:
                w.writerow([name, row["ph"], row["z"], row.get("z_lo", ""), row.get("z_hi", "")])


def write_svg(path: str, result: dict) -> None:
    """Dependency-free vector plot: Z vs pH, one polyline per engine, envelope band for the rule engine."""
    W, H, L, R, T, B = 720, 440, 70, 20, 40, 60
    blocks = [(n, b) for n, b in result["engines"].items() if b.get("ok")]
    if not blocks:
        return
    xs = [r["ph"] for r in blocks[0][1]["curve"]]
    ys = []
    for _, b in blocks:
        ys += [r["z"] for r in b["curve"]] + [r.get("z_lo", r["z"]) for r in b["curve"]] + [r.get("z_hi", r["z"]) for r in b["curve"]]
    x0, x1 = min(xs), max(xs)
    y0, y1 = math.floor(min(ys) - 0.05), math.ceil(max(ys) + 0.05)
    if y1 == y0:
        y1 = y0 + 1
    sx = lambda x: L + (x - x0) / (x1 - x0) * (W - L - R)  # noqa: E731
    sy = lambda y: T + (y1 - y) / (y1 - y0) * (H - T - B)  # noqa: E731
    colors = {"dimorphite": "#1f77b4"}
    ident = result.get("identity", {})
    title = f"{ident.get('inchikey', '?')}  ·  {ident.get('formula', '')}"
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif" font-size="12">',
             f'<rect width="{W}" height="{H}" fill="white"/>',
             f'<text x="{L}" y="22" font-size="13">{title}</text>']
    # grid + axes
    for gx in range(int(math.ceil(x0)), int(math.floor(x1)) + 1):
        parts.append(f'<line x1="{sx(gx):.1f}" y1="{T}" x2="{sx(gx):.1f}" y2="{H-B}" stroke="#eee"/>')
        parts.append(f'<text x="{sx(gx):.1f}" y="{H-B+16}" text-anchor="middle">{gx}</text>')
    for gy in range(y0, y1 + 1):
        parts.append(f'<line x1="{L}" y1="{sy(gy):.1f}" x2="{W-R}" y2="{sy(gy):.1f}" stroke="{"#999" if gy == 0 else "#eee"}"/>')
        parts.append(f'<text x="{L-6}" y="{sy(gy)+4:.1f}" text-anchor="end">{gy:+d}</text>')
    parts.append(f'<rect x="{L}" y="{T}" width="{W-L-R}" height="{H-T-B}" fill="none" stroke="#333"/>')
    parts.append(f'<text x="{(L+W-R)/2:.1f}" y="{H-14}" text-anchor="middle">pH</text>')
    parts.append(f'<text x="16" y="{(T+H-B)/2:.1f}" text-anchor="middle" transform="rotate(-90 16 {(T+H-B)/2:.1f})">net charge Z</text>')
    legend_y = T + 14
    for name, b in blocks:
        col = colors.get(name, "#333")
        if "z_lo" in b["curve"][0]:
            up = " ".join(f"{sx(r['ph']):.1f},{sy(r['z_hi']):.1f}" for r in b["curve"])
            dn = " ".join(f"{sx(r['ph']):.1f},{sy(r['z_lo']):.1f}" for r in reversed(b["curve"]))
            parts.append(f'<polygon points="{up} {dn}" fill="{col}" fill-opacity="0.12" stroke="none"/>')
        pts = " ".join(f"{sx(r['ph']):.1f},{sy(r['z']):.1f}" for r in b["curve"])
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="2"/>')
        if b.get("pI") is not None:
            parts.append(f'<circle cx="{sx(b["pI"]):.1f}" cy="{sy(0):.1f}" r="4" fill="{col}"/>')
        label = f"{name} v{b.get('version', '?')}" + (f"  pI {b['pI']}" if b.get("pI") is not None else "")
        parts.append(f'<line x1="{W-R-190}" y1="{legend_y}" x2="{W-R-170}" y2="{legend_y}" stroke="{col}" stroke-width="2"/>')
        parts.append(f'<text x="{W-R-164}" y="{legend_y+4}">{label}</text>')
        legend_y += 16
    for p in result.get("at", []):
        parts.append(f'<line x1="{sx(p):.1f}" y1="{T}" x2="{sx(p):.1f}" y2="{H-B}" stroke="#555" stroke-dasharray="4 3"/>')
    parts.append("</svg>")
    with open(path, "w") as fh:
        fh.write("\n".join(parts))


def write_png(path: str, result: dict) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=150)
    for name, b in result["engines"].items():
        if not b.get("ok"):
            continue
        xs = [r["ph"] for r in b["curve"]]
        ax.plot(xs, [r["z"] for r in b["curve"]], label=f"{name} v{b.get('version', '?')}")
        if "z_lo" in b["curve"][0]:
            ax.fill_between(xs, [r["z_lo"] for r in b["curve"]], [r["z_hi"] for r in b["curve"]], alpha=0.12)
        if b.get("pI") is not None:
            ax.plot([b["pI"]], [0], "o")
    for p in result.get("at", []):
        ax.axvline(p, ls="--", lw=0.8, color="0.4")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("pH"); ax.set_ylabel("net charge Z")
    ident = result.get("identity", {})
    ax.set_title(f"{ident.get('inchikey', '?')} · {ident.get('formula', '')}", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    return True


# ----------------------------------------------------------------------------------------------
# main computation
# ----------------------------------------------------------------------------------------------
def _quiet_rdkit() -> None:
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.error")
    except Exception:
        pass


def compute(smiles: str, ph: str = "0:14:0.1", at: list[float] | None = None, engine: str = "auto",
            precision: float = 1.0) -> dict:
    _quiet_rdkit()
    at = list(at or [])
    try:
        lo, hi, step = parse_grid(ph)
    except ValueError as e:
        return {"ok": False, "error": str(e)}
    if engine not in ("auto", "all") + ENGINES:
        return {"ok": False, "error": f"unknown engine {engine!r}; choose auto|all|{'|'.join(ENGINES)}"}
    ident = identity_block(smiles)
    if not ident.get("ok") or not ident.get("valid", True):
        return {"ok": False, "error": "invalid SMILES", "identity": ident, "input": smiles,
                "note": "flag-not-fabricate: no curve for a string RDKit cannot parse"}
    grid = make_grid(lo, hi, step)
    engines: dict[str, dict] = {}
    want = list(ENGINES) if engine in ("auto", "all") else [engine]
    for name in want:
        if name == "dimorphite":
            engines[name] = engine_dimorphite(smiles, grid, at, precision)
    ok_blocks = [b for b in engines.values() if b.get("ok")]
    result = {"ok": bool(ok_blocks), "tool": "molcharge", "tool_version": VERSION, "input": smiles,
              "identity": ident, "grid": {"ph_min": lo, "ph_max": hi, "step": step, "n": len(grid)}, "at": at,
              "engines": engines}
    if len(ok_blocks) >= 2:
        a, b = ok_blocks[0], ok_blocks[1]
        dz = max(abs(x["z"] - y["z"]) for x, y in zip(a["curve"], b["curve"]))
        result["engine_agreement"] = {"engines": [a["engine"], b["engine"]], "max_abs_dz": round(dz, 4),
                                      "pI": [a.get("pI"), b.get("pI")]}
    if not ok_blocks:
        result["error"] = "no engine produced a curve"
    result["note"] = ("Z(pH) = population-weighted mean charge. Each engine block names its method + version + "
                      "rule-table hash: a predicted pKa is a function of (molecule, method). Compare compounds by "
                      "InChIKey; the reference state is the neutralised parent.")
    return result


def write_outputs(result: dict, stem: str) -> dict:
    os.makedirs(os.path.dirname(os.path.abspath(stem)) or ".", exist_ok=True)
    files = {}
    with open(stem + ".json", "w") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=1)
    files["json"] = stem + ".json"
    if result.get("ok"):
        write_csv(stem + ".csv", result); files["csv"] = stem + ".csv"
        write_svg(stem + ".svg", result); files["svg"] = stem + ".svg"
        if write_png(stem + ".png", result):
            files["png"] = stem + ".png"
        else:
            files["png"] = None  # matplotlib absent in this venv — SVG is the vector deliverable
    return files


# ----------------------------------------------------------------------------------------------
# self-test — invariants + golden molecules (identity codes generated with RDKit, never typed)
# ----------------------------------------------------------------------------------------------
GOLDEN = {  # name: (smiles, expected InChIKey)
    "glycine": ("NCC(=O)O", "DHMQDGOQFOQNFH-UHFFFAOYSA-N"),
    "acetic acid": ("CC(=O)O", "QTBSBXVTEAMEQO-UHFFFAOYSA-N"),
    "benzoic acid": ("O=C(O)c1ccccc1", "WPYMKLBDIGXBTP-UHFFFAOYSA-N"),
    "phenol": ("Oc1ccccc1", "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"),
    "aniline": ("Nc1ccccc1", "PAYRUJLWNCNPSJ-UHFFFAOYSA-N"),
    "L-lysine": ("NCCCC[C@H](N)C(=O)O", "KDXKERNSBIXSRK-YFKPBYRVSA-N"),
    "L-glutamic acid": ("N[C@@H](CCC(=O)O)C(=O)O", "WHUUTDBJXJRKMK-VKHMYHEASA-N"),
}
SAM_PUBCHEM = "C[S+](CC[C@H](N)C(=O)[O-])C[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O"   # CID 34755
SAM_KEY = "MEFKEPWMEQBLKI-AIRLBKTGSA-N"


def selftest() -> int:
    import tempfile
    passed, failed = 0, 0

    def gate(desc: str, cond: bool, detail: str = ""):
        nonlocal passed, failed
        if cond:
            passed += 1; print(f"  PASS  {desc}")
        else:
            failed += 1; print(f"  FAIL  {desc}  {detail}")

    print("== molcharge.charge selftest ==")
    r = compute("this is not smiles")
    gate("1 invalid SMILES → ok:false, no curve", r.get("ok") is False and "engines" not in r)
    res = {k: compute(v[0], at=[7.4], engine="dimorphite") for k, v in GOLDEN.items()}
    mono = all(all(c["z"][i + 1] <= c["z"][i] + 1e-9 for i in range(len(c["z"]) - 1))
               for c in ({"z": [row[key] for row in blk["engines"]["dimorphite"]["curve"]]}
                         for blk in res.values() for key in ("z", "z_lo", "z_hi")))
    gate("3 monotone non-increasing Z(pH) for every golden molecule (z, z_lo, z_hi)", mono)
    lim_ok = True
    for k, blk in res.items():
        d = blk["engines"]["dimorphite"]
        z_lo_ph = hh_charge(d["sites"], d["q0"], -60.0); z_hi_ph = hh_charge(d["sites"], d["q0"], 60.0)
        lim_ok &= abs(z_lo_ph - d["limits"]["z_low_ph"]) < 1e-9 and abs(z_hi_ph - d["limits"]["z_high_ph"]) < 1e-9
    gate("4 limits Z(−∞)=q0+N_base, Z(+∞)=q0−N_acid (evaluated analytically)", lim_ok)
    zac = res["acetic acid"]["engines"]["dimorphite"]["z_at"][0]["z"]
    gate("5 acetic acid Z(7.4) in [−1.00, −0.99]", -1.0 <= zac <= -0.99, f"got {zac}")
    pig = res["glycine"]["engines"]["dimorphite"]["pI"]
    gate("6 glycine pI within 0.3 of 5.97 (textbook)", pig is not None and abs(pig - 5.97) <= 0.3, f"got {pig}")
    sam = compute(SAM_PUBCHEM, engine="dimorphite")
    d = sam["engines"]["dimorphite"]
    tma = compute("C[N+](C)(C)C", engine="dimorphite")["engines"]["dimorphite"]
    nb = compute("O=[N+]([O-])c1ccccc1", engine="dimorphite")["engines"]["dimorphite"]
    gate("7 permanent charge: SAM q0=+1 (sulfonium), Me4N+ Z≡+1, nitrobenzene Z≡0",
         d["q0"] == 1 and any(a["element"] == "S" for a in d["fixed_charge_atoms"])
         and all(abs(r["z"] - 1) < 1e-9 for r in tma["curve"]) and all(abs(r["z"]) < 1e-9 for r in nb["curve"]),
         f"q0={d['q0']} tma={tma['curve'][0]['z']} nb={nb['curve'][70]['z']}")
    a = compute("OC(=O)c1ccccc1", engine="dimorphite"); b = compute("c1ccccc1C(O)=O", engine="dimorphite")
    gate("8 two spellings of benzoic acid → same InChIKey, identical curve",
         a["identity"]["inchikey"] == b["identity"]["inchikey"] == GOLDEN["benzoic acid"][1]
         and a["engines"]["dimorphite"]["curve"] == b["engines"]["dimorphite"]["curve"])
    gate("9 identity block: golden InChIKeys reproduce; SAM flags 1 undefined stereocentre (sulfonium)",
         all(res[k]["identity"].get("inchikey") == v[1] for k, v in GOLDEN.items())
         and sam["identity"].get("inchikey") == SAM_KEY and sam["identity"].get("n_stereocentres_undefined") == 1)
    gate("10 method label on the engine block (engine/version/table_sha256/rdkit_version)",
         all(d.get(k) for k in ("engine", "version", "table_sha256", "rdkit_version")))
    gate("11 bad grid rejected (step≤0, hi≤lo)", compute("CCO", ph="0:14:0").get("ok") is False and compute("CCO", ph="7:7").get("ok") is False)
    env_ok = all(r["z_lo"] - 1e-9 <= r["z"] <= r["z_hi"] + 1e-9 for blk in res.values() for r in blk["engines"]["dimorphite"]["curve"])
    gate("12 envelope order z_lo ≤ z ≤ z_hi everywhere", env_ok)
    with tempfile.TemporaryDirectory() as td:
        files = write_outputs(res["glycine"], os.path.join(td, "gly"))
        back = json.load(open(files["json"]))
        n_csv = sum(1 for _ in open(files["csv"])) - 1
        svg_ok = open(files["svg"]).read().startswith("<svg")
        gate("13 outputs: JSON round-trips, CSV rows = grid points, SVG written",
             back["identity"]["inchikey"] == GOLDEN["glycine"][1] and n_csv == res["glycine"]["grid"]["n"] and svg_ok,
             f"csv rows {n_csv} vs grid {res['glycine']['grid']['n']}")
    zat = res["glycine"]["engines"]["dimorphite"]["z_at"][0]
    gly = res["glycine"]["engines"]["dimorphite"]
    gate("14 --at point equals the function value (to the 6-dp rounding of the output)",
         abs(zat["z"] - hh_charge(gly["sites"], gly["q0"], 7.4)) < 5e-7)
    ph_site = res["phenol"]["engines"]["dimorphite"]["sites"]
    gly_kinds = sorted(e["kind"] for e in gly["sites"])
    gate("15 site table read live from Dimorphite (phenol pKa 7.065 ± 3.277; glycine = 1 acid + 1 base)",
         len(ph_site) == 1 and abs(ph_site[0]["pka"] - 7.065) < 1e-3 and abs(ph_site[0]["sigma"] - 3.277) < 1e-3
         and gly_kinds == ["acid", "base"], f"{ph_site} {gly_kinds}")
    # 17 — site COUNTING. Gates 3/4/12 cannot catch a double-counted site: a duplicated HH term is still monotone
    # and still meets its own limits. Before the 2026-10-05 fix triethylamine read Z@7.4 +2.56, H3PO4 six terms.
    amines = {s: compute(s, at=[7.4], engine="dimorphite")["engines"]["dimorphite"] for s in ("CCN", "CCNCC", "CCN(CC)CC")}
    z_am = [b["z_at"][0]["z"] for b in amines.values()]
    acet, h3po4 = detect_sites_dimorphite("CC(N)=O"), detect_sites_dimorphite("OP(=O)(O)O")

    def _one_site_per_site(det):
        keys = [(e["atom_idx"], e["kind"], e["pka"]) for e in det["sites"]]
        per = {}
        for e in det["sites"]:
            per.setdefault((e["atom_idx"], e["kind"]), set()).add(e["pka"])
        return len(keys) == len(set(keys)) and all(len(v) == 1 for v in per.values())
    probe = [detect_sites_dimorphite(v[0]) for v in GOLDEN.values()] + [acet, h3po4] + [
        detect_sites_dimorphite(s) for s in ("CCNCC", "CCN(CC)CC", "C1CCNCC1", "CS(N)(=O)=O", "CN1C=CCC(C(N)=O)=C1")]
    gate("17 one entry per site: 1°/2°/3° amine = 1 base each with equal Z(7.4); NH2 amide = 1 acid; H3PO4 = 2 "
         "ionizations; no repeated (atom,kind,pKa) and no atom with two pKa of one kind",
         all(b["n_base"] == 1 for b in amines.values()) and max(z_am) - min(z_am) < 1e-6
         and acet["n_acid"] == 1 and len(h3po4["sites"]) == 2 and all(_one_site_per_site(d) for d in probe),
         f"z={z_am} amide_acid={acet['n_acid']} h3po4={len(h3po4['sites'])}")
    # 18 — the other direction: dedup must never MERGE distinct sites. A dedup keyed too coarsely (pattern name
    # alone, or (kind, pKa) without the atom) passes gate 17 while cutting lysine to one base and pyrophosphoric
    # acid to two acids (mutation review 2026-10-05). (n_acid, n_base) expected per molecule:
    lower = {"NCCN": (0, 2), "OC(=O)CC(=O)O": (2, 0), "NC(=O)CC(N)=O": (2, 0), "NC(N)=O": (2, 0),
             "OP(=O)(O)OP(=O)(O)O": (4, 0), "OC(=O)CN(CCN(CC(=O)O)CC(=O)O)CC(=O)O": (4, 2),
             "O=C1NN1": (2, 0)}       # diaziridinone: two distinct N–H on one heavy-atom set (small-ring case)
    got = {s: (dd["n_acid"], dd["n_base"]) for s, dd in ((s, detect_sites_dimorphite(s)) for s in lower)}
    lys, glu = (res[k]["engines"]["dimorphite"] for k in ("L-lysine", "L-glutamic acid"))
    gate("18 distinct sites stay distinct: lysine 1+2, glutamate 2+1, ethylenediamine 2 base, malonic acid / "
         "malonamide / urea 2 acid, pyrophosphoric acid 4, EDTA 4+2, diaziridinone 2; terms collapsed Et3N 2, H3PO4 4",
         (lys["n_acid"], lys["n_base"]) == (1, 2) and (glu["n_acid"], glu["n_base"]) == (2, 1) and got == lower
         and detect_sites_dimorphite("CCN(CC)CC")["duplicate_terms_collapsed"] == 2
         and h3po4["duplicate_terms_collapsed"] == 4,
         f"{got} lys={(lys['n_acid'], lys['n_base'])} glu={(glu['n_acid'], glu['n_base'])}")
    # 19 no silent truncation: Dimorphite's default 50-hit cap cut a 41-N oligo-ethyleneimine to 25 base sites
    pei = detect_sites_dimorphite("N" + "CCN" * 40)
    gate("19 no silent site cap: 41-N oligo-ethyleneimine → 41 base sites, nothing truncated",
         pei["n_base"] == 41 and not pei["truncated_patterns"], f"n_base={pei['n_base']} truncated={pei['truncated_patterns']}")
    # 22 table reorder: a phosphoramidate monoester P–OH is no longer lost to the amide/amine patterns
    agr = detect_sites_dimorphite("CC(C)[C@@H](O)[C@H](O)C(=O)NP(=O)(O)OC[C@@H]1C[C@H](O)[C@H](n2cnc3c(NP(=O)(O)O"
                                  "[C@@H]4O[C@@H]([C@H](O)CO)[C@@H](O)[C@H]4O)ncnc32)O1")   # agrocin 84 (published)
    mini = detect_sites_dimorphite("CNP(=O)(O)OC")
    n_pe = sum(1 for e in agr["sites"] if e["moiety"] == "Phosphonate_ester")
    z_agr, z_mini = hh_charge(agr["sites"], agr["q0"], 7.4), hh_charge(mini["sites"], mini["q0"], 7.4)
    gate("22 table reorder: agrocin 84 keeps both phosphoramidate P–OH (Z@7.4 ≈ −2); CNP(=O)(O)OC is an acid, not an amine",
         n_pe == 2 and abs(z_agr + 2.0) < 0.05 and abs(z_mini + 1.0) < 0.01 and mini["n_base"] == 0,
         f"Phosphonate_ester sites {n_pe}, Z agrocin {z_agr:+.3f}, Z mini {z_mini:+.3f}, mini bases {mini['n_base']}")
    print(f"selftest: {passed}/{passed + failed} passed")
    return 0 if failed == 0 else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="net charge Z vs pH from a SMILES (local; Dimorphite-DL rule table, Henderson–Hasselbalch sum)")
    ap.add_argument("--smiles")
    ap.add_argument("--ph", default="0:14:0.1", help="grid lo:hi[:step] (default 0:14:0.1)")
    ap.add_argument("--at", default="", help="comma-separated pH values to report explicitly, e.g. 7.4,8.0")
    ap.add_argument("--engine", default="auto", choices=("auto", "all") + ENGINES)
    ap.add_argument("--precision", type=float, default=1.0, help="σ multiplier for the rule-table envelope (default 1.0)")
    ap.add_argument("--out", help="output STEM → STEM.json/.csv/.svg (+ .png when matplotlib is importable)")
    ap.add_argument("--json", action="store_true", help="print the full JSON result to stdout")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if not a.smiles:
        ap.error("--smiles is required (or --selftest)")
    try:
        at = [float(x) for x in a.at.split(",") if x.strip()]
    except ValueError:
        print(json.dumps({"ok": False, "error": f"--at must be comma-separated numbers, got {a.at!r}"})); return 2
    result = compute(a.smiles, ph=a.ph, at=at, engine=a.engine, precision=a.precision)
    if a.out:
        result["files"] = write_outputs(result, a.out)
    if a.json:
        print(json.dumps(result, ensure_ascii=False))
    else:
        ident = result.get("identity", {})
        print(f"input: {a.smiles}")
        print(f"identity: {ident.get('canonical_smiles')}  InChIKey {ident.get('inchikey')}  {ident.get('formula')}  MW {ident.get('mol_weight')}")
        if ident.get("stereo_warning"):
            print(f"  ⚠ {ident['stereo_warning']}")
        if not result.get("ok"):
            print(f"ERROR: {result.get('error')}"); return 1
        for name, blk in result["engines"].items():
            if not blk.get("ok"):
                print(f"[{name}] not run: {blk.get('error')}  {blk.get('hint', '')}"); continue
            print(f"[{name}] {blk['model']}  (v{blk.get('version')}; hash {str(blk.get('table_sha256'))[:12]}…)")
            if name == "dimorphite":
                print(f"  reference state {blk['reference_state_smiles']}   q0={blk['q0']}   sites: {blk['n_acid']} acid / {blk['n_base']} base")
                for e in blk["sites"]:
                    print(f"    atom {e['atom_idx']:>3} {e['element']:<2} {e['moiety']:<34} {e['kind']:<4} pKa {e['pka']:6.2f} ± {e['sigma']:.2f}")
            print(f"  pI: {blk.get('pI')}  ({blk.get('pI_note')})")
            for p in blk.get("z_at", []):
                band = f"  [{p['z_lo']:+.3f}, {p['z_hi']:+.3f}]" if "z_lo" in p else ""
                print(f"  Z(pH {p['ph']}) = {p['z']:+.3f}{band}")
            step = max(1, len(blk["curve"]) // 14)
            print("  curve: " + "  ".join(f"{r['ph']:g}:{r['z']:+.2f}" for r in blk["curve"][::step]))
        if "engine_agreement" in result:
            print(f"engine agreement: {result['engine_agreement']}")
        if result.get("files"):
            print("files: " + ", ".join(f"{k}={v}" for k, v in result["files"].items()))
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
