#!/usr/bin/env python3
"""molcharge.properties — small-molecule physicochemical properties from a SMILES. Local, offline, pure RDKit.

The identity + descriptor half of molcharge.
Its sibling `molcharge.charge` owns the pH-dependent charge / pKa half; this file
optionally folds in that module's pI and Z@pH (--with-charge) so one call gives the whole property block.

Every value names its METHOD (a computed descriptor is a function of (molecule, method), never a measured
property): logP/MR = Crippen; TPSA = Ertl 2000; H-bond counts = RDKit Lipinski definitions; QED = Bickerton
2012. Rule-of-5 / Veber are reported as VIOLATION COUNTS, never a pass/fail verdict. Invalid SMILES → ok:false.

Usage:
  molcharge-props --smiles "CC(=O)Oc1ccccc1C(=O)O" [--with-charge] [--ph 7.4] [--json] [--out STEM]
  molcharge-props --selftest
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

VERSION = "1.0.0"


def _quiet_rdkit() -> None:
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass


def properties(smiles: str, *, with_charge: bool = False, ph: float = 7.4) -> dict:
    _quiet_rdkit()
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, inchi
    except Exception as e:  # pragma: no cover
        return {"ok": False, "error": f"rdkit unavailable: {e}", "note": "pip install rdkit"}
    m = Chem.MolFromSmiles(smiles or "")
    if m is None:
        return {"ok": False, "valid": False, "input": smiles,
                "note": "invalid SMILES — flagged, not coerced (flag-not-fabricate)"}

    ident = {"canonical_smiles": Chem.MolToSmiles(m), "formula": rdMolDescriptors.CalcMolFormula(m),
             "mol_weight": round(Descriptors.MolWt(m), 2), "exact_mass": round(Descriptors.ExactMolWt(m), 4),
             "n_heavy_atoms": m.GetNumHeavyAtoms(), "formal_charge": Chem.GetFormalCharge(m)}
    try:
        key = inchi.MolToInchiKey(m) if inchi.INCHI_AVAILABLE else None
        ident["inchikey"] = key
        ident["inchikey_skeleton"] = key.split("-")[0] if key else None
    except Exception:
        ident["inchikey"] = None
    try:
        unspec = [i for i, t in Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False) if t == "?"]
        ident["n_stereocentres_undefined"] = len(unspec)
        if unspec:
            ident["stereo_warning"] = (f"{len(unspec)} undefined stereocentre(s) at atom idx {unspec} — this string does "
                                       "NOT identify a single stereoisomer; pin them before using it as an identity")
    except Exception:
        pass

    logp = round(Crippen.MolLogP(m), 3)
    tpsa = round(rdMolDescriptors.CalcTPSA(m), 2)
    hbd = Lipinski.NumHDonors(m); hba = Lipinski.NumHAcceptors(m)
    rotb = Lipinski.NumRotatableBonds(m)
    desc = {
        "logP_crippen": logp, "molar_refractivity_crippen": round(Crippen.MolMR(m), 3),
        "tpsa_ertl": tpsa, "h_bond_donors": hbd, "h_bond_acceptors": hba,
        "rotatable_bonds": rotb, "rings": rdMolDescriptors.CalcNumRings(m),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(m),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(m), 3),
        "heavy_atoms": m.GetNumHeavyAtoms(),
    }
    try:
        from rdkit.Chem import QED
        desc["qed_bickerton"] = round(QED.qed(m), 3)
    except Exception:
        desc["qed_bickerton"] = None

    ro5 = {"MW>500": ident["mol_weight"] > 500, "logP>5": logp > 5, "HBD>5": hbd > 5, "HBA>10": hba > 10}
    veber = {"rotatable_bonds>10": rotb > 10, "TPSA>140": tpsa > 140}
    druglike = {
        "lipinski_ro5_violations": sum(ro5.values()), "lipinski_ro5_detail": {k: v for k, v in ro5.items() if v} or None,
        "veber_violations": sum(veber.values()), "veber_detail": {k: v for k, v in veber.items() if v} or None,
        "note": "violation COUNTS only — Ro5/Veber are heuristics for oral drug-likeness, NOT a pass/fail gate; "
                "natural-product substrates routinely and legitimately violate them",
    }

    out = {"ok": True, "valid": True, "tool": "mol_properties", "tool_version": VERSION, "input": smiles,
           "identity": ident, "descriptors": desc, "druglikeness": druglike,
           "methods": {"logP": "Crippen (Wildman-Crippen 1999)", "MR": "Crippen", "TPSA": "Ertl 2000",
                       "HBD/HBA/rotatable": "RDKit Lipinski definitions", "QED": "Bickerton 2012",
                       "formal_charge": "from the SMILES as written (NOT the charge at any pH — see charge block)"},
           "note": "computed descriptors (each a function of molecule+method). For the actual charge in solution use "
                   "the molcharge charge block — formal_charge above is only the drawn state."}

    if with_charge:
        try:
            from . import charge as cvp
            cr = cvp.compute(smiles, ph="0:14:0.1", at=[ph], engine="dimorphite")
            blk = (cr.get("engines") or {}).get("dimorphite") or {}
            if blk.get("ok"):
                out["charge"] = {"engine": "dimorphite (rule table, ±σ)", "version": blk.get("version"),
                                 "pI": blk.get("pI"), "q0_permanent": blk.get("q0"),
                                 "Z_at_ph": {str(ph): (blk["z_at"][0]["z"] if blk.get("z_at") else None)},
                                 "n_acid_sites": blk.get("n_acid"), "n_base_sites": blk.get("n_base"),
                                 "note": "rule-table engine — coarse; for phosphoramidates, nucleotides and unusual scaffolds "
                                         "check the ionisable-group list by hand (see the README limitations)"}
            else:
                out["charge"] = {"ok": False, "error": blk.get("error", "charge engine unavailable")}
        except Exception as e:
            out["charge"] = {"ok": False, "error": f"charge_vs_ph unavailable: {e}"}
    return out


def _fmt(r: dict) -> str:
    if not r.get("ok"):
        return f"ERROR: {r.get('error') or 'invalid SMILES'}"
    i, d, dl = r["identity"], r["descriptors"], r["druglikeness"]
    L = [f"identity : {i['canonical_smiles']}",
         f"           InChIKey {i.get('inchikey')}  {i['formula']}  MW {i['mol_weight']}  exact {i['exact_mass']}  q(drawn) {i['formal_charge']}"]
    if i.get("stereo_warning"):
        L.append(f"  ⚠ {i['stereo_warning']}")
    L.append(f"lipophilicity : logP {d['logP_crippen']} (Crippen)   MR {d['molar_refractivity_crippen']}")
    L.append(f"polarity/H-bond : TPSA {d['tpsa_ertl']} Å²   HBD {d['h_bond_donors']}   HBA {d['h_bond_acceptors']}")
    L.append(f"shape/flex : rotatable {d['rotatable_bonds']}   rings {d['rings']} (aromatic {d['aromatic_rings']})   fsp3 {d['fraction_csp3']}   QED {d['qed_bickerton']}")
    L.append(f"drug-likeness : Ro5 violations {dl['lipinski_ro5_violations']} {dl['lipinski_ro5_detail'] or ''}   Veber violations {dl['veber_violations']} {dl['veber_detail'] or ''}")
    if "charge" in r:
        c = r["charge"]
        if "Z_at_ph" in c:   # pI is legitimately absent when Z never crosses 0 on pH 0-14 (a plain acid or base)
            pi = c["pI"] if c.get("pI") is not None else "none on pH 0-14"
            L.append(f"charge@pH : pI {pi}   Z(pH {list(c['Z_at_ph'])[0]}) {list(c['Z_at_ph'].values())[0]}   (engine {c['engine']}; {c['n_acid_sites']} acid/{c['n_base_sites']} base sites)")
        else:
            L.append(f"charge@pH : not computed — {c.get('error')}")
    return "\n".join(L)


def selftest() -> int:
    _quiet_rdkit()
    passed = failed = 0

    def gate(desc, cond, extra=""):
        nonlocal passed, failed
        if cond:
            passed += 1; print(f"  PASS  {desc}")
        else:
            failed += 1; print(f"  FAIL  {desc}  {extra}")

    print("== mol_properties selftest ==")
    r = properties("not a smiles at all")
    gate("1 invalid SMILES → ok:false, no descriptors", r.get("ok") is False and "descriptors" not in r)

    asp = properties("CC(=O)Oc1ccccc1C(=O)O")   # aspirin
    gate("2 aspirin: formula C9H8O4, MW ~180.16, InChIKey BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
         asp["identity"]["formula"] == "C9H8O4" and abs(asp["identity"]["mol_weight"] - 180.16) < 0.1
         and asp["identity"]["inchikey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    gate("3 aspirin descriptors in known ranges (logP 1.0-1.5, TPSA ~63.6, HBD 1, HBA 3-4)",
         1.0 <= asp["descriptors"]["logP_crippen"] <= 1.5 and abs(asp["descriptors"]["tpsa_ertl"] - 63.6) < 0.5
         and asp["descriptors"]["h_bond_donors"] == 1 and asp["descriptors"]["h_bond_acceptors"] in (3, 4),
         str(asp["descriptors"]))

    caf = properties("Cn1cnc2c1c(=O)n(C)c(=O)n2C")   # caffeine
    gate("4 caffeine: MW ~194.19, HBD 0, Ro5 violations 0",
         abs(caf["identity"]["mol_weight"] - 194.19) < 0.1 and caf["descriptors"]["h_bond_donors"] == 0
         and caf["druglikeness"]["lipinski_ro5_violations"] == 0, str(caf["identity"]["mol_weight"]))

    benz = properties("c1ccccc1")   # benzene: TPSA 0, logP ~1.7-2.0
    gate("5 benzene: TPSA 0, logP 1.5-2.2, fsp3 0", benz["descriptors"]["tpsa_ertl"] == 0.0
         and 1.5 <= benz["descriptors"]["logP_crippen"] <= 2.2 and benz["descriptors"]["fraction_csp3"] == 0.0)

    big = properties("CC(C)[C@H]([C@@H](C(=O)NP(=O)(O)OC[C@@H]1C[C@@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)NP(=O)(O)O[C@H]4[C@@H]([C@@H]([C@@H](O4)[C@@H](CO)O)O)O)O)O)O")  # agrocin 84
    gate("6 agrocin 84: MW>700 → multiple Ro5 violations; 2 undefined stereocentres flagged",
         big["identity"]["mol_weight"] > 700 and big["druglikeness"]["lipinski_ro5_violations"] >= 2
         and big["identity"]["n_stereocentres_undefined"] == 2, str(big["identity"]["mol_weight"]))

    wc = properties("NCC(=O)O", with_charge=True, ph=7.4)   # glycine + charge fold-in
    gate("7 --with-charge folds in the charge block (glycine pI ~5.8, Z@7.4 near 0)",
         "charge" in wc and wc["charge"].get("pI") is not None and abs(wc["charge"]["pI"] - 5.8) < 0.5,
         str(wc.get("charge")))

    ac = properties("CC(=O)O", with_charge=True, ph=7.4)   # acetic acid: no pI on 0-14, Z@7.4 ~ -1
    gate("10 charge line is printed even without a pI (acetic acid: Z@7.4 shown, pI 'none')",
         "charge@pH" in _fmt(ac) and "none on pH 0-14" in _fmt(ac) and ac["charge"].get("pI") is None, _fmt(ac).splitlines()[-1])

    gate("8 every descriptor names its method", set(asp["methods"]) >= {"logP", "TPSA", "QED"}
         and "Crippen" in asp["methods"]["logP"])

    two = (properties("OC(=O)c1ccccc1")["identity"]["inchikey"], properties("c1ccccc1C(O)=O")["identity"]["inchikey"])
    gate("9 two SMILES spellings of benzoic acid → same InChIKey", two[0] == two[1] and two[0] == "WPYMKLBDIGXBTP-UHFFFAOYSA-N")
    print(f"selftest: {passed}/{passed + failed} passed")
    return 0 if failed == 0 else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="small-molecule physicochemical properties from a SMILES (local, RDKit)")
    ap.add_argument("--smiles")
    ap.add_argument("--with-charge", action="store_true", help="fold in pI and Z@pH from the molcharge.charge module")
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--out", help="write STEM.json")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if not a.smiles:
        ap.error("--smiles is required (or --selftest)")
    r = properties(a.smiles, with_charge=a.with_charge, ph=a.ph)
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
        json.dump(r, open(a.out + ".json", "w"), ensure_ascii=False, indent=1)
    if a.json:
        print(json.dumps(r, ensure_ascii=False))
    else:
        print(_fmt(r))
    return 0 if r.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
