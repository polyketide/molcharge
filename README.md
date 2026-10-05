# molcharge

Net charge versus pH, ionisable sites with pKa, and the isoelectric point of a small molecule, from a SMILES string.
Runs locally and offline; the SMILES never leaves the machine.

The engine is the [Dimorphite-DL](https://github.com/durrantlab/dimorphite_dl) 2.0.2 rule table (41 ionisable
moieties, each with a mean pKa and its standard deviation) combined as an independent-site Henderson–Hasselbalch
sum. Every result names the engine, its version and the SHA-256 of the rule table, because a predicted pKa is a
function of the molecule *and* the method.

## Install

```bash
pip install git+https://github.com/polyketide/molcharge
```

Requires Python ≥ 3.10. Dependencies: `dimorphite-dl==2.0.2` (pinned, see "Rule-table order" below) and RDKit.
`matplotlib` is optional (PNG output; an SVG is always written).

## Use

```bash
molcharge --smiles "NCC(=O)O" --at 7.4            # glycine: site table, pI, Z at pH 7.4, coarse curve
molcharge --smiles "NCC(=O)O" --out glycine       # writes glycine.json, .csv, .svg (+ .png)
molcharge-props --smiles "CC(=O)Oc1ccccc1C(=O)O" --with-charge   # RDKit descriptors + pI and Z
molcharge --selftest                               # 18 gates, about one second
```

```python
from molcharge import compute
r = compute("NCC(=O)O", at=[7.4])
blk = r["engines"]["dimorphite"]
blk["pI"], blk["z_at"][0]          # isoelectric point; Z with its ±σ band at pH 7.4
blk["sites"]                       # one entry per ionisable site: atom, moiety, acid/base, pKa ± σ
```

## Method

    Z(pH) = q0 + Σ_base 1/(1 + 10^(pH − pKa)) − Σ_acid 1/(1 + 10^(pKa − pH))

`q0` is the permanent formal charge (sulfonium, quaternary ammonium). The `z_lo`/`z_hi` band shifts every pKa by
±σ of its moiety (`--precision` scales it). The curve is monotone non-increasing, so the pI, when one exists in the
range, is unique and is found by bisection.

Three corrections are applied over calling Dimorphite-DL directly. Each is guarded by a self-test gate.

1. **One site per site.** Dimorphite-DL returns one record per SMARTS hit. A secondary or tertiary amine N is hit
   once per C–N bond, so summing the records counts it two or three times: triethylamine reads Z(7.4) = +2.56
   instead of +0.85, and phosphoric acid gets six terms for two ionisations. Records are merged per site.
2. **No silent truncation.** Dimorphite-DL stops at 50 hits per pattern by default. Because of (1) that cut any
   molecule with more than about 25 aliphatic N to 25 base sites. The cap is raised to 1000, and reaching it is
   reported in the output.
3. **Rule-table order.** Dimorphite-DL matches patterns in file order and protects every matched atom. The pattern
   for a phosphoramidate or phosphonate monoester, P(=O)(OR)(X)–OH, comes after the amide, aniline and amine
   patterns. When the N on P also carries a carbon, those patterns take the N first and the P–OH is lost; the N is
   then read as an amine or a weak acid. Here that pattern is matched immediately after its sibling `Phosphonate`.
   Agrocin 84 (InChIKey `FIMRCGIHIAIVOL-LQIBLKNOSA-N`) moves from Z(7.4) = −0.995 to −1.996, in line with its two
   phosphoramidate P–O⁻ groups. This uses a private method of Dimorphite-DL 2.0.2, which is why the version is
   pinned; a missing pattern raises an error rather than falling back to the original order.

## Validation

Measured with `molcharge-validate` against experimental data that is **not** redistributed here (fetch it yourself,
point `MOLCHARGE_DATA_DIR` at it, and check each dataset's licence): the IUPAC Digitized pKa Dataset (Zheng & Lafontant-Joseph),
high-confidence file, release v2.3e (doi:10.5281/zenodo.21533589; all releases doi:10.5281/zenodo.7236452; CC BY-NC 4.0,
reproduced by permission of IUPAC; `molcharge-validate` reports which release it read, by MD5), and the SAMPL6 pKa challenge measurements. Nucleotide literature
values ship in `reference_literature.json` with a note on what was checked: its one DOI resolves and agrees with
Crossref, PubMed and the publisher page on title, journal, volume, pages and year; the Alberty pK is traced to
Alberty & Goldberg, Biochemistry 31, 10610 (1992), at I = 0.25 M, and the
ATP values at I = 0.1 M to Sigel et al., Inorg. Chem. 26, 2149 (1987); no value was re-checked against the full texts.

| reference set | n | result |
|---|---|---|
| IUPAC, curve level (median 20–30 °C macro pKa → reference Z(pH)) | 23 molecules | median \|ΔZ(7.4)\| 0.133; median max \|ΔZ\| 0.755; median \|ΔpI\| 0.7 |
| IUPAC + literature, pKa read-off | 43 pKa | RMSE 1.627 |
| SAMPL6, pKa read-off | 31 pKa (24 molecules) | RMSE 1.274, MAE 0.994 |
| ATP / ADP / AMP literature pKa | 7 pKa | RMSE 0.177; Z(7.4) −3.88 / −2.88 / −1.88 |
| secondary/tertiary amines, amides, sulfonamides, imides (IUPAC, reported separately) | 19 molecules | median \|ΔZ(7.4)\| 0.148; pKa RMSE 2.105 |

The last row is a coverage panel added after correction (1) left every headline number unchanged: the headline sets
contain no secondary or tertiary amines. On that panel the median |ΔZ(7.4)| is 0.704 without the correction and
0.148 with it. The pKa read-off RMSE moves the other way, from 1.764 to 2.105. The duplicated terms happened to put
a crossing nearer the true pKa, so read-off RMSE cannot detect double counting; |ΔZ| can.

## Limitations

- **Rule-table pKa, not per-molecule predictions.** σ is 1–3 pKa units for many moieties; read the band.
- **Aliphatic amines.** All aliphatic amines (primary, secondary, tertiary) get the same pKa, 8.16. Measured
  values are mostly 10–11, so Z at physiological pH reads about 0.15 low per amine.
- **Sulfonamides** are assigned pKa 7.92, about two units too acidic. **Phenols** are assigned 7.07 ± 3.28, against
  about 10 for phenol itself.
- **Aromatic N** is counted generously. Every aromatic N without H is a base (pKa 4.35) and every aromatic N–H an
  acid, so adenine-type rings over-count below pH 5.
- **Reduced nicotinamide** (1,4-dihydropyridine N-glycoside, as in NADH) is read as an aliphatic amine, about
  +0.85 at pH 7.4, although that N does not protonate.
- **Not modelled:** coupled sites (zwitterions, polyphosphates), tautomers, temperature and ionic strength (25 °C,
  I = 0 assumed), and N–H deprotonation of phosphoramidates.
- For anything that matters, list the ionisable groups by hand and treat the output as a starting point. A
  computed pKa does not replace a measured one.

## Citation

If you use this, cite Dimorphite-DL: Ropp, P. J. et al. *J. Cheminform.* **11**, 14 (2019),
doi:10.1186/s13321-019-0336-9.

## Licence

Apache-2.0. See `LICENSE` and `NOTICE`.
