# Solvent-accessible surface area

Ensemble accessibility of functionally annotated residues in the restrained
MPXV A29 model, used to test the model against published mutagenesis that was
never supplied to the prediction.

**Date:** 24–25 September 2026

---

## Why the ensemble rather than a single model

Accessibility was first measured on the selected model alone
(`absolute_area_A29_s2m0_pymol.txt`, PyMOL `get_area`). That measurement is
retained for reference but is **not** the reported result.

Single-model accessibility is unreliable in the C-terminal half of this
protein. The ensemble RMSF rises from ~1 Å in the coiled-coil core to 3–4 Å by
residue 100 and 12–16 Å in the N-terminal region, so a chain whose terminus
happens to fold across the bundle in one model will bury residues that are
exposed in another. The preliminary measurement showed exactly this: Thr88 read
91 Å² in chain A but 3 and 5 Å² in chains B and C — an order-of-magnitude
difference between chemically identical chains.

The ensemble analysis therefore reports the distribution across all accepted
models instead.

---

## Method

| | |
| --- | --- |
| Script | `ensemble_sasa.py` |
| SASA algorithm | Shrake–Rupley, via Biopython `Bio.PDB.SASA` |
| Relative values | measured area ÷ theoretical maximum, Tien et al. (2013) *PLoS ONE* 8:e80635 |
| Buried | < 20 % relative accessibility |
| Exposed | > 25 % relative accessibility |
| Models | 100 generated (Boltz-2 v2.2.1, restraint set S1, seeds 1–20) |
| Accepted | 97 — antiparallel topology and ipTM ≥ 0.60 |

Excluded: seed 17 sample 4 (all-parallel topology, restraints violated by
26 Å); seed 14 sample 4 and seed 16 sample 4 (ipTM 0.595 and 0.564).

Every accepted model is relabelled before measurement so that **C is the
antiparallel chain, A presents the heptad *g* face to C, and B presents the
*e* face**. Boltz assigns chain identifiers arbitrarily, so without this the
per-chain columns would not be comparable across models.

---

## The two runs

### `full_length/` — residues 1–110

```bash
python3 ensemble_sasa.py <results_dir> --out .
```

The model as built. **These are the reportable values for the hydrophobic
core, the cysteine pair and Asn75.**

### `truncated_44-110/` — residues 44–110

```bash
python3 ensemble_sasa.py <results_dir> --keep 44 110 --out .
```

An occlusion control. Residues 1–43 are removed before measurement to test
whether the flexible N-terminal region covers the A18-binding surface.
**These are the reportable values for the A18-binding surface (residues
81–99).**

> **Do not read core values from the truncated run.** Truncation creates an
> artificial chain terminus at residue 44, so residues within roughly two
> helical turns of the cut read as more exposed than they are. Leu47 rises
> from 3.8 % to 14.1 %, Asn75 from 13.6 % to 26.7 % and Cys72 from 2.7 % to
> 25.7 % — artefacts of the cut, not structural features.

---

## Results

### Hydrophobic core — *a* and *d* positions (full-length run)

| Residue | Heptad | Mean relative SASA | Verdict |
| --- | --- | --- | --- |
| Leu47 | d | 3.8 % | buried |
| Leu51 | a | 3.7 % | buried |
| Leu54 | d | 3.4 % | buried |
| Ile58 | a | 2.9 % | buried |
| Ile61 | d | 2.5 % | buried |
| Phe65 | a | 5.4 % | buried |

All six are buried in all three chains across all 97 models, with maxima never
exceeding 17 %. Alanine substitution at these positions abolishes self-assembly
in VACV A27; the model places every one of them in the core without any such
information having been supplied.

### A18-binding surface (truncated run)

| Residue | Chain A | Chain B | Chain C | Verdict |
| --- | --- | --- | --- | --- |
| Arg81 | 53.7 % | 48.3 % | 53.2 % | exposed |
| Glu87 | 37.4 % | 41.8 % | 41.9 % | exposed |
| Thr88 | 49.7 % | 48.7 % | 49.9 % | exposed |
| Ile94 | 40.9 % | 41.6 % | 41.3 % | exposed |
| Ser95 | 46.0 % | 45.6 % | 45.7 % | exposed |
| Lys98 | 41.1 % | 40.6 % | 40.7 % | exposed |
| Lys99 | 46.9 % | 34.5 % | 47.3 % | exposed |

**The control resolves the apparent asymmetry.** In the full-length run these
residues appeared systematically less accessible in chain C than chain A
(Thr88: 48/31/20 %; Ser95: 44/30/20 %). With the N-terminal region removed the
values converge (Thr88: 50/49/50 %; Ser95: 46/46/46 %) and the standard
deviations collapse — for Thr88 from 21 to 2 percentage points. The asymmetry
was caused by the flexible N-terminal tails of the partner chains lying across
the surface, not by any structural difference. All three chains present an
equivalent binding face.

### Cysteine pair (full-length run)

| Residue | Heptad | Mean relative SASA | Max | Verdict |
| --- | --- | --- | --- | --- |
| Cys71 | g | 15.0 % | 39.5 % | partially accessible |
| Cys72 | a | 2.7 % | 26.1 % | buried |

Consistent with the register: Cys71 at a *g* position faces the surface and is
accessible enough to form the intermolecular disulfide to A28, while Cys72 at
an *a* position is core-facing, forming the anomalous polar core layer with
Asn75 described in Section 3.2.

### Asn75 (full-length run)

13.6 % mean relative accessibility — buried, confirming its assignment as the
only buried polar position in the core.

---

## Caveats

The models were not energy-minimised, so side-chain rotamers are as generated.
This affects absolute values more than the buried/exposed classification.

Every truncation creates an artificial terminus. Residues within about seven
positions of a cut should not be reported from that run. To test the cysteines
free of this effect, a run with `--keep 30 110` would place the cut 40 residues
away from Cys71 and Cys72.

Accessibility is computed on the isolated trimer. In the virion A29 is
membrane-associated and bound to A18 and A28, so absolute exposure in vivo will
differ; the measurement tests the internal consistency of the model, not the
biological state.

---

## Files

| | |
| --- | --- |
| `ensemble_sasa.py` | Analysis script |
| `full_length/` | Residues 1–110 — core, cysteines, Asn75 |
| `truncated_44-110/` | Residues 44–110 — A18-binding surface |
| `absolute_area_A29_s2m0_pymol.txt` | Preliminary single-model measurement; superseded |
| `surface_check_A29_s2m0.pse` | PyMOL session used for visual inspection |

Each run folder contains `sasa_summary*.tsv` (per residue, per chain),
`sasa_per_model*.tsv` (every model), `sasa*.xlsx`, and the relabelled PDB files
the measurements were taken from.
