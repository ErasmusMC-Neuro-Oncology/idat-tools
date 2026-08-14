# idat-tools

[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.9-blue?logo=python&logoColor=white)](https://www.python.org/)
[![License: GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.xcrm.2026.102682-orange)](https://doi.org/10.1016/j.xcrm.2026.102682)

**idat-tools** is a command-line toolkit for reading, inspecting, and manipulating Illumina IDAT files — the raw intensity files produced by Illumina methylation arrays (e.g. EPIC, 850k, 450k).

Supported operations:

| Command | Description |
|---------|-------------|
| `idat-tools view` | Inspect IDAT metadata and probe intensity table |
| `idat-tools mix`  | Create in-silico mixed samples at a controlled ratio |
| `idat-tools subtract` | Remove the contribution of one sample from another (e.g. normal from tumour) |

---

## Installation

Requires Python ≥ 3.9.

```bash
git clone https://github.com/yhoogstrate/idat-tools.git
cd idat-tools
python3 -m venv .venv
source .venv/bin/activate
pip install .
idat-tools --version
```

---

## `idat-tools view`

Prints the IDAT metadata header followed by a summary of the probe intensity table.

```bash
idat-tools view [OPTIONS] IDAT_FILE

Options:
  -n INTEGER        Number of rows to display  [default: 10]
  --no-header       Print only the probe table, omitting the metadata header
  --version         Show version and exit
  --help            Show this message and exit
```

### Example

```bash
idat-tools view 207513420108_R01C01_Grn.idat
```

Output:

```
# array_n_probes:       1052641
# total intensity:      1398487731
# manifest:             ''
# manifest (old style): ''
# unknown #1:           [1][0][0][0]
# sample id:            ''
# description:          ''
# plate:                ''
# well:                 ''
# unknown #2:           ''
# run info:
# 1. [04/09/2023 2:07:30 PM] [Decoding] [...] [AutoDecode] [3.0.1.0]
# 2. [4/9/2024 3:54:42 PM]   [Scan]     [...] [iScan Control Software] [4.2.1.729]
# ...

IDAT v3: 207513420108_R01C01 (R/G: 0, BeadChip 8x5)
         probe_ids  probe_std_devs  probe_mean_intensities  probe_n_beads  probe_mid_block
0          1600101             264                    1103             18          1600101
1          1600111             185                     956             11          1600111
...            ...             ...                     ...            ...              ...
1052639   99810990             424                    1815             11         99810990
1052640   99810992              97                     400             12         99810992

[1052641 rows x 5 columns]
```

Use `-n` to control how many rows are shown, and `--no-header` to suppress the metadata block (useful for piping into downstream tools):

```bash
idat-tools view -n 20 sample_Grn.idat
idat-tools view --no-header sample_Grn.idat | head -5
```

---

## `idat-tools mix`

Creates an artificial mixed IDAT file by blending probe intensities from two input files at a configurable ratio. Useful for generating in-silico tumour purity gradients or benchmarking deconvolution methods.

```bash
idat-tools mix [OPTIONS] IDAT_FILE_REFERENCE IDAT_FILE_MIXED_IN IDAT_FILE_OUTPUT

Options:
  -r, --mix-ratio FLOAT RANGE  Fraction of the mixed-in file  [default: 0.5; 0<=x<=1]
  --geometric-mean             Mix in log space instead of linear space
  --help                       Show this message and exit
```

### Mixing models

#### Linear average (default)

Intensities are combined as a **weighted linear average**:

```
I_out = (1 − r) × I_ref  +  r × I_mix
```

This is the physically correct model for a mixture of cell populations: each bead reports a signal that averages over the cells in its vicinity, so the expected intensity is a linear combination of the two population signals weighted by their cell fractions.

#### Geometric mean (`--geometric-mean`)

Mixes in log space:

```
I_out = I_ref^(1−r) × I_mix^r
      = exp( (1−r) × log(I_ref) + r × log(I_mix) )
```

Appropriate when intensities are modelled as log-normally distributed and multiplicative noise dominates. Equal steps in log space correspond to equal fold-changes. Note that for a physical cell mixture, the linear model is more accurate; the geometric mean systematically underweights high-intensity probes relative to low-intensity ones. Probes with zero intensity are clipped to 1 before log-transformation.

### Examples

```bash
# 75% reference, 25% mixed-in (linear, default)
idat-tools mix -r 0.25 \
    207513420108_R01C01_Grn.idat \
    207513420108_R02C01_Grn.idat \
    207513420108_R99C01_Grn.idat

# Same ratio, geometric mean
idat-tools mix -r 0.25 --geometric-mean \
    207513420108_R01C01_Grn.idat \
    207513420108_R02C01_Grn.idat \
    207513420108_R98C01_Grn.idat
```

The output filename should follow standard Sentrix ID nomenclature (`<barcode>_<position>_<channel>.idat`); if it does not, a random barcode and chip label are generated and a warning is printed.

---

## `idat-tools subtract`

Removes the estimated contribution of one sample from another — for example subtracting a matched normal from a tumour sample to approximate the profile of the tumour cells alone. It is the algebraic inverse of `idat-tools mix`.

```bash
idat-tools subtract [OPTIONS] IDAT_FILE_OBSERVED IDAT_FILE_SUBTRACTED IDAT_FILE_OUTPUT

Options:
  -r, --mix-ratio FLOAT RANGE  Estimated fraction of the subtracted file present
                               in the observed file  [default: 0.5; 0<=x<1]
  --help                       Show this message and exit
```

### Subtraction model

The observed sample is assumed to follow the same linear mixing model that `mix` implements:

```
I_obs = (1 − r) × I_pure  +  r × I_sub
```

Solving for the signal of interest gives:

```
I_pure = ( I_obs  −  r × I_sub ) / (1 − r)
```

Here `r` is the estimated fraction of the subtracted sample in the observed file — for a tumour biopsy, one minus the tumour purity. The division by `1 − r` rescales the remaining signal back to full intensity, so subtraction is an exact inverse of mixing rather than merely a difference.

Consequences worth being aware of:

- **`r` must be below 1.** At `r = 1` the observed file consists entirely of the subtracted sample and there is nothing left to recover; the value is rejected.
- **Negative results are clipped to 0.** Intensities are stored as unsigned 16-bit integers, so a probe where the subtracted signal exceeds the observed signal cannot be represented. Such probes are set to 0 and the number of affected probes is reported as a warning. A high proportion is a sign that `r` is set too high.
- **Standard deviations are propagated, not subtracted.** Removing a component from the mean does not remove its noise, so the uncertainties are combined in quadrature and rescaled: `sqrt( sd_obs² + (r × sd_sub)² ) / (1 − r)`.
- **Bead counts are carried over unchanged.** `probe_n_beads` is a physical count of beads on the array and is unaffected by the arithmetic.
- There is no `--geometric-mean` equivalent; subtraction is defined only for the linear model.

### Examples

```bash
# Remove an estimated 25% normal contamination from a tumour sample
idat-tools subtract -r 0.25 \
    207513420108_R01C01_Grn.idat \
    207513420108_R02C01_Grn.idat \
    207513420108_R97C01_Grn.idat
```

Because `subtract` inverts `mix`, applying both at the same ratio returns the original reference file, up to the rounding of the two intermediate integer casts:

```bash
idat-tools mix      -r 0.25 ref_Grn.idat other_Grn.idat mixed_Grn.idat
idat-tools subtract -r 0.25 mixed_Grn.idat other_Grn.idat recovered_Grn.idat
# recovered_Grn.idat matches ref_Grn.idat to within 1 intensity unit
```

The rounding error grows with `r`, since the `1 / (1 − r)` rescaling amplifies it: it stays within 1 unit up to `r = 0.5` and reaches roughly 5 units at `r = 0.9`.

As with `mix`, the output filename should follow standard Sentrix ID nomenclature.

### Interpretation

Subtraction operates on **raw probe intensities**, which still contain background fluorescence, dye bias and array-specific batch effects. The linear mixing model is therefore only approximately true of real data: the round-trip against `mix` is exact by construction, but subtracting a real matched normal from a real tumour will not yield a perfectly pure tumour profile. Treat the output as an enrichment of the signal of interest rather than a clean deconvolution, and normalise downstream as usual.

---

## Testing

The test suite downloads a small set of public IDAT files from GEO (cached in `cache/`) and exercises all three commands:

```bash
source .venv/bin/activate
make test
```

To run against your own files instead:

```bash
IDAT_REF=my_R01C01_Grn.idat IDAT_MIX=my_R02C01_Grn.idat bash tests.sh
```

Both files must originate from the same array type; mixing an EPIC with a 450k file is rejected.

---

## Citation

If you use idat-tools in your research, please cite:

> Youri Hoogstrate, Santoesha A. Ghisai, Levi van Hijfte, Rania Head, Iris de Heer,
> Marta Padovan, Maurice de Wit, Wies R. Vallentgoed, Angelo Dipasquale,
> Maarten M.J. Wijnenga, Bas Weenink, Rosa Luning, Sybren L.N. Maas,
> Adela Brzobohata, Michael Weller, Tobias Weiss, Maximilian J. Mair,
> Anna S. Berghoff, Adelheid Wöhrer, Albert Jeltsch, Johan A.F. Koekkoek,
> Hans M. Hazelbag, Mathilde C.M. Kouwenhoven, Yongsoo Kim, Bart A. Westerman,
> Bauke Ylstra, Johanna M. Niers, Kevin C. Johnson, Frederick S. Varn,
> Roel G.W. Verhaak, Mustafa Khasraw, Martin J. van den Bent, Pieter Wesseling,
> Pim J. French.
> **TET CpG sequence-context-specific DNA demethylation shapes progression of IDH-mutant gliomas.**
> *Cell Reports Medicine*, Volume 7, Issue 3, 2026, 102682. ISSN 2666-3791.
> https://doi.org/10.1016/j.xcrm.2026.102682

---

## License

idat-tools is released under the [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.html).
Copyright (C) 2024 Youri Hoogstrate.
