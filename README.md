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
