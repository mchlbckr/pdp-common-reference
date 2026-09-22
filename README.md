# Partial Dependence Beyond the Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22904255.svg)](https://doi.org/10.5281/zenodo.22904255)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Research code, simulation pipeline, and Quarto manuscript sources for
*Partial Dependence Beyond the Data: Support Diagnostics and Common-Reference
Feature Effects*.

Partial-dependence plots may evaluate a fitted model at feature combinations
that the data represent poorly. This project turns that problem into an
explicit support audit and characterizes when a single reference population
remains valid over a displayed grid.

## Scope

The repository contains:

* in-support-coverage (ISC), hard common-support partial-dependence (CSPD), and
  relaxed common-reference (RCPD) estimators;
* an exact empirical solver for interval validity, together with ALE,
  conditional-curve, and subgroup-PDP baselines;
* an independent audit split that certifies simultaneous query validity, with
  explicit abstention when no candidate reference can be certified;
* unit tests, a `targets` pipeline, and a JMLR-formatted Quarto manuscript.

The main experiment matrix crosses six prediction-model classes, three
correlations, three grid widths, and five calibrated validity rules. Further
experiments study the validity–fidelity trade-off, pointwise-trimming
composition, competing constraint formulations, grid resolution,
validity-learning error, and application stability.

The manuscript reports two applications: a confirmatory **California Housing**
analysis over 30 independently seeded splits, and a controlled **Bank
Marketing** off-support stress test. Wine Quality, Bank Marketing, and Online
Shoppers additionally serve as held-out negative controls; see
[`results/extended-results-index.md`](results/extended-results-index.md).

## Layout

| Path | Contents |
|---|---|
| `R/` | Estimators and supporting functions. |
| `simulations/` | Data-generating processes and simulation studies. |
| `exploration/` | Development scripts that the confirmatory California analysis (`s33`) loads into a private environment. |
| `tests/` | Unit tests. |
| `paper/article.qmd` | Main manuscript and mathematical appendices. |
| `paper/supplement.qmd` | Empirical online supplement. |
| `paper/_extensions/jmlr/` | Local Quarto format using the JMLR style file. |
| `paper/references.bib` | Manuscript bibliography. |
| `results/` | Figure panels and machine-readable result objects. |
| `_targets.R` | Reproducible computational pipeline. |
| `_targets/` | Stored target objects, so the manuscript renders without rerunning the pipeline. |
| `data/README.md` | Acquisition commands, source URLs, licences, and checksums. |

## Data

No raw data files are redistributed here. Every dataset is obtained from its
original source; [`data/README.md`](data/README.md) gives the download
commands, licences, and SHA-256 checksums needed to reconstruct the exact
inputs.

## Reproduction

Restore the locked R environment, regenerate all targets, and run the tests:

```r
renv::restore()
targets::tar_make()
testthat::test_dir("tests/testthat")
```

Rendering the manuscript does **not** require rerunning the pipeline, because
the target objects under `_targets/` are included:

```sh
quarto render paper
```

This creates `paper/_output/pdp-common-reference-jmlr.pdf` and
`paper/_output/pdp-common-reference-jmlr-supplement.pdf`.

The rendered PDFs and the intermediate TeX files are **not** tracked here: this
repository is the computational archive, not a copy of the article. The
manuscript sources are included because they carry the inline computations that
produce the reported numbers from the stored target objects.

See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for the full checklist,
including expected runtimes and the environment the results were produced in.

## Citation

If you use this code or the accompanying results, please cite the article and
this archive; see [`CITATION.cff`](CITATION.cff).

The archive is deposited on Zenodo:

* Concept DOI [10.5281/zenodo.22904255](https://doi.org/10.5281/zenodo.22904255)
  always resolves to the most recent version and is the one to cite in general.
* Version DOI [10.5281/zenodo.22904256](https://doi.org/10.5281/zenodo.22904256)
  pins release v1.0.0 exactly.

## Licence

MIT, see [`LICENSE`](LICENSE). The JMLR style file `paper/jmlr2e.sty` is
distributed by the Journal of Machine Learning Research under its own terms and
is included only to make the manuscript reproducible.
