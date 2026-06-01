[![Build](https://github.com/johan-stph/PureCLIP/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/johan-stph/PureCLIP/actions/workflows/build.yml) [![GitHub release](https://img.shields.io/github/release/johan-stph/PureCLIP.svg)](https://github.com/johan-stph/PureCLIP/releases/latest)

PureCLIP is a tool to detect protein-RNA interaction footprints from
single-nucleotide CLIP-seq data, such as iCLIP and eCLIP.

> This is a fork of [skrakau/PureCLIP](https://github.com/skrakau/PureCLIP)
> with CI, platform packaging, Apple Silicon support, and ongoing maintenance
> by [@johan-stph](https://github.com/johan-stph).

---

## Installation

See **[INSTALLATION.md](INSTALLATION.md)** — covers pre-built release binaries,
building from source (macOS & Linux), and upcoming package manager support.

---

## Galaxy: use PureCLIP online

PureCLIP has been integrated into the European Galaxy server
https://usegalaxy.eu/, an open, web-based platform for accessible,
reproducible, and transparent computational biological research and is available
[here](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/pureclip/pureclip/1.0.4)
(currently not using latest PureCLIP version).

Thanks to the Freiburg Galaxy Team!

---

## Quick start

```bash
# Create a config file (edit paths to your data)
cp pureclip_defaults.toml my_run.toml
# Edit: set bam, bai, genome paths

# Run
pureclip -c my_run.toml

# Override specific settings on the command line
pureclip -c my_run.toml --threads 8 --bandwidth 100
```

---

## Usage

PureCLIP has many parameters, but most have sensible defaults.  Parameters
fall into three tiers:

### 🔄 Change every run

| Parameter | Config key | CLI flag | Description |
|-----------|-----------|----------|-------------|
| Input BAM  | `bam`       | `-i`     | One or more BAM files (array for replicates) |
| BAM index  | `bai`       | `-bai`   | One BAI per BAM |
| Reference  | `genome`    | `-g`     | FASTA reference genome |
| Output prefix | `output_prefix` | `-o` | Prefix for output files |
| Chromosomes to learn on | `learn_on` | `-iv` | e.g. `"chr1;chr2;chr3"` — reduces runtime |
| Threads    | `threads`   | `-nt`    | Parallel threads for learning (default: 1) |

### 🎯 Adjust per protein

| Parameter | Config key | CLI flag | Default | Notes |
|-----------|-----------|----------|---------|-------|
| Bandwidth | `bandwidth` | `-bdw` | 50 | Larger = smoother enrichment signal |
| N bandwidth | `bandwidth_n` | `-bdwn` | same as bdw | For estimating binomial N |
| Binding mode | `binding_characteristics` | `-bc` | 0 | 0=short motifs, 1=larger clusters |
| Binomial p init | `bin1_p_init` / `bin2_p_init` | `-b1p` / `-b2p` | 0.01 / 0.15 | Crosslink probability priors |
| Scoring scheme | `scoring_scheme` | `-st` | 0 | 0=log-ratio, 1=crosslink, 2=enrichment, 3=balanced |
| Merge distance | `merge_distance` | `-dm` | 8 bp | Max gap to merge sites into regions |

### 🔒 Rarely changed

| Parameter | Config key | Default | Description |
|-----------|-----------|---------|-------------|
| `max_iter_baumwelch` | `max_iter_baumwelch` | 50 | Outer Baum-Welch iterations |
| `max_iter_simplex` | `max_iter_simplex` | 500 | GSL simplex iterations |
| `g1_k_min` / `g2_k_max` | gamma shape bounds | 1.0 / 10.0 | Shape constraints |
| `gamma_k_convergence` | `gamma_k_convergence` | 1e-4 | Convergence tolerance |
| `min_trans_prob_crosslink` | `min_trans_prob_crosslink` | 0.0001 | Floor for transition probability |
| `verbosity` | `verbosity` | 1 | 0=quiet, 2=verbose, 3=debug |

See [`pureclip_defaults.toml`](pureclip_defaults.toml) for all available
parameters with descriptions and defaults.

**CLI flags override TOML values.**  This lets you keep a base config and
override per-run:

```bash
# Base config for PUM2 experiments
pureclip -c pum2_base.toml

# Same base, but 4× faster with 8 threads
pureclip -c pum2_base.toml -nt 8

# Test a wider bandwidth
pureclip -c pum2_base.toml --bandwidth 100
```

---

## Documentation

Please have a look at PureCLIP's [documentation](http://pureclip.readthedocs.io/en/latest/).

---

## Citation

Krakau S, Richard H, Marsico A: PureCLIP: Capturing target-specific
protein-RNA interaction footprints from single-nucleotide CLIP-seq data.
Genome Biology 2017; 18:240; https://doi.org/10.1186/s13059-017-1364-2
