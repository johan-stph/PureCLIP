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

## Configuration file

Parameters can be given in a TOML file instead of on the command line:

```bash
pureclip -c my_config.toml
```

`pureclip_defaults.toml` documents every supported key and can be copied as a
starting point. Precedence is **command line > config file > built-in
defaults**, so a config can hold your standing setup while flags override it
per run:

```bash
pureclip -c base.toml -bw 100 -nt 8
```

`output_prefix` derives the three output paths (`<prefix>_sites.bed`,
`<prefix>_regions.bed`, `<prefix>_params.txt`), so `-o`/`-or`/`-p` can be
omitted. Inputs (`bam`, `bai`, `genome`) may come from either source; unknown
keys are ignored so a config written for a newer version still loads.

## Testing

See **[TESTING.md](TESTING.md)** — `make test` runs a fast correctness check;
`make test-all` adds the real-data tier.

## Galaxy: use PureCLIP online

PureCLIP has been integrated into the European Galaxy server
https://usegalaxy.eu/, an open, web-based platform for accessible,
reproducible, and transparent computational biological research and is available
[here](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/pureclip/pureclip/1.0.4)
(currently not using latest PureCLIP version).

Thanks to the Freiburg Galaxy Team!

---

## Documentation

Please have a look at PureCLIP's [documentation](http://pureclip.readthedocs.io/en/latest/).

---

## Citation

Krakau S, Richard H, Marsico A: PureCLIP: Capturing target-specific
protein-RNA interaction footprints from single-nucleotide CLIP-seq data.
Genome Biology 2017; 18:240; https://doi.org/10.1186/s13059-017-1364-2
