[![Build](https://github.com/johan-stph/PureCLIP/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/johan-stph/PureCLIP/actions/workflows/build.yml) [![GitHub release](https://img.shields.io/github/release/johan-stph/PureCLIP.svg)](https://github.com/johan-stph/PureCLIP/releases/latest)

PureCLIP is a tool to detect protein-RNA interaction footprints from single-nucleotide CLIP-seq data, such as iCLIP and eCLIP.

# Installation

You can install PureCLIP from the release tarballs or from source.

## Release Binaries

Pre-built binaries for **macOS arm64**, **Linux x86\_64** and **Linux arm64** are available on the [Releases](https://github.com/johan-stph/PureCLIP/releases/latest) page.
Download the archive for your platform and extract it:

    $ tar -xzf pureclip-<version>-<platform>.tar.gz
    $ ./pureclip --version

> **TODO:** publish packages to [Bioconda](http://bioconda.github.io) and [Homebrew](https://brew.sh) so that PureCLIP can be installed via
> `conda install pureclip` / `brew install pureclip`.

# Galaxy: use PureCLIP online

PureCLIP has also been integrated into the European Galaxy server https://usegalaxy.eu/, an open, web-based platform for accessible, reproducible, and transparent computational biological research and is available [here](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/pureclip/pureclip/1.0.4) (currently not using latest PureCLIP version).

Thanks to the Freiburg Galaxy Team!

# Installation

For full installation instructions — including **Apple Silicon (arm64)** native
builds and Linux — see **[INSTALLATION.md](INSTALLATION.md)**.

### Quick start (Linux / macOS)

    $ git clone https://github.com/johan-stph/PureCLIP.git
    $ cd PureCLIP
    $ mkdir build && cd build
    $ cmake ../src -DCMAKE_BUILD_TYPE=Release
    $ make -j$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)

> **macOS (Apple Silicon):** install `libomp` first (`brew install cmake gsl libomp`).
> The build system detects Apple Clang automatically — no extra flags needed.
> See [INSTALLATION.md](INSTALLATION.md) for details.

Requirements

 - C++14 compliant compiler (GCC ≥ 5 or Clang ≥ 3.4)
 - GSL
 - cmake 3.0 or newer
 - OpenMP (`libomp` via Homebrew on macOS)


# Documentation

Please have a look at PureCLIPs [documentation](http://pureclip.readthedocs.io/en/latest/).

# Citation

Krakau S, Richard H, Marsico A: PureCLIP: Capturing target-specific protein-RNA interaction footprints from single-nucleotide CLIP-seq data. Genome Biology 2017; 18:240; https://doi.org/10.1186/s13059-017-1364-2
