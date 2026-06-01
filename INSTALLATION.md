# Installing PureCLIP

This guide covers all installation methods: pre-built binaries, build from
source, and package managers (coming soon).

---

## Table of Contents

- [Pre-built release binaries](#pre-built-release-binaries)
- [Build from source](#build-from-source)
  - [macOS (Intel + Apple Silicon)](#macos-intel--apple-silicon)
  - [Linux](#linux)
- [Verify the build](#verify-the-build)
- [Quick sample run](#quick-sample-run)
- [Coming soon](#coming-soon)

---

## Pre-built release binaries

Pre-built binaries for **macOS arm64**, **Linux x86\_64** and **Linux arm64** are
available on the [Releases](https://github.com/johan-stph/PureCLIP/releases/latest) page.

```bash
# Download the archive for your platform
tar -xzf pureclip-<version>-<platform>.tar.gz
./pureclip --version
```

---

## Build from source

### macOS (Intel + Apple Silicon)

```bash
# Prerequisites
brew install cmake gsl libomp

# Build
git clone https://github.com/johan-stph/PureCLIP.git
cd PureCLIP
mkdir build && cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
make -j$(sysctl -n hw.logicalcpu)
```

| Package | Purpose |
|---------|---------|
| `cmake`  | Build system |
| `gsl`    | GNU Scientific Library |
| `libomp` | OpenMP runtime for Apple Clang |

> **Apple Silicon:** The build system detects Apple Clang automatically, calls
> `brew --prefix libomp`, and sets `-Xpreprocessor -fopenmp` +
> correct include/library paths. No extra flags needed.

### Linux

```bash
# Prerequisites (Debian/Ubuntu)
sudo apt install cmake g++ libgsl-dev zlib1g-dev

# Build
git clone https://github.com/johan-stph/PureCLIP.git
cd PureCLIP
mkdir build && cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

SeqAn and Boost are **downloaded automatically** by CMake on first configure.

Optionally install system-wide:

```bash
sudo make install   # installs to /usr/local/bin by default
```

---

## Verify the build

```bash
./pureclip --version
./winextract --version
```

---

## Quick sample run

The following uses a preprocessed eCLIP BAM from ENCODE
([Van Nostrand et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27018577))
and the hg19 reference genome to run PureCLIP in basic mode.

> **Disk space:** ~4 GB for the BAM + genome.
> **Runtime:** ~10–20 min depending on CPU count.

```bash
# 1. Download data (requires samtools on PATH)
wget -O aligned.prepro.bam \
    https://www.encodeproject.org/files/ENCFF280ONP/@@download/ENCFF280ONP.bam
samtools view -hb -f 130 aligned.prepro.bam -o aligned.prepro.R2.bam
samtools index aligned.prepro.R2.bam

# 2. Download reference genome
wget -O ref.hg19.fa.gz \
    https://www.encodeproject.org/files/female.hg19/@@download/female.hg19.fasta.gz
gunzip ref.hg19.fa.gz

# 3. Run PureCLIP (learn on chr1–3, 10 threads)
./pureclip \
    -i  aligned.prepro.R2.bam \
    -bai aligned.prepro.R2.bam.bai \
    -g  ref.hg19.fa \
    -iv 'chr1;chr2;chr3;' \
    -nt 10 \
    -o  PureCLIP.crosslink_sites.bed
```

Output files:
- `PureCLIP.crosslink_sites.bed` — individual crosslink sites
- `PureCLIP.crosslink_clusters.bed` — merged binding regions (generated automatically)

For more options see the [full documentation](http://pureclip.readthedocs.io/en/latest/).

---

## Coming soon

| Method | Install command | Status |
|--------|----------------|--------|
| **Homebrew** | `brew install pureclip2` | 🚧 tap ready, formula in review |
| **Bioconda** | `conda install -c bioconda pureclip2` | 🚧 recipe scaffolded, submission pending |

Both will auto-resolve dependencies and keep PureCLIP up to date with each release.
