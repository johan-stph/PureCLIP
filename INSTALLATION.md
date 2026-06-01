# Installing PureCLIP

This guide covers installing PureCLIP from source on **Linux** and **macOS**
(including Apple Silicon / arm64), as well as via Bioconda.

---

## Table of Contents

- [Bioconda (easiest)](#bioconda-easiest)
- [Linux – build from source](#linux--build-from-source)
- [macOS – build from source](#macos--build-from-source)
  - [Intel (x86\_64)](#intel-x86_64)
  - [Apple Silicon (arm64)](#apple-silicon-arm64)
- [Verify the build](#verify-the-build)
- [Quick sample run](#quick-sample-run)

---

## Bioconda (easiest)

```bash
conda install -c bioconda pureclip
```

> Requires an active [Bioconda channel](https://bioconda.github.io).
> Pre-built bottles are available for Linux (x86\_64) and macOS (x86\_64).
> For Apple Silicon use the native build from source described below.

---

## Linux – build from source

### Prerequisites

| Tool | Version | Install |
|------|---------|---------|
| CMake | ≥ 3.0 | `apt install cmake` / `dnf install cmake` |
| GCC | ≥ 5 (C++14) | `apt install g++` |
| GSL | any | `apt install libgsl-dev` |
| zlib | any | `apt install zlib1g-dev` |

SeqAn 2.2.0 and Boost 1.64.0 are **downloaded automatically** by CMake on
first configure.

### Build

```bash
git clone https://github.com/skrakau/PureCLIP.git
cd PureCLIP
mkdir build && cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Binaries are written to `build/pureclip` and `build/winextract`.
Optionally install system-wide:

```bash
sudo make install   # installs to /usr/local/bin by default
```

---

## macOS – build from source

### Prerequisites (all architectures)

Install [Homebrew](https://brew.sh) if you do not have it, then:

```bash
brew install cmake gsl libomp
```

| Package | Purpose |
|---------|---------|
| `cmake` | Build system |
| `gsl` | GNU Scientific Library |
| `libomp` | OpenMP runtime for Apple Clang (both Intel and arm64) |

### Intel (x86\_64)

```bash
git clone https://github.com/skrakau/PureCLIP.git
cd PureCLIP
mkdir build && cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
make -j$(sysctl -n hw.logicalcpu)
```

### Apple Silicon (arm64)

Apple's bundled `clang` does not ship with an OpenMP runtime.
The CMakeLists.txt in this repo detects Apple Clang automatically and picks up
`libomp` from Homebrew — **no extra flags are needed**.

```bash
git clone https://github.com/skrakau/PureCLIP.git
cd PureCLIP
mkdir build && cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
make -j$(sysctl -n hw.logicalcpu)
```

> **What happens under the hood:**  CMake detects `Apple Clang`, calls
> `brew --prefix libomp`, and sets `-Xpreprocessor -fopenmp` together with the
> correct include / library paths automatically.
> If `libomp` is not found cmake will abort with a clear install instruction.

#### Troubleshooting

| Problem | Fix |
|---------|-----|
| `brew: command not found` during cmake | Install Homebrew or set `LIBOMP_PREFIX` manually: `cmake ../src -DOpenMP_omp_LIBRARY=/path/to/libomp.dylib` |
| `ld: library 'omp' not found` | Run `brew install libomp` and re-run cmake |
| SeqAn / Boost download fails | Check internet access; or pass `-DSEQAN_ROOT=` / `-DBOOST_ROOT=` pointing to local copies |

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
