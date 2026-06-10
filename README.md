# telomere

Tools for locating telomeric repeat motifs in genome assemblies and summarizing telomere-associated windows, developed in the context of the [Vertebrate Genomes Project (VGP)](https://vertebrategenomesproject.org/) / [Darwin Tree of Life](https://www.darwintreeoflife.org/) pipelines at the Wellcome Sanger Institute.

It has been run on **more than 3,000** Darwin Tree of Life assemblies to support genome curation workflows.

## Author

**Sergey Koren** — [github.com/skoren](https://github.com/skoren) · [sergey.koren@nih.gov](mailto:sergey.koren@nih.gov)

## Maintainer

**Yumi Sims** — [yy5@sanger.ac.uk](mailto:yy5@sanger.ac.uk)

## Overview

The pipeline:

1. **`find_telomere`** — Scans a FASTA file for a telomeric repeat (default **vertebrate**: `TTAGGG`) on both strands and reports runs of the motif and its reverse complement.
2. **`SizeFasta`** — Emits per-scaffold sequence lengths (`.lens`).
3. **`sdust`** (via minimap2) — Low-complexity / dust masking (`.sdust`).
4. **`FindTelomereWindows`** — Slides windows along each scaffold and flags regions where telomeric hits exceed a density threshold (optionally adjusted by sequencing identity).
5. **`FindTelomereBreaks`** — Combines lens, dust, and telomere hit tables to refine candidate telomere stretches (filtering short hits and overlap with low-complexity runs).

Shell wrappers tie these steps together for batch and cluster use.

## Requirements

- **C++** compiler (`g++`) for `find_telomere.c`
- **zlib** development headers/library (for reading `.gz` FASTA in `find_telomere`)
- **Java** (JDK) to compile/run the `.java` tools and build `telomere.jar`
- **minimap2** (for `sdust` — `find_telomere.sh` uses `module load minimap2`)
- **bedtools** (for `telomere_analysis.sh`)

## Build

From the repository root:

```bash
chmod +x build.sh find_telomere.sh telomere_analysis.sh _submit_telomere.sh
./build.sh
```

This runs `javac *.java`, packs **`telomere.jar`**, and compiles **`find_telomere`** from `find_telomere.c`.

## Usage

### `find_telomere` (standalone)

```text
./find_telomere <assembly.fasta|assembly.fasta.gz> [repeat]
```

- Second argument is optional; default repeat is **`TTAGGG`**.
- Input FASTA can be plain text or gzip-compressed (detected by gzip magic bytes; not just the filename extension).
- Matching is case-insensitive (both input sequence and motif are normalized to uppercase).
- Output lines to **stdout** are tab-separated: scaffold name, length, strand (`0` = forward motif, `1` = reverse complement), start, end, run length.
  - Scaffold names are printed **without** a leading `>` from the FASTA header.
- In addition, `find_telomere` writes **BED6** files alongside the input FASTA:
  - `<input>.fwd.telomere.bed` (strand `+`, forward motif)
  - `<input>.rev.telomere.bed` (strand `-`, reverse-complement motif)

### `find_telomere.sh` (full preprocessing)

```bash
./find_telomere.sh /path/to/assembly.fasta
```

| Output | Description |
|--------|-------------|
| `prefix.telomere` | Motif hits (from `find_telomere`, reformatted) |
| `prefix.sdust` | `sdust` low-complexity intervals |
| `prefix.lens` | Scaffold sizes (`SizeFasta`) |
| `prefix.windows` | Telomere density windows (`FindTelomereWindows`, example args `99.9 0.1`) |
| `prefix.breaks` | Refined breaks (`FindTelomereBreaks`) |

### Java entry points (classpath)

```bash
java -cp telomere.jar SizeFasta <assembly.fa> > prefix.lens
java -cp telomere.jar FindTelomereWindows prefix.telomere <identity_percent> [threshold] > prefix.windows
java -cp telomere.jar FindTelomereWindows --split prefix.telomere <identity_percent> [threshold] > prefix.windows
java -cp telomere.jar FindTelomereBreaks prefix.lens prefix.sdust prefix.telomere > prefix.breaks
```

- **`FindTelomereWindows`**: identity is typically given as e.g. `99.9` (percent); threshold defines minimum window occupancy and is scaled by identity (`threshold * identity^6`).
- **Output format**: bedGraph-like, tab-separated `chrom  start  end  score`, where `score = ceil(1000 * density)` (capped at 1000).
- **Default mode** (no `--split`): combines both strands in one window track (printed to stdout); default base threshold is `0.4`.
- **Split mode** (`--split`): writes two stranded output files (default base threshold `0.05`):
  - if input is `prefix.telomere`: `prefix.fwd.windows` and `prefix.rev.windows`
  - otherwise: `<input>.fwd.windows` and `<input>.rev.windows`
  - scaffold names are unchanged in both files (no `_fwd` / `_rev` suffix)
  - combined (both-strand) windows are still emitted to **stdout** (so you can redirect to `prefix.windows` as in the example above)
- **`FindTelomereBreaks`**: filters telomere runs shorter than 24 bp and uses dust masks to avoid calling repeats inside low-complexity-only regions.

### `telomere_analysis.sh` (BED / assembly QC)

Higher-level script that runs `find_telomere.sh`, merges windows with **bedtools**, intersects scaffold ends, and optionally filters by a chromosome-assignment CSV. Intended for Slurm (`$SLURM_CPUS_PER_TASK`) and requires **`$VGP_PIPELINE`**, **bedtools**, and **Java**. See script comments for positional arguments (genome label, threshold, end distance, FASTA, optional CSV).

## Repository layout

| Path | Role |
|------|------|
| `find_telomere.c` | Motif scanner (C++) |
| `find_telomere`, `find_telomere.sh` | Binary and pipeline wrapper |
| `SizeFasta.java`, `FindTelomereWindows.java`, `FindTelomereBreaks.java`, `Utils.java` | Java analysis |
| `telomere.jar` | Pre-built classes (rebuild with `build.sh`) |
| `telomere_analysis.sh` | BED/summary workflow |
| `_submit_telomere.sh` | Cluster submission helper |
| `LICENSE` | BSD-3-Clause |

## Want to contribute?

Contributions are welcome. Open an issue to discuss a bug, feature idea, or larger change, or send a pull request. If you are not sure where to start, contact the [maintainer](#maintainer). Please follow the [code of conduct](#code-of-conduct) in issues and pull requests.

## License

This project is licensed under the **BSD 3-Clause License** — see [`LICENSE`](LICENSE).

## Citation

If you use this software, please cite the relevant VGP / Darwin Tree of Life publications and acknowledge the [Wellcome Sanger Institute](https://www.sanger.ac.uk/) and project funders as appropriate for your work.

## Code of conduct

This repository includes the Darwin Tree of Life **Code of Professional Conduct**; contributors are expected to follow it in issues, pull requests, and discussions.
