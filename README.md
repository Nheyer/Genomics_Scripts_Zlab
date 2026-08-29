# Genomics_Scripts_Zlab

Small genomics command line tools used in the Zabel Lab at Colorado State
University, plus a few external tools vendored as submodules.

## Contents

| Path | What it is |
|---|---|
| `CppSrc/` | C++ source only, one directory per program |
| `CppSrc/MSA-to-consensus/MSA_fasta_to_Consensus.cpp` | `MSA-to-consensus` — collapses a multiple sequence alignment into a single consensus sequence |
| `CppSrc/Enzyme-digest/Restriction_Enzyme_Digest.cpp` | `Enzyme-digest` — in silico restriction digest, fragment sizes and an ASCII gel |
| `CppSrc/Enzyme-digest/enzyme_data.hpp` | IUPAC codes, complement table and ladder — lookup tables only, no logic |
| `CppSrc/Enzyme-digest/enzyme_csv_embedded.hpp.in` | Template CMake fills with the enzyme CSV at configure time; required to build |
| `Data/restriction_enzymes.csv` | The enzyme database, in NEB `^`/`_` notation. **The source of truth** — see [Enzyme data](#enzyme-data) |
| `LICENSE.restriction-digest` | MIT notice for the code `Enzyme-digest` was ported from, see [License](#license) |
| `Tests/test_consensus.cpp` | CUnit unit tests, see [Running the tests](#running-the-tests) |
| `Tests/test_digest.cpp` | CUnit unit tests for `Enzyme-digest`, same |
| `External_tools/argparse` | [p-ranav/argparse](https://github.com/p-ranav/argparse), header only CLI parsing (build dependency) |
| `External_tools/cunit` | [cunity/cunit](https://gitlab.com/cunity/cunit), unit test framework (build dependency) |
| `Data/test_files/` | Known truth fixtures, see [Test files](#test-files) |

## Building

Requires CMake >= 3.15 and a C++17 compiler.

All dependencies are vendored as submodules, so clone recursively:

```bash
git clone --recurse-submodules https://github.com/Nheyer/Genomics_Scripts_Zlab.git
cd Genomics_Scripts_Zlab
```

If you already cloned without `--recurse-submodules`:

```bash
git submodule update --init --recursive
```

Then configure and build:

```bash
cmake .
cmake --build .
```

The binaries land in `bin/MSA-to-consensus` and `bin/Enzyme-digest`; build one on
its own by naming it as the `--target`. The build uses the vendored
`External_tools/argparse` headers, not any system wide install, so it does not
matter whether argparse is installed on the machine.

Configuring also embeds `Data/restriction_enzymes.csv` into a generated
header under `generated/`, which is why `Enzyme-digest` needs no data file at
runtime. See [Enzyme data](#enzyme-data).

### Installing

```bash
cmake --install .
```

(`make install` does the same thing.) That puts both tools on your PATH and
drops the README and the licence files next to them:

```
<prefix>/bin/MSA-to-consensus
<prefix>/bin/Enzyme-digest
<prefix>/share/doc/Genomics_Scripts_Zlab/
```

The installed docs include `LICENSE.restriction-digest`. That notice has to
travel with the binary, so do not drop it from the install list — see
[License](#license).

The prefix defaults to `/usr/local`, which needs `sudo`. To install somewhere
you own instead, no `sudo` required:

```bash
cmake --install . --prefix ~/.local
```

The test binary and the vendored CUnit library are deliberately not installed —
they are only there to build and check the tools.

There is no `uninstall` target, but installing writes `install_manifest.txt`
listing every file it placed, so removing them is:

```bash
xargs rm -f < install_manifest.txt
```

The digest simulator no longer needs installing separately: `Enzyme-digest` is a
C++ port of it and installs with everything else — see
[The digest simulator](#the-digest-simulator).

### The digest simulator

`Enzyme-digest` is a C++ port of the Python restriction digest simulator, taken
from commit `82994c5` on `main` of
[our fork](https://github.com/Nheyer/restriction-enzyme-digest-simulator).
Upstream (`wl5e/restriction-enzyme-digest-simulator`) is archived and read-only,
so the fork's `main` is where that project's work actually lives.

The Python is **not** vendored here. It was carried as a submodule while the
port was being written and has been removed now that the port stands on its
own — `External_tools/` holds only the two build dependencies. Record the
commit above rather than the working copy: that hash is what the port
corresponds to, and it is the provenance behind the attribution in
[License](#license). Removing the copy changes nothing about that obligation.

The two implementations still share a CLI surface, so a differential run is
possible by cloning the Python separately and feeding both the same FASTA. The
one deliberate difference is the program name: ours reports itself as
`Enzyme-digest`, matching the binary, where the Python calls itself
`enzyme_digest.py`. That shows up in the usage line and in the
required-argument error, so those two do not compare byte for byte. Everything
else does, including the digest output and the other error messages.

### Running the tests

`Tests/` unit tests both tools with CUnit:

```bash
cmake --build . --target test_consensus test_digest
./bin/test_consensus
./bin/test_digest
```

or through CTest, which is what CI would use:

```bash
ctest
```

Either way the runner exits non zero if anything fails. The consensus suite
covers `clean_nucliotide`, `encode_nucliotide`, `decode_nucliotide` and
`make_consensus` — ambiguity codes, strict masking, gap handling, and the
accumulator compression path.

Because the tests build the alignments in memory and call `make_consensus`
directly, they can cover things a FASTA fixture makes awkward: an alignment
long enough to force compression is `repeated("G", 1000)` rather than a
thousand line file.

The digest suite leads with invariants, which is what the upstream ROADMAP asks
of a port — they are language independent and they are what catches a sign
error or an off-by-one in the cut coordinate frame:

- **Fragments sum to the sequence length**, for all 262 enzymes, linear and
  circular, on concrete and on ambiguous sequence. Whatever the cuts are, the
  pieces have to add back up to the molecule.
- **The 15-code definite-matching table**, written out by hand rather than
  derived from `IUPAC_BASES`, so changing the matching rule fails here instead
  of being mirrored by a re-derivation of itself.
- **Overhang mirroring** for palindromic enzymes that cut inside their own site,
  which is forced by the symmetry of the molecule and so is independent of
  REBASE as well as of the code.

Beyond that it covers both NEB notations, the CSV `^`/`_` form, both strands,
overlapping sites, origin-spanning circular sites, gap stripping, and the
`definite`/`possible` split.

### Cleaning up

The build happens in source, so it leaves generated files lying around the
repo. To put the tree back the way you found it:

```bash
cmake --build . --target clean
```

That removes the binaries plus everything configuring generated —
`CMakeCache.txt`, `CMakeFiles/`, `Makefile`, `cmake_install.cmake`,
`compile_commands.json`, the CTest leftovers, and what building CUnit dropped
inside `External_tools/cunit` — so `git status` comes back clean and there is
no need for a `.gitignore`. Anything else you happen to have parked in `bin/`
is left alone.

Since this also deletes the `Makefile`, run `cmake .` again before your next
build.

## MSA-to-consensus

Reads a multiple sequence aligned FASTA file and writes a single consensus
sequence. Every input sequence must be the same length (it is an *aligned*
file); mismatched lengths are rejected.

```bash
./bin/MSA-to-consensus -i Data/test_files/TEST_IUPAC_R.msa.fna -a
```

```
Opening: Data/test_files/TEST_IUPAC_R.msa.fna as multiply aligned Nucliotide Sequence input!
Found: 2 sequense(s) in 4 lines
>Consensus_Sequence_of:
RRRRRRRRRRRR
```

Progress messages go to stderr, the FASTA itself to stdout (or to `-o`), so
`> out.fna` captures just the sequence.

### Options

| Flag | Meaning |
|---|---|
| `-i`, `--input` | Input aligned FASTA (required) |
| `-o`, `--output` | Output file, defaults to `stdout` |
| `-a`, `--use-ambiguity` | Collapse disagreeing columns to an IUPAC ambiguity code instead of masking them |
| `-fna`, `--force-fna` | Treat input as nucleotides regardless of file extension |
| `-faa`, `--force-faa` | Treat input as amino acids regardless of file extension |
| `--verbose` | Be loud |
| `-dm`, `--disagree-mask` | **Accepted but not yet implemented** — the mask is always `N` |

Input type is otherwise inferred from the file ending (`.fna` nucleotide,
`.faa` amino acid); anything else needs `-fna`/`-faa`.

### How a column is resolved

For each alignment column:

- **Any `-` in the column wins.** A gap trumps everything else and the
  consensus gets `-` for that position, in both modes.
- **Default (strict) mode:** if every sequence agrees, that base is used;
  any disagreement is masked to `N`.
- **`-a` (ambiguity) mode:** the set of observed bases is collapsed to the
  IUPAC code covering exactly that set, e.g. `A` + `G` becomes `R`. Ambiguity
  codes already present in the input are expanded and merged too, so `R` + `C`
  becomes `V` (`{A,G}` ∪ `{C}` = `{A,C,G}`).

Encoding uses one prime per base (A=3, T=5, C=7, G=11); a column's bases are
multiplied together and the product is matched back to an IUPAC letter. `U` is
read as `T`.

## Enzyme-digest

Simulates a restriction digest: finds where the given enzymes cut a DNA
sequence, reports the fragment sizes, and optionally draws an ASCII gel.

```bash
./bin/Enzyme-digest -f Data/test_files/TEST.fna -e EcoRI,BamHI,BsaI
```

```
=== Digest of NM_018947.6 Homo sapiens cytochrome c, somatic (CYCS), mRNA ===
Sequence length: 5432 bp, linear
Enzymes: EcoRI, BamHI, BsaI
  EcoRI (G^AATTC) cut sites at: [475, 1053, 4002, 4177]
  BamHI (G^GATCC): no cut sites found
  BsaI (GGTCTC(1/5)) cut sites at: [178, 1738, 2358, 3777, 4587]
    3 site(s) on the top strand, 2 on the bottom strand

--- Fragment Summary ---
  # Size (bp)
--------------
  1      1419
  2       845
  3       685
  4       620
  5       578
  6       410
  7       297
  8       225
  9       178
 10       175
--------------
Total      5432
```

BsaI is non-palindromic, so it is searched on both strands and the per-strand
counts are reported; EcoRI is symmetric, so each of its sites is reported once.

### Options

| Flag | Meaning |
|---|---|
| `-f`, `--fasta` | Input FASTA (required) |
| `-e`, `--enzymes` | Comma separated enzyme names or custom specs (required) |
| `-c`, `--circular` | Treat the DNA as circular; default is linear |
| `-m`, `--min-fragment` | Smallest fragment to report, in bp (default 1) |
| `-o`, `--output` | `table` (default), `gel`, or `both` |
| `-a`, `--ambiguity` | `definite` or `possible` — also opens the input alphabet, see below |
| `--gel-height` | Height of the ASCII gel (default 30) |
| `--enzyme-db` | Read enzymes from a CSV instead of the built in table |
| `--list-enzymes` | List the built in enzymes and exit |
| `--convert` | Show one enzyme in both notations and exit |

### Enzyme specifications

Built in enzymes go by name, case insensitively (`EcoRI`, `bsai`, `HincII`).
`--list-enzymes` prints all 262 with their sites and end types. Custom enzymes
can be written four ways:

| Form | Example | Meaning |
|---|---|---|
| `NAME:RECOGNITION` | `PmeI:GTTT^AAAC` | NEB caret notation, cut marked inside the site |
| `NAME:RECOGNITION` | `MyI:GGTCTC(1/5)` | NEB offset notation, cuts 1 nt (top) and 5 nt (bottom) past the site |
| `NAME:SEQ:OFFSET` | `PmeI:GTTTAAAC:4` | Top strand cut; the bottom strand is assumed to mirror it |
| `NAME:SEQ:TOP:BOTTOM` | `MyI:GGTCTC:7:11` | Both cut offsets, from the start of the site |

### Ambiguity

Without `--ambiguity` the input must be plain `ACGT` and anything else is a hard
error. That is deliberate: an ambiguous base that silently fails to match is how
a sequence gets wrongly reported as uncut.

Passing the flag opens the alphabet to every IUPAC code plus `-` gaps, and its
value says how to read them:

- `definite` — count a site only if it is cut *however* the ambiguous bases
  resolve. "Where will this definitely cut?"
- `possible` — count a site if it is cut under *at least one* resolution.
  "Could this cut at all?", which is the question behind trusting that an
  enzyme leaves a sequence intact. Speculative cuts are marked `?`.

Every definite site is also a possible site, never the other way round. Gaps are
stripped before digestion — a gap is not a base, the molecule reads through it,
so `AAT-ATT` really is an SspI site — and cut positions are reported in ungapped
coordinates.

### Cut coordinates

Offsets are 0-based from the first base of the recognition site, both measured
along the top strand, with the cut falling immediately before that offset.
`cut_bottom - cut_top` is the overhang: positive 5', zero blunt, negative 3'.
Offsets may fall outside the site — type IIS enzymes cut downstream, `^GATC`
(MboI) cuts before. A site found on the bottom strand mirrors: an offset `p` in
the enzyme's frame lands at `start + len(site) - p`, so a reverse orientation hit
cuts *upstream* of where it was found.

On linear DNA a cut falling past either end is reported as recognised-but-uncut.
On circular DNA it wraps, and the scan window is extended so sites straddling the
origin are found too.

### Enzyme data

`Data/restriction_enzymes.csv` is the single source of truth for the enzyme
table — 262 usable enzymes, plus 13 nicking and 12 type IIB enzymes that are
quarantined rather than dropped, so asking for one gives a named error instead of
"unknown enzyme". It uses NEB bench notation, `^` for the top strand cut and `_`
for the bottom:

```
EcoRI,G^AATTC        AatII,G_ACGT^C        BsaI,GGTCTCN^NNNN_
```

CMake embeds that file verbatim into a generated header at configure time, and
the built in table is parsed from it by the same loader `--enzyme-db` uses on a
real file. So there is exactly one copy of the data and one parser: **edit the
CSV, never a table in the source.** Reconfigure (`cmake .`) to pick up an edit.

Specificities come from REBASE, never from recall — wrong cut data is silent and
biological, the worst failure mode here.

### Known limitations

Carried over from the Python original unchanged:

- Enzymes that cut on both sides of their site (BaeI, BcgI, CspCI) are rejected
  rather than partly parsed — two cuts per site do not fit the fragment model.
- Fragment lengths are top strand lengths, so overhangs are not reflected in the
  reported sizes.
- Methylation sensitivity, star activity, and enzymes needing two sites are not
  modelled.

## Test files

`Data/test_files/` holds known truth fixtures — the expected consensus is obvious by
construction, so the tool can be checked against it. Each `TEST_IUPAC_*` file
contains one row of each base making up that ambiguity code (per the
[IUPAC table](https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation)),
so the whole consensus should be that one letter repeated.

| File | Contents | Expected with `-a` | Expected strict |
|---|---|---|---|
| `TEST_IUPAC_R.msa.fna` | A, G | `R` | `N` |
| `TEST_IUPAC_Y.msa.fna` | C, T | `Y` | `N` |
| `TEST_IUPAC_S.msa.fna` | G, C | `S` | `N` |
| `TEST_IUPAC_W.msa.fna` | A, T | `W` | `N` |
| `TEST_IUPAC_K.msa.fna` | G, T | `K` | `N` |
| `TEST_IUPAC_M.msa.fna` | A, C | `M` | `N` |
| `TEST_IUPAC_B.msa.fna` | C, G, T | `B` | `N` |
| `TEST_IUPAC_D.msa.fna` | A, G, T | `D` | `N` |
| `TEST_IUPAC_H.msa.fna` | A, C, T | `H` | `N` |
| `TEST_IUPAC_V.msa.fna` | A, C, G | `V` | `N` |
| `TEST_IUPAC_N.msa.fna` | A, C, G, T | `N` | `N` |
| `TEST_ENCODE_AMBIG_INPUT.msa.fna` | `R`, C — ambiguity code already in the input | `V` | `N` |
| `TEST_N_PASSTHROUGH_MID.msa.fna` | A, G, `N` — literal N mid alignment | `N` | `N` |
| `TEST_N_PASSTHROUGH_SEED.msa.fna` | `N`, A, G — literal N in the first row | `N` | `N` |
| `TEST_GAP_MIDDLE.msa.fna` | A, `-`, G — gap mid alignment | `-` | `-` |

The other fixtures are real sequence: `TEST.fna` (single CYCS mRNA),
`TEST_IDNT.fna` (that same sequence four times over, so the consensus should
come back as the input sequence unchanged — output is rewrapped at 80 columns),
`Sheep_products.fna` (1000 sequences, long enough that the accumulator has to
compress), `TEST.faa` / `NP_061820.1` (protein).

These fixtures exercise the tool end to end through its CLI. The finer grained
checks on the individual functions live in the CUnit suite, see
[Running the tests](#running-the-tests).

Check everything at once:

```bash
for code in R Y S W K M B D H V N; do
  printf '%s -> %s\n' "$code" \
    "$(./bin/MSA-to-consensus -i Data/test_files/TEST_IUPAC_$code.msa.fna -a -o stdout 2>/dev/null | tail -n1)"
done
```

## License

**This project is licensed under the GPL v3 — see [LICENSE](LICENSE).** That is
the licence you receive it under and the one that binds anything you build on
it. There is no second licence to choose from and no part of the repository is
offered under different terms.

Copyright © 2026 The Zabel Lab at Colorado State University, for the work in
this repository — including the C++ port described below, which is a substantial
body of original work in its own right even though it is a derivative of the
code it was ported from.

### Third party attribution

`CppSrc/Enzyme-digest/Restriction_Enzyme_Digest.cpp` and
`Data/restriction_enzymes.csv`
are derived from the Python restriction digest simulator, © Collins Amatu
Gorgerat, which is MIT licensed —
[LICENSE.restriction-digest](LICENSE.restriction-digest).

This is attribution, not an alternative licence. MIT is GPL compatible, which is
exactly why that code can be absorbed here: the combined work is distributed
under the GPL v3 above, and the MIT notice travels with the portion it covers.
Keeping the notice is a condition of the permission that lets us use the code at
all, so it stays with those files, with the source header, and with anything
built or installed from them — regardless of how far the port diverges from the
original. A translation into another language is a derivative work; the extent
of the changes does not retire the obligation.
