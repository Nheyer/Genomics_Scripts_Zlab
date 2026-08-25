# Genomics_Scripts_Zlab

Small genomics command line tools used in the Z lab, plus a few external tools
vendored as submodules.

## Contents

| Path | What it is |
|---|---|
| `Scripts/MSA_fasta_to_Consensus.cpp` | `MSAtoCONSENSUS` — collapses a multiple sequence alignment into a single consensus sequence |
| `Tests/test_consensus.cpp` | CUnit unit tests, see [Running the tests](#running-the-tests) |
| `External_tools/argparse` | [p-ranav/argparse](https://github.com/p-ranav/argparse), header only CLI parsing (build dependency) |
| `External_tools/cunit` | [cunity/cunit](https://gitlab.com/cunity/cunit), unit test framework (build dependency) |
| `External_tools/ReDtool` | [CBL205NIPGR/ReDtool](https://github.com/CBL205NIPGR/ReDtool), restriction digest tool |
| `External_tools/restriction-enzyme-digest-simulator` | [wl5e/restriction-enzyme-digest-simulator](https://github.com/wl5e/restriction-enzyme-digest-simulator) |
| `test_files/` | Known truth fixtures, see [Test files](#test-files) |

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
cmake --build . --target MSAtoCONSENSUS
```

The binary lands in `bin/MSAtoCONSENSUS`. The build uses the vendored
`External_tools/argparse` headers, not any system wide install, so it does not
matter whether argparse is installed on the machine.

### Running the tests

`Tests/test_consensus.cpp` unit tests the consensus code with CUnit:

```bash
cmake --build . --target test_consensus
./bin/test_consensus
```

or through CTest, which is what CI would use:

```bash
ctest
```

Either way the runner exits non zero if anything fails. The suite covers
`clean_nucliotide`, `encode_nucliotide`, `decode_nucliotide` and
`make_consensus` — ambiguity codes, strict masking, gap handling, and the
accumulator compression path.

Because the tests build the alignments in memory and call `make_consensus`
directly, they can cover things a FASTA fixture makes awkward: an alignment
long enough to force compression is `repeated("G", 1000)` rather than a
thousand line file.

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

## MSAtoCONSENSUS

Reads a multiple sequence aligned FASTA file and writes a single consensus
sequence. Every input sequence must be the same length (it is an *aligned*
file); mismatched lengths are rejected.

```bash
./bin/MSAtoCONSENSUS -i test_files/TEST_IUPAC_R.msa.fna -a
```

```
Opening: test_files/TEST_IUPAC_R.msa.fna as multiply aligned Nucliotide Sequence input!
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

## Test files

`test_files/` holds known truth fixtures — the expected consensus is obvious by
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
    "$(./bin/MSAtoCONSENSUS -i test_files/TEST_IUPAC_$code.msa.fna -a -o stdout 2>/dev/null | tail -n1)"
done
```

## License

GPL v3, see [LICENSE](LICENSE).
