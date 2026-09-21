# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository.

## Project status

Three C++17 command line tools for the Zabel Lab at Colorado State University,
built with CMake:

- `MSA-to-consensus` — collapses a multiple sequence alignment into one
  consensus sequence.
- `Enzyme-digest` — in silico restriction digest, ported from Python in
  August 2026 and verified byte-for-byte against the original.
- `PCR-protocol` — reaction setup and thermocycler program for a primer pair
  and polymerase, with primer Tm and primer3 `thal` hairpin/dimer checks.

`README.md` is the user-facing documentation and is kept accurate; prefer
reading it over re-deriving how a tool behaves.

## Commands

```bash
cmake .                                   # configure (also embeds the enzyme CSV)
cmake --build .                           # all three tools + test binaries
cmake --build . --target Enzyme-digest    # one target
ctest                                     # all suites, non-zero on failure
./bin/test_digest                         # one suite, verbose per-test output
cmake --build . --target clean            # see "No .gitignore" below

./bin/MSA-to-consensus -i Data/test_files/TEST_IUPAC_R.msa.fna -a
./bin/Enzyme-digest -f Data/test_files/TEST.fna -e EcoRI,BamHI,BsaI
./bin/Enzyme-digest --list-enzymes        # 262 enzymes + 25 quarantined
./bin/PCR-protocol -p Q5 -f Data/test_files/TEST_PCR_LAMBDA.fna \
    -F GTCACCAGTGCAGTGCTTGATAACAGG -R GATGACGCATCCTCACGATAATATCCGG
```

**Binaries built with CLion's bundled MinGW need its runtime on PATH or they
die with a bare segfault and no diagnostic — even a hello-world does.** If a
freshly built binary segfaults instantly, this is why, not your change:

```bash
export PATH="/c/Program Files/JetBrains/CLion 2026.2.1/bin/mingw/bin:$PATH"
```

## Architecture

**Each tool is a single translation unit that owns its own `main()`.** The
CUnit tests in `Tests/` `#define main disabled_main` and `#include` the `.cpp`
directly, so they call the real functions with no seam and the program source
needs no changes to be testable. Keep it that way: splitting a tool into
`.hpp`/`.cpp` breaks that pattern.

**The one exception is header-only code shared between tools, in `CppSrc/common/`.**
`fasta.hpp` is the single FASTA reader for all three; it is not a target and has
no `main()`, so the rule above is untouched. Each tool's `parse_fasta` is a thin
adapter over it: it maps the reader's typed errors (`fasta::problem`) onto the
tool's own messages and adds the tool's own validation through
`options::check_line`. That split is what let `Enzyme-digest` keep its
Python-parity messages, and it is also what keeps the MIT lineage out of the
shared file, so do not move that wording into `fasta.hpp`. `Tests/test_fasta.cpp`
includes the header directly. Adding a file here means adding a `test_*` for it.

**Directory name = CMake target name = binary name.** `CppSrc/Enzyme-digest/`
builds target `Enzyme-digest` into `bin/Enzyme-digest`. `CppSrc/` holds source
only; data lives in `Data/`.

**`Data/restriction_enzymes.csv` is the single source of truth for the enzyme
table.** CMake reads it at configure time and embeds it verbatim into a
generated header, and the built-in table is parsed from that string by the same
loader `--enzyme-db` uses on a real file — one copy of the data, one parser.
Never transcribe enzymes into a C++ table; two copies would drift. Re-run
`cmake .` after editing the CSV.

`Data/polymerases.csv` follows the same pattern for `PCR-protocol`, and the
same rule: **every polymerase number comes from the manufacturer's datasheet**,
cited in its `source` column. Anything not straight off the page goes in the
`notes` column, which the tool prints with every protocol.

**`PCR-protocol` computes Tm itself but takes hairpins and dimers from primer3.**
Our SantaLucia 1998 Tm matches primer3's `oligotm()` to 1e-6 C, and the tests
call both side by side. Structures come from `thal`, built from the
`External_tools/primer3` submodule (v2.6.1, parameters compiled in) as the
static library `primer3_thal`. A home-grown perfect-match dimer model was tried
first and missed ~80% of the self-dimers `thal` flags; do not reintroduce one.
`thal.h` has no `extern "C"` guard, so it is included inside one.

**Enzyme specificities come from REBASE, never from recall.** Wrong cut data is
silent and biological — the worst failure mode here.

**Cut coordinates are the thing to get right.** Both offsets are 0-based from
the first base of the recognition site, both measured along the *top* strand,
the cut falling immediately before that offset. `cut_bottom - cut_top` is the
overhang: positive 5', zero blunt, negative 3'. Offsets may fall outside the
site — type IIS cuts downstream, `^GATC` (MboI) cuts before. A bottom-strand
hit mirrors: enzyme-frame offset `p` lands at `start + len(site) - p`, so a
reverse-orientation hit cuts *upstream* of where it was found.

**Matching is bitmask, not regex.** Each IUPAC code is a 4-bit mask over
{A,C,G,T}. `definite` requires the target's base set to be a subset of the
site's; `possible` requires only an intersection. The direction is asymmetric
and trivially invertible, so it is pinned by a test. A non-IUPAC character has
mask 0 and must match under neither — otherwise it would be a subset of
everything.

## Conventions that matter here

**No `.gitignore` in the repo, deliberately.** The build happens in-source and
`cmake --build . --target clean` removes everything configuring generated, so
`git status` comes back clean without one. Consequence: **stage by explicit
path, never `git add -A`** — an earlier commit swept ~33k lines of `CMakeFiles/`
into the repo that way. `.idea/` is untracked and should stay that way.

A local `.gitignore` may sit in the working tree; it lists itself on its first
line, so it never shows in `git status` and is never committed. Leave it
untracked. It also hides `cmake-build-debug/` internals and `.cmake/`, so do not
trust a clean `git status` on its own as proof that `clean` worked — use
`git status --short --ignored`.

**Do not derive test expectations from the implementation.** Several tests in
`Tests/test_digest.cpp` exist specifically to be independent of the code they
check, and re-deriving them with the same comprehension would make them test
themselves:

- the hand-written 15-code definite-matching table,
- fragments summing to sequence length across all 262 enzymes, linear and
  circular — this is what catches sign errors and off-by-ones in bulk,
- overhang mirroring for palindromic in-site cutters, which is forced by the
  symmetry of the molecule and so is independent of REBASE too.

**Verify a port or a matching change differentially.** The Python original is
no longer vendored; clone it separately at commit `82994c5` of
`github.com/Nheyer/restriction-enzyme-digest-simulator` and run both over the
same corpus. The CLI surfaces match apart from the program name, which differs
in the usage line and the required-argument error. Include `B/D/H/V` in any
ambiguity corpus — that is where subset and intersection diverge most. For
anything touching input handling add CRLF, blank-only, empty, headerless and
invalid-character files: after the shared reader went in, 140 such runs matched
the pre-change binary byte for byte, and matched the Python on all but the
control-character message below.

**Licensing: GPL v3 governs, and the MIT notice stays.** `Enzyme-digest` is a
derivative of MIT-licensed work by Collins Amatu Gorgerat. MIT is
GPL-compatible, so the combined work ships under GPL v3, but keeping the notice
is a condition of the permission that allows the code's use at all. Do not
remove the source header block or `LICENSE.restriction-digest`, and keep both in
the CMake install list. The extent of divergence from the original does not
retire the obligation — a translation is a derivative work. This was raised and
settled; do not relitigate it.

**Scope of the MIT notice** (audited September 2026; the README has the table).
It bites on the port `.cpp` — `digest_linear`, `digest_circular`,
`print_fragment_table`, `draw_ascii_gel`, the shape of `main()`, and the wording
of `parse_fasta`'s errors — and on the 100 bp ladder in `enzyme_data.hpp`. The
FASTA *reading* is not in that list: it is `CppSrc/common/fasta.hpp`, which grew
out of the `MSA-to-consensus` reader and takes nothing from the Python. `Data/restriction_enzymes.csv` is **not**
Collins's: it is the fork's NEB compilation, every line by Nheyer, and Collins's
original was a 20-enzyme table in another format. Everything else is plain
GPL v3, tagged `SPDX-License-Identifier: GPL-3.0-only`. The CSVs cannot carry a
header (the loader would read it as data), so the README covers them.
When re-checking against the fork, blame with `git blame -C -C -C -M 82994c5`.
Plain blame credits `enzyme_data.py` wholly to whoever *moved* the tables and
shows no Collins line in it, which hides the ladder.

## Known warts

- **`Enzyme-digest`'s `py_repr` prints a control character raw**, where Python
  prints its escape (`'\x0c'` for a form feed). Both tools refuse the input; only
  the invalid-character message differs, so "byte for byte" is true for every
  input except one holding a control character. Found September 2026 with a
  form feed and a vertical tab in a sequence line. Unfixed, and unrelated to the
  shared FASTA reader.
- **`encode_nucliotide` has no return for non-IUPAC input**
  (`CppSrc/MSA-to-consensus/`). It warns `control reaches end of non-void
  function` on every build and is UB if ever reached. Unfixed.
- **"nucliotide" is misspelled in ~67 places**, including three public function
  names and the tests. **This was offered and explicitly declined** — do not
  mass-rename it unprompted.
- **`write_fasta` only catches a failed *open***, not a write that fails partway
  (full disk, broken pipe); nothing checks stream state after the loop.
- **`MSA-to-consensus` return codes are negative** (`rt_READ_FAILED` = -1 etc.),
  which the shell cannot tell apart — they all surface as one number. Build with
  `DEBUG > 1` to have `describe_return()` name the failure on stderr.
- **`Enzyme-digest` scans once per enzyme**, so a full-database digest passes
  over the sequence 262 times. A multi-pattern single pass is the remaining
  algorithmic win; the port already beat the Python by ~23-62x on constant
  factor alone.
