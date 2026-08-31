# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository.

## Project status

Two C++17 command line tools for the Zabel Lab at Colorado State University,
built with CMake:

- `MSA-to-consensus` — collapses a multiple sequence alignment into one
  consensus sequence.
- `Enzyme-digest` — in silico restriction digest, ported from Python in
  August 2026 and verified byte-for-byte against the original.

`README.md` is the user-facing documentation and is kept accurate; prefer
reading it over re-deriving how a tool behaves.

## Commands

```bash
cmake .                                   # configure (also embeds the enzyme CSV)
cmake --build .                           # both tools + both test binaries
cmake --build . --target Enzyme-digest    # one target
ctest                                     # both suites, non-zero on failure
./bin/test_digest                         # one suite, verbose per-test output
cmake --build . --target clean            # see "No .gitignore" below

./bin/MSA-to-consensus -i Data/test_files/TEST_IUPAC_R.msa.fna -a
./bin/Enzyme-digest -f Data/test_files/TEST.fna -e EcoRI,BamHI,BsaI
./bin/Enzyme-digest --list-enzymes        # 262 enzymes + 25 quarantined
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

**Directory name = CMake target name = binary name.** `CppSrc/Enzyme-digest/`
builds target `Enzyme-digest` into `bin/Enzyme-digest`. `CppSrc/` holds source
only; data lives in `Data/`.

**`Data/restriction_enzymes.csv` is the single source of truth for the enzyme
table.** CMake reads it at configure time and embeds it verbatim into a
generated header, and the built-in table is parsed from that string by the same
loader `--enzyme-db` uses on a real file — one copy of the data, one parser.
Never transcribe enzymes into a C++ table; two copies would drift. Re-run
`cmake .` after editing the CSV.

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

**No `.gitignore`, deliberately.** The build happens in-source and
`cmake --build . --target clean` removes everything configuring generated, so
`git status` comes back clean without one. Consequence: **stage by explicit
path, never `git add -A`** — an earlier commit swept ~33k lines of `CMakeFiles/`
into the repo that way. `.idea/` is untracked and should stay that way.

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
ambiguity corpus — that is where subset and intersection diverge most.

**Licensing: GPL v3 governs, and the MIT notice stays.** `Enzyme-digest` is a
derivative of MIT-licensed work by Collins Amatu Gorgerat. MIT is
GPL-compatible, so the combined work ships under GPL v3, but keeping the notice
is a condition of the permission that allows the code's use at all. Do not
remove the source header block or `LICENSE.restriction-digest`, and keep both in
the CMake install list. The extent of divergence from the original does not
retire the obligation — a translation is a derivative work. This was raised and
settled; do not relitigate it.

## Known warts

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
