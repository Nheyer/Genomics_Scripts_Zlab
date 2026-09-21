# Genomics_Scripts_Zlab

Small genomics command line tools used in the Zabel Lab at Colorado State
University.

## Install

Requires CMake >= 3.15 and a C++17 compiler (plus a C compiler, for the
vendored primer3). The three build dependencies are vendored as submodules, so
clone recursively:

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

The binaries land in `bin/MSA-to-consensus`, `bin/Enzyme-digest` and
`bin/PCR-protocol`; build one on its own by naming it as the `--target`. The build uses the vendored
`External_tools/argparse` headers, not any system wide install, so it does not
matter whether argparse is installed on the machine.

Configuring also embeds `Data/restriction_enzymes.csv` and
`Data/polymerases.csv` into generated headers under `generated/`, which is why
neither `Enzyme-digest` nor `PCR-protocol` needs a data file at runtime. See
[Enzyme data](#enzyme-data) and [Polymerase data](#polymerase-data).

To put the tools on your PATH:

```bash
cmake --install .
```

(`make install` does the same thing.) That drops the binaries and the docs:

```
<prefix>/bin/MSA-to-consensus
<prefix>/bin/Enzyme-digest
<prefix>/bin/PCR-protocol
<prefix>/share/doc/Genomics_Scripts_Zlab/
```

The prefix defaults to `/usr/local`, which needs `sudo`. To install somewhere
you own instead, no `sudo` required:

```bash
cmake --install . --prefix ~/.local
```

The installed docs include `LICENSE.restriction-digest` and `LICENSE.primer3`.
Those notices have to travel with the binaries, so do not drop them from the
install list — see [License](#license). The test binaries and the vendored CUnit and primer3
libraries are deliberately not installed; primer3 is linked statically into
`PCR-protocol`, and CUnit is only there to check the tools.

There is no `uninstall` target, but installing writes `install_manifest.txt`
listing every file it placed, so removing them is:

```bash
xargs rm -f < install_manifest.txt
```

## The tools

| Tool | What it does |
|---|---|
| [`MSA-to-consensus`](#msa-to-consensus) | Collapses a multiple sequence alignment into one consensus sequence, masking disagreements or coding them as IUPAC ambiguity |
| [`Enzyme-digest`](#enzyme-digest) | In silico restriction digest — cut sites, fragment sizes and an ASCII gel, over a 262 enzyme database |
| [`PCR-protocol`](#pcr-protocol) | A reaction setup and thermocycler program for a primer pair and polymerase, with primer Tm and hairpin/dimer checks |

All three write their result to stdout (or, for the first two, to a file), with
progress and errors on stderr, so redirecting stdout captures just the output.

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

## PCR-protocol

Takes two primers and a polymerase and writes out the reaction: primer Tms,
hairpin and dimer checks, a reaction setup, and a thermocycler program with the
annealing temperature and extension time filled in.

```bash
./bin/PCR-protocol -p Q5 -f Data/test_files/TEST_PCR_LAMBDA.fna \
    -F GTCACCAGTGCAGTGCTTGATAACAGG -R GATGACGCATCCTCACGATAATATCCGG
```

```
PCR protocol: Q5 High-Fidelity DNA Polymerase (NEB M0491)

Primers
  Forward 5'-GTCACCAGTGCAGTGCTTGATAACAGG-3'
          27 nt, GC 51.9%, Tm 69.0 C
  Reverse 5'-GATGACGCATCCTCACGATAATATCCGG-3'
          28 nt, GC 50.0%, Tm 68.0 C
  Amplicon 1347 bp (NC_001416.1:29951-40100, bases 56..1402)

Primer structures (primer3 thal, Tm of the most stable; flagged above 47 C)
                     3' end      anywhere
  Forward hairpin                35.4 C
  Forward self-dimer 18.5 C      25.1 C
  Reverse hairpin                35.5 C
  Reverse self-dimer none        none
  Cross-dimer        none        none

Reaction setup (50 uL)
  Component                           Volume      Final
  5X Q5 Reaction Buffer               10 uL       1X
  10 mM dNTPs                         1 uL        200 uM each
  10 uM Forward primer                2.5 uL      0.5 uM
  10 uM Reverse primer                2.5 uL      0.5 uM
  Template DNA                        variable
  Q5 polymerase (2 U/uL)              0.5 uL      1 U
  Nuclease-free water                 to 50 uL

Thermocycler program
  Initial denaturation    98 C    30 s
  30 cycles of:
    Denaturation          98 C    10 s
    Annealing             71 C    30 s
    Extension             72 C    41 s
  Final extension         72 C    2 min
  Hold                    4 C     forever

Notes
  Annealing: lower Tm 68.0 C + 3 C per the Q5 datasheet. A gradient around it is the surest check.
  ...

Warnings
  ! Reverse primer has 4 G/C in its last 5 bases (more than 3 invites mispriming)
```

These are NEB's own lambda control primers, and the 1347 bp product lands on
the lambda coordinates NEB gives for them (30,006–31,352).

### Options

| Flag | Meaning |
|---|---|
| `-F`, `--forward` | Forward primer, 5'→3' (required, unless `--primers`) |
| `-R`, `--reverse` | Reverse primer, 5'→3' **as ordered** — not the top strand sequence (required, unless `--primers`) |
| `-P`, `--primers` | Both primers from a FASTA instead of `-F`/`-R` — see below |
| `-p`, `--polymerase` | `Q5`, `Taq` or `Phusion`, case insensitively (required) |
| `-f`, `--fasta` | Template FASTA; both primers are located on it to get the amplicon length |
| `-l`, `--amplicon-length` | Amplicon length in bp, instead of `--fasta` |
| `-c`, `--cycles` | Number of cycles (default: the polymerase's, 30 for all three) |
| `--volume` | Reaction volume in µL (default 50); the setup scales from the datasheet's 50 µL |
| `-n`, `--replicates` | Number of reactions (default 1); above 1, a master mix is added — see below |
| `--template-volume` | Template per tube in µL (default 1), left out of the master mix |
| `--simple-template` | Plasmid, lambda or *E. coli* template: use the datasheet's faster extension rate |
| `--polymerase-db` | Read polymerases from a CSV instead of the built in table |
| `--list-polymerases` | List the known polymerases and exit |

`--primers` takes a FASTA holding exactly two records: the forward primer first,
the reverse second, both 5'→3' as ordered. Order is the only thing that says
which is which, so any other count is an error. Names are free, sequences may be
lower case or wrapped, and it cannot be combined with `-F`/`-R`:

```bash
./bin/PCR-protocol -p Q5 -P Data/test_files/TEST_PCR_PRIMERS.fna \
    -f Data/test_files/TEST_PCR_LAMBDA.fna
```

Exactly one of `--fasta` and `--amplicon-length` is needed — extension time
depends on the product length, and two primers alone do not give it.

With `--fasta`, the primers must match the template exactly. The template may
be either strand and may hold several records; the pair has to make exactly one
product across all of them, or the tool lists what it found and stops rather
than guess. If nothing is found it says why — most usefully when the reverse
primer was given as the top strand sequence instead of 5'→3' as ordered.

### Master mix

With `--replicates` above 1 the output adds a master mix for that many
reactions plus 10% for pipetting loss. The template goes into each tube on its
own, so the mix holds everything else, with the water worked out from
`--template-volume`: aliquot the mix, then add the template.

```bash
./bin/PCR-protocol -p Q5 -P Data/test_files/TEST_PCR_PRIMERS.fna -l 1347 -n 8
```

```
Master mix (8 reactions + 10% = 8.8, 431.2 uL)
  Component                           Per rxn     Mix
  5X Q5 Reaction Buffer               10 uL       88 uL
  10 mM dNTPs                         1 uL        8.8 uL
  10 uM Forward primer                2.5 uL      22 uL
  10 uM Reverse primer                2.5 uL      22 uL
  Q5 polymerase (2 U/uL)              0.5 uL      4.4 uL
  Nuclease-free water                 32.5 uL     286 uL
  Aliquot 49 uL per tube, then add 1 uL template.
```

A template too big to fit beside the reagents is an error rather than a
negative water volume.

### How the numbers are worked out

- **Primer Tm** — the SantaLucia 1998 unified nearest-neighbour model with its
  salt correction, magnesium counted as sodium equivalents
  (`[Na+] + 120·√([Mg²+] − [dNTP])`), at the polymerase's buffer and primer
  concentrations. These are primer3's default formulas, and the tests check the
  two agree to 10⁻⁶ °C.
- **Annealing** — the datasheet rule applied to the lower Tm, rounded to a whole
  degree: Q5 is Tm + 3; Phusion is Tm + 3 when both primers are over 20 nt and
  Tm otherwise; Taq is Tm − 5, which is where NEB says to start a gradient (the
  Taq datasheet gives no fixed rule). If that reaches the extension temperature
  the program goes 2-step, annealing folded into extension.
- **Extension** — amplicon length times the datasheet rate, rounded up to a
  whole second. Q5 and Phusion use their rate for complex (genomic) templates
  unless `--simple-template`; Q5 goes to 50 s/kb over 6 kb; Taq never drops
  below 45 s.
- **Primer structures** — hairpins, self-dimers and the cross-dimer, from
  primer3's thermodynamic aligner (`thal`), which handles mismatches, bulges,
  loops and dangling ends. Each is the melting temperature of the most stable
  structure, and anything above 47 °C is flagged — primer3's own default limit.
  "3' end" is the structure anchored on a primer's 3' end, the one the
  polymerase can extend.
- **Other checks** — GC 40–60% and length 20–40 nt (the NEB datasheets); the
  two Tms within 5 °C, at most 3 G/C in the last 5 bases, no runs over 4 and no
  more than 4 dinucleotide repeats (Premier Biosoft's guidelines).

### Polymerase data

`Data/polymerases.csv` is the single source of truth for the polymerases, in
the same arrangement as the enzyme table: embedded at configure time and parsed
by the loader `--polymerase-db` uses. **Every number in it comes from the
manufacturer's datasheet**, cited in the `source` column:

| Polymerase | Source |
|---|---|
| Q5 (NEB M0491) | NEB protocol *PCR Using Q5 High-Fidelity DNA Polymerase (M0491)* |
| Taq (NEB M0273) | NEB *PCR with Taq DNA Polymerase (M0273)*, [protocols.io](https://dx.doi.org/10.17504/protocols.io.ch7t9m); buffer composition and the 45 s floor from the Taq PCR Kit (E5000) manual |
| Phusion (NEB M0530) | NEB *Phusion High Fidelity PCR Kit (E0553)* manual, v4.0 |

Where a datasheet gives a range, the table takes the longer or more
conservative end — 10 s denaturation for the high fidelity enzymes, the
genomic extension rate, the longer annealing and final extension times — and
the `notes` column records anything that is not straight off the page. Those
notes are printed with every protocol. Add a polymerase by adding a row, from
its datasheet; the loader refuses a row with a missing or non-numeric column
rather than let it become a zero second step.

### Known limitations

- **Tm is not NEB's Tm.** The datasheets tune their annealing rules against
  NEB's Tm calculator, whose model is not published, and it gives different
  numbers. NEB's Phusion manual quotes Tms for its own control primers; under
  the Phusion buffer this tool gives:

  | Primer | NEB | This tool |
  |---|---|---|
  | `GTCACCAGTGCAGTGCTTGATAACAGG` | 71.0 | 68.0 |
  | `GATGACGCATCCTCACGATAATATCCGG` | 73.2 | 67.0 |
  | `CAGTGCAGTGCTTGATAACAGG` | 63.0 | 62.4 |
  | `GTAGTGCGCGTTTGATTTCC` | 62.7 | 60.9 |

  So the annealing temperature can come out several degrees under what NEB's
  calculator would suggest, and that can change the program's shape, not just
  a number: NEB runs the 1.3 kb pair as a 2-step program at 72 C, where this
  tool, from its lower Tms, anneals separately at 70 C. Run a gradient around
  it.
- **Monovalent salt for Q5 and Phusion is assumed.** NEB does not publish those
  buffer compositions, so 50 mM is used, and the output says so. Taq's 50 mM
  KCl is published.
- **Primers must match the template exactly.** Tailed primers (restriction
  sites, Gibson overlaps) and degenerate primers are not located or scored; for
  a tailed pair, give `--amplicon-length` with the tails included.
- **Linear templates only.** A product spanning the origin of a circular
  template is not found.

## Repository layout

| Path | What it is |
|---|---|
| `CppSrc/` | C++ source only, one directory per program |
| `CppSrc/MSA-to-consensus/MSA_fasta_to_Consensus.cpp` | `MSA-to-consensus` |
| `CppSrc/Enzyme-digest/Restriction_Enzyme_Digest.cpp` | `Enzyme-digest` |
| `CppSrc/Enzyme-digest/enzyme_data.hpp` | IUPAC codes, complement table and ladder — lookup tables only, no logic |
| `CppSrc/Enzyme-digest/enzyme_csv_embedded.hpp.in` | Template CMake fills with the enzyme CSV at configure time; required to build |
| `CppSrc/PCR-protocol/PCR_Protocol.cpp` | `PCR-protocol` |
| `CppSrc/PCR-protocol/pcr_data.hpp` | Nearest-neighbour parameters and primer limits, each with its source — lookup tables only, no logic |
| `CppSrc/PCR-protocol/polymerase_csv_embedded.hpp.in` | Template CMake fills with the polymerase CSV at configure time; required to build |
| `CppSrc/common/fasta.hpp` | The FASTA reader all three tools share. Header only, not a target, no `main()` — each tool is still one translation unit |
| `Data/restriction_enzymes.csv` | The enzyme database, in NEB `^`/`_` notation. **The source of truth** — see [Enzyme data](#enzyme-data) |
| `Data/polymerases.csv` | The polymerase table, from the datasheets. **The source of truth** — see [Polymerase data](#polymerase-data) |
| `Data/test_files/` | Known truth fixtures, see [Test files](#test-files) |
| `Tests/test_consensus.cpp` | CUnit unit tests, see [Running the tests](#running-the-tests) |
| `Tests/test_digest.cpp` | CUnit unit tests for `Enzyme-digest`, same |
| `Tests/test_pcr.cpp` | CUnit unit tests for `PCR-protocol`, same |
| `Tests/test_fasta.cpp` | CUnit unit tests for the shared reader; it includes the header directly, there is no `main()` to rename |
| `LICENSE.restriction-digest` | MIT notice for the code `Enzyme-digest` was ported from, see [License](#license) |
| `External_tools/argparse` | [p-ranav/argparse](https://github.com/p-ranav/argparse), header only CLI parsing (build dependency) |
| `External_tools/cunit` | [cunity/cunit](https://gitlab.com/cunity/cunit), unit test framework (build dependency) |
| `External_tools/primer3` | [primer3-org/primer3](https://github.com/primer3-org/primer3) v2.6.1; `PCR-protocol` builds its `thal` hairpin/dimer aligner (build dependency) |

## Development

### Running the tests

`Tests/` unit tests all three tools with CUnit:

```bash
cmake --build . --target test_consensus test_digest test_pcr test_fasta
./bin/test_consensus
./bin/test_digest
./bin/test_pcr
./bin/test_fasta
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

The PCR suite pins every expected number to something outside the code:

- **The paper against itself.** SantaLucia 1998 gives the nearest-neighbour
  parameters twice, as ΔH/ΔS (Table 2, what we compute from) and as ΔG₃₇
  (Table 1, typed in separately in the test). Each row has to reproduce the
  other, so a transcription slip in either fails.
- **Tm against primer3**, both by calling the vendored `oligotm()` side by side
  and against values from a separate primer3-py build, for all three buffers.
- **thal against the primer3 manual's worked examples**, which it reproduces to
  four decimals.
- **NEB's lambda control primers**, located at the genome coordinates NEB
  states, including the reverse primer's orientation — given as the top strand
  sequence it must find nothing.
- **Annealing and extension by hand**, including every boundary the datasheets
  draw: 20 vs 21 nt for Phusion, "over 6 kb" for Q5, "above 65 °C" for Taq.

### The digest simulator

`Enzyme-digest` is a C++ port of the Python restriction digest simulator, taken
from commit `82994c5` on `main` of
[our fork](https://github.com/Nheyer/restriction-enzyme-digest-simulator).
Upstream (`wl5e/restriction-enzyme-digest-simulator`) is archived and read-only,
so the fork's `main` is where that project's work actually lives.

The Python is **not** vendored here. It was carried as a submodule while the
port was being written and has been removed now that the port stands on its
own — `External_tools/` holds only build dependencies. Record the
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

`TEST_PCR_LAMBDA.fna` is bases 29,951–40,100 of phage lambda (NC_001416.1), the
stretch holding both of NEB's lambda control amplicons from the Phusion manual:
the 1.3 kb pair (30,006–31,352, 1347 bp) and the 10 kb pair (30,011–40,043,
10,033 bp). Position `p` in the file is genome position `p + 29,950`.
`TEST_PCR_PRIMERS.fna` is the 1.3 kb pair as a primer FASTA, the reverse primer
lower case and wrapped over two lines to exercise the parsing.

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
offered under different terms. Source files carry an
`SPDX-License-Identifier: GPL-3.0-only` line; the CSV and fixture files cannot
hold a header and are covered by this section instead.

Copyright © 2026 The Zabel Lab at Colorado State University, for the work in
this repository — including the C++ port of the digest simulator, which is a
substantial body of original work in its own right even though it is a
derivative of the code it was ported from.

### Third party attribution

`CppSrc/Enzyme-digest/Restriction_Enzyme_Digest.cpp` is a translation of the
Python restriction digest simulator, © Collins Amatu Gorgerat, which is MIT
licensed — [LICENSE.restriction-digest](LICENSE.restriction-digest). What
descends from Collins's original first commit (`8661193`) is `digest_linear`,
`digest_circular`, `print_fragment_table`, `draw_ascii_gel`, the shape of `main()`
with its original options, the 100 bp ladder, and the wording of `parse_fasta`'s
error messages. Reading the FASTA itself is now the shared
`CppSrc/common/fasta.hpp`, which grew out of the `MSA-to-consensus` reader and
takes nothing from the Python; `parse_fasta` in the port is only what wraps it.
The enzyme matching, NEB notation, IUPAC handling, CSV loading and the rest of
the CLI were added to the Python in the fork after that commit, or are new in the
port.

That is the whole of the MIT lineage; nothing else in the repository descends
from that code:

| Path | Where it comes from | MIT notice |
|---|---|---|
| `Restriction_Enzyme_Digest.cpp` | The translation above | Yes — the header block stays with the whole file |
| `enzyme_data.hpp` | The ladder is Collins's `DNA_LADDER_100BP`; the IUPAC codes and complements are the standard nomenclature | The ladder only, noted in the file |
| `CppSrc/common/fasta.hpp` | The FASTA reader all three tools share. Grew out of the `MSA-to-consensus` reader; nothing taken from the port's `parse_fasta` | No |
| `Tests/test_digest.cpp` | Case names and layout follow Collins's tests; the inputs and assertions are rewritten, apart from the two trivial no-cuts cases (`digest_linear` and `digest_circular` of 100 with no cuts, giving `[100]`), which are identical | A provenance note in the file |
| `Data/restriction_enzymes.csv` | NEB's commercially available specificities in NEB notation, compiled in the fork (`a2e48c4`) and checked against REBASE. Collins's original table was 20 enzymes in another format and none of it is in this file | **No** |
| Everything else — `MSA-to-consensus`, `PCR-protocol`, their tests, `polymerases.csv`, the fixtures, the build files | Written here, or from the datasheets and papers cited in the files | No |

This is attribution, not an alternative licence. MIT is GPL compatible, which is
exactly why that code can be absorbed here: the combined work is distributed
under the GPL v3 above, and the MIT notice travels with the portion it covers.
Keeping the notice is a condition of the permission that lets us use the code at
all, so it stays with those files, with the source header, and with anything
built or installed from them — regardless of how far the port diverges from the
original. A translation into another language is a derivative work; the extent
of the changes does not retire the obligation.

`PCR-protocol` links `thal.c`, `thal_parameters.c` and `oligotm.c` from
[primer3](https://github.com/primer3-org/primer3), © Whitehead Institute for
Biomedical Research, Steve Rozen, Andreas Untergasser and Helen Skaletsky
(years per file, 1996–2018), distributed under the GPL version 2 or (at your
option) any later version — which is what lets it be combined into this GPL v3 work. Its licence
and copyright headers are unchanged in the `External_tools/primer3`
submodule, and its licence is installed alongside the binaries as
`LICENSE.primer3`.
