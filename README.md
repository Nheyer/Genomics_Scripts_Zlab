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

The binaries land in `bin/MSA-to-consensus`, `bin/Enzyme-digest`,
`bin/PCR-protocol` and `bin/Digest-protocol`; build one on its own by naming it as the `--target`. The build uses the vendored
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
<prefix>/bin/Digest-protocol
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
| [`Digest-protocol`](#digest-protocol) | The bench side of a restriction digest — reaction setup, master mix, incubation and how to stop it, for one or more enzymes |

All four write their result to stdout (or, for the first, to a file), with
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

## Digest-protocol

The bench side of a restriction digest. `Enzyme-digest` says where an enzyme
cuts; this says what to pipette. Give it one or more enzymes and an amount of
DNA and it writes the reaction setup, a master mix for several tubes, the
incubation, and how to stop it.

```bash
./bin/Digest-protocol -e EcoRI
```

```
Restriction digest protocol: EcoRI

Enzymes
  EcoRI         G^AATT_C          10X NEBuffer EcoRI/SspI, 37 C

Reaction setup (50 uL, 1 ug DNA, 10 units of each enzyme)
  Component                           Volume      Final
  Nuclease-free water                 to 50 uL
  10X NEBuffer EcoRI/SspI             5 uL        1X
  DNA                                 variable    1 ug (20 ng/uL)
  EcoRI (10 U/uL)                     1 uL        10 U
  Enzyme(s) 1 uL of the 50 uL, 2.0% (NEB: at most 10%).

Procedure
  1. Add the water, buffer and DNA, and the enzyme last. Keep the enzyme on ice when it is not in the freezer.
  2. Mix by pipetting up and down, or by flicking the tube. Give it a quick spin in a microcentrifuge. Do not vortex.
  3. Incubate at 37 C for 60 min.
  4. Stop the reaction: with 10 uL of stop solution (NEB uses 10 uL per 50 uL reaction) if the DNA needs no further handling; if it does, heat inactivate or remove the enzyme with a spin column or phenol/chloroform extraction.
     Heat inactivate at 65 C for 20 min.

Notes
  Amounts: NEB's Restriction Digest protocol - buffer at 1X, 10 units per ug of DNA (20 for genomic), enzyme at most 10% of the volume. Source: New England Biolabs, Restriction Digest v2, protocols.io, DOI 10.17504/protocols.io.isycefw.
  DNA should be free of phenol, chloroform, alcohol, EDTA, detergents and excess salt; methylation can block some enzymes.
  Time-Saver Qualified (EcoRI): NEB says 5-15 min is enough; the incubation here is 60 min.
  EcoRI source: NEB NEBuffer Activity/Performance Chart with Restriction Enzymes (neb.com) - saved 2026-09-22; NEB Heat Inactivation (neb.com) - saved 2026-09-22

Warnings
  ! stock concentration of EcoRI not on file, so 10 U/uL is assumed (NEB's "10 units, generally 1 uL"); if your tube says otherwise, give --stock EcoRI=UNITS
```

Buffer, incubation temperature and full heat inactivation (temperature and
minutes) come from [Conditions data](#conditions-data) whenever there is a
row — 266 enzymes have one, so most digests need nothing extra. Stock
concentration is never on that row (see the warning). A double digest of two
enzymes that both have rows, and that agree on buffer, works the same way,
with a master mix for several tubes:

```bash
./bin/Digest-protocol -e EcoRI-HF,BamHI-HF -n 6 --dna-ug 2 --dna-conc 100 --volume 100
```

```
Reaction setup (100 uL, 2 ug DNA, 20 units of each enzyme)
  Component                           Volume      Final
  Nuclease-free water                 66 uL
  10X rCutSmart                       10 uL       1X
  DNA                                 20 uL       2 ug (20 ng/uL)
  EcoRI-HF (10 U/uL)                  2 uL        20 U
  BamHI-HF (10 U/uL)                  2 uL        20 U
  Enzyme(s) 4 uL of the 100 uL, 4.0% (NEB: at most 10%).

Master mix (6 reactions + 10% = 6.6, 528 uL)
  Component                           Per rxn     Mix
  10X rCutSmart                       10 uL       66 uL
  EcoRI-HF (10 U/uL)                  2 uL        13.2 uL
  BamHI-HF (10 U/uL)                  2 uL        13.2 uL
  Nuclease-free water                 66 uL       435.6 uL
  Aliquot 80 uL per tube, then add 20 uL DNA.
```

`BamHI-HF` cannot be heat inactivated at all (NEB lists none for it), which the
procedure step says; `EcoRI-HF`'s row gives a full 65 C for 20 min. When the
enzymes in a tube disagree on how, or whether, they can be heat inactivated,
each gets its own line rather than one line speaking for both. Give conditions
on the command line — `NAME:BUFFER:TEMP_C[:INACTIVATION]` — to use an enzyme
with no row, or to override one that has it:

```bash
./bin/Digest-protocol -e SomeNewEnzyme:rCutSmart:37:65/20 --stock SomeNewEnzyme=20
```

### Options

| Flag | Meaning |
|---|---|
| `-e`, `--enzymes` | Comma separated enzymes, each `NAME` or `NAME:BUFFER:TEMP_C[:INACTIVATION]` (required) — see below |
| `-d`, `--dna-ug` | DNA in µg per reaction (default: the most NEB pairs with the volume, 1 µg at 50 µL) |
| `--dna-conc` | DNA concentration in ng/µL; works out the volume of DNA per tube |
| `--dna-volume` | DNA per tube in µL, left out of the master mix, used when there is no `--dna-conc` (default 1) |
| `--volume` | Reaction volume in µL (default 50) |
| `-n`, `--replicates` | Number of reactions (default 1); above 1, a master mix is added with 10% extra |
| `--genomic` | Genomic DNA: 20 units per µg of DNA instead of 10 |
| `--units-per-ug` | Units of each enzyme per µg, instead of NEB's 10 (20 with `--genomic`) |
| `--stock` | An enzyme's concentration in U/µL as `NAME=UNITS`, repeatable; read it off the tube |
| `-t`, `--minutes` | Incubation time in minutes (default 60) |
| `--buffer` | The one 10X buffer for the whole reaction, when the enzymes' own do not settle it |
| `--incubate-c` | The one incubation temperature for the whole reaction |
| `--conditions-db` | Read per-enzyme conditions from a CSV instead of the built in table |
| `--list-enzymes` | List the known enzymes and which have conditions on file, and exit |

### Enzyme specifications

An enzyme is checked against the same 287 enzymes `Enzyme-digest` reads, case
insensitively, so a typo is an error and not a protocol for something that does
not exist. There are three forms:

| Form | Example | Meaning |
|---|---|---|
| `NAME` | `EcoRI` | Conditions from the table if there is a row, otherwise unknown |
| `NAME:BUFFER:TEMP_C` | `EcoRI:rCutSmart:37` | Buffer and incubation temperature as you read them off NEB's chart |
| `NAME:BUFFER:TEMP_C:INACTIVATION` | `EcoRI:rCutSmart:37:65/20` | Also heat inactivation: `TEMP_C/MINUTES`, `TEMP_C` alone if the time is not known, or `no` if NEB lists none |

The second and third forms also let you use an enzyme neither table lists, such
as a `-HF` version. Conditions given on the command line win over a table row,
and if they differ from it the output says so, since the row's source no longer
describes them.

Several enzymes make a double digest, and it needs one buffer and one
temperature. If the enzymes disagree the tool stops rather than pick one:
digest them one after the other, or give `--buffer` / `--incubate-c` once you
have checked a single choice works for all of them.

### How the numbers are worked out

Everything that is not per enzyme comes from one document, NEB's *Restriction
Digest* protocol (version 2, [protocols.io, DOI
10.17504/protocols.io.isycefw](https://dx.doi.org/10.17504/protocols.io.isycefw),
CC BY), read in full; `digest_data.hpp` cites each constant to it.

- **Buffer** at 1X from a 10X stock: a tenth of the volume.
- **Enzyme** at 10 units for every µg of DNA — the upper end of NEB's "5–10 units
  per µg", and 20 for genomic DNA, the upper end of its 10–20. Each enzyme of a
  double digest gets its own full dose. The volume comes from the stock
  concentration; where that is not known the tool takes NEB's "10 units,
  generally 1 µL", which is 10 U/µL, and says it has.
- **At most 10% of the volume is enzyme**, because the glycerol it is stored in
  causes star activity. More than that is an error, not a warning.
- **DNA** against volume follows NEB's table for smaller reactions: 0.1 µg in
  10 µL, 0.5 µg in 25 µL, 1 µg in 50 µL, and 1 µg per 50 µL above that. That
  table is what `--dna-ug` defaults to, and more than it is a warning.
- **Stop solution** at 10 µL per 50 µL of reaction, as NEB does.
- **Master mix** takes 10% extra for pipetting loss. That is this repo's
  convention, the same as `PCR-protocol`, and not NEB's.

The three rows NEB tabulates — 1 unit / 0.1 µg / 1 µL buffer at 10 µL, 5 / 0.5 /
2.5 at 25 µL, 10 / 1 / 5 at 50 µL — come out exactly, and a test holds them
there.

### Conditions data

`Data/digest_conditions.csv` holds what depends on the enzyme, in the same
arrangement as the enzyme and polymerase tables: embedded at configure time and
parsed by the loader `--conditions-db` uses. It holds 266 rows, built from two
of NEB's own reference pages the person saved from their own browser and
handed over as files — NEB's site refuses automated access (HTTP 403), and a
page pasted into chat isn't a file a parser can read reliably, so a saved page
was the way in both times. The rule for this repo is that every number comes
from the manufacturer and is cited, never typed from memory: `source` names
every page a row's data came from, and the save date.

The **NEBuffer Performance Chart** (saved 2026-09-22) gives buffer, incubation
temperature, heat-inactivation temperature and Time-Saver status, but never a
heat-inactivation *time*. The **Heat Inactivation** page (also saved
2026-09-22) gives a temperature and a time together, per enzyme. Its
temperature was checked against the chart's for all 266 enzymes the two pages
share before either was trusted to fill in the other — zero disagreements — so
the merge fills `inactivate_min` from a second, independent NEB source, not a
recalled default; every non-"No" time on that page reads exactly "20 minutes",
with no exception among its 280 rows. Fourteen enzymes it names use older,
non-versioned spellings with no match in the table (`BsaI-HF` where the chart
has `BsaI-HFv2`, and similar) and were left out rather than guess the mapping.

Add a row for an enzyme neither page lists and it is picked up on the next
`cmake .`. The loader refuses a row that is missing a column, has a
non-numeric number, or has no source. Not every chart NEB publishes carries
every column — stock concentration is on neither of these two — so three
columns have their own "not known" value rather than forcing a row to invent
one or be refused. Every row on file today has `inactivate_c`/`inactivate_min`
fully known (a real temperature and minutes, or `0,0` where NEB lists no heat
inactivation) and `0` for `stock_u_per_ul` (not on either page); the
"temperature known, time not" state below is what a single-source row would
use, and is what the tool falls back to if you add one from the chart alone:

| Column | Meaning |
|---|---|
| `enzyme` | Name, as NEB spells it; may be one the enzyme table lacks, such as `EcoRI-HF` |
| `buffer` | The 10X buffer NEB lists |
| `incubation_c` | Incubation temperature |
| `inactivate_c` | Heat inactivation temperature: `-1` not known, `0` NEB lists none, else the temperature |
| `inactivate_min` | Minutes to go with it: `-1` not known (the temperature above may still be), `0` only paired with `inactivate_c` `0`, else the minutes |
| `time_saver` | `1` if NEB lists it as Time-Saver Qualified, `0` if NEB lists it as not, `-1` if this source does not say |
| `stock_u_per_ul` | Concentration of the tube NEB supplies, U/µL; `0` if this source does not say |
| `notes` | Anything not straight off the page; printed with the protocol |
| `source` | Where NEB lists it. Required |

The same `-1`/temperature-alone forms work on `--enzymes`: `EcoRI:rCutSmart:37:65`
means the temperature is known and the time is not, same as `65/20` means both
are and `no` means NEB lists no heat inactivation.

### Known limitations

- **No stock concentration ships for any enzyme**, because neither of the two
  NEB pages the table is built from carries it. The output assumes 10 U/µL and
  warns whenever it does.
- **Buffer, temperature and heat inactivation are on file for 266 enzymes, not
  all of them.** An enzyme newer than the two pages, or a variant neither
  names, needs its conditions given on the command line; the output says so
  and names what is missing rather than guess.
- **The stock concentration is assumed** at 10 U/µL unless a row or `--stock`
  says otherwise, so the enzyme volume can be off by a factor of two for a
  20 U/µL stock. The output warns whenever it is doing this.
- **No methylation or star activity check per enzyme**, and no SAM supplement
  for the enzymes that need it. Those go in a row's `notes`; none of the rows
  on file today carry one.
- **It does not look at the DNA.** Whether the enzymes cut it, and where, is
  `Enzyme-digest`.
- **No plan for a sequential digest.** It refuses enzymes that disagree on buffer
  or temperature; it does not schedule the two steps.

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
| `CppSrc/Digest-protocol/Digest_Protocol.cpp` | `Digest-protocol` |
| `CppSrc/Digest-protocol/digest_data.hpp` | NEB's protocol constants, each with its source — lookup constants only, no logic |
| `CppSrc/Digest-protocol/digest_csv_embedded.hpp.in` | Template CMake fills with the conditions CSV at configure time; required to build |
| `CppSrc/common/fasta.hpp` | The FASTA reader shared by the three tools that read FASTA. Header only, not a target, no `main()` — each tool is still one translation unit |
| `Data/restriction_enzymes.csv` | The enzyme database, in NEB `^`/`_` notation. **The source of truth** — see [Enzyme data](#enzyme-data) |
| `Data/polymerases.csv` | The polymerase table, from the datasheets. **The source of truth** — see [Polymerase data](#polymerase-data) |
| `Data/digest_conditions.csv` | Per-enzyme digest conditions for 266 enzymes, from two of NEB's reference pages. **The source of truth** — see [Conditions data](#conditions-data) |
| `Data/test_files/` | Known truth fixtures, see [Test files](#test-files) |
| `Tests/test_consensus.cpp` | CUnit unit tests, see [Running the tests](#running-the-tests) |
| `Tests/test_digest.cpp` | CUnit unit tests for `Enzyme-digest`, same |
| `Tests/test_pcr.cpp` | CUnit unit tests for `PCR-protocol`, same |
| `Tests/test_digest_protocol.cpp` | CUnit unit tests for `Digest-protocol`, same |
| `Tests/test_fasta.cpp` | CUnit unit tests for the shared reader; it includes the header directly, there is no `main()` to rename |
| `LICENSE.restriction-digest` | MIT notice for the code `Enzyme-digest` was ported from, see [License](#license) |
| `External_tools/argparse` | [p-ranav/argparse](https://github.com/p-ranav/argparse), header only CLI parsing (build dependency) |
| `External_tools/cunit` | [cunity/cunit](https://gitlab.com/cunity/cunit), unit test framework (build dependency) |
| `External_tools/primer3` | [primer3-org/primer3](https://github.com/primer3-org/primer3) v2.6.1; `PCR-protocol` builds its `thal` hairpin/dimer aligner (build dependency) |

## Development

### Running the tests

`Tests/` unit tests all four tools with CUnit:

```bash
cmake --build . --target test_consensus test_digest test_pcr test_fasta test_digest_protocol
./bin/test_consensus
./bin/test_digest
./bin/test_pcr
./bin/test_fasta
./bin/test_digest_protocol
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

The digest protocol suite pins the amounts to NEB's own table for 10, 25 and
50 µL, typed in the test separately from `digest_data.hpp`, and works everything
else out by hand: the water in a typical digest (50 − 5 − 1 − 1 = 43 µL), the
exact boundary of the 10% enzyme limit, and a master mix down to 189.2 µL of
water. It also holds the conditions loader to its rules — no row without a
source — and checks the enzyme specs and what a double digest can and cannot
agree on.

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
| `CppSrc/common/fasta.hpp` | The FASTA reader shared by the three tools that read FASTA. Grew out of the `MSA-to-consensus` reader; nothing taken from the port's `parse_fasta` | No |
| `Tests/test_digest.cpp` | Case names and layout follow Collins's tests; the inputs and assertions are rewritten, apart from the two trivial no-cuts cases (`digest_linear` and `digest_circular` of 100 with no cuts, giving `[100]`), which are identical | A provenance note in the file |
| `Data/restriction_enzymes.csv` | NEB's commercially available specificities in NEB notation, compiled in the fork (`a2e48c4`) and checked against REBASE. Collins's original table was 20 enzymes in another format and none of it is in this file | **No** |
| Everything else — `MSA-to-consensus`, `PCR-protocol`, `Digest-protocol`, their tests, `polymerases.csv`, `digest_conditions.csv`, the fixtures, the build files | Written here, or from the datasheets and papers cited in the files | No |

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
