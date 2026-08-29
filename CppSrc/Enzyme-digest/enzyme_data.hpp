// Static reference data for the restriction enzyme digest simulator.
//
// Lookup tables only - no logic, and nothing here depends on the simulator, so
// this header can be read, diffed and regenerated on its own. That mirrors the
// enzyme_data.py / enzyme_digest.py split in the Python original, which has a
// test asserting the data module holds no functions and does not import the
// logic. Keep that direction of dependency.
//
// The enzyme table itself is deliberately NOT here: restriction_enzymes.csv is
// the source of truth for it, embedded verbatim at configure time (see
// enzyme_csv_embedded.hpp.in) and parsed by load_enzyme_csv_string(). Encoding
// the enzymes as C++ literals as well would give two copies to drift apart.
#ifndef ENZYME_DATA_HPP
#define ENZYME_DATA_HPP

namespace enzyme_data {

// The 15 IUPAC nucleotide codes, and the set of bases each one stands for.
// Base strings are in ACGT order; expand_ambiguity() emits variants in this
// order, so changing it changes that output.
struct iupac_row {
    char code;
    const char *bases;
};
inline constexpr iupac_row IUPAC_BASES[] = {
    {'A', "A"},    {'C', "C"},    {'G', "G"},    {'T', "T"},
    {'R', "AG"},   {'Y', "CT"},   {'S', "CG"},   {'W', "AT"},
    {'K', "GT"},   {'M', "AC"},
    {'B', "CGT"},  {'D', "AGT"},  {'H', "ACT"},  {'V', "ACG"},
    {'N', "ACGT"},
};
inline constexpr int IUPAC_COUNT = 15;

// Iteration order for the codes. The Python original derives this from the
// insertion order of its IUPAC_BASES dict, and site_to_regex() picks the first
// matching code, so the order is load bearing - keep it as is.
inline constexpr char IUPAC_ALPHABET[] = "ACGTRYSWKMBDHVN";

// Order used when telling the user which characters are allowed. Cosmetic.
inline constexpr char IUPAC_DISPLAY_ORDER[] = "ACGTRYKMSWBDHVN";

// Without --ambiguity this is the whole accepted input alphabet.
inline constexpr char CONCRETE_BASES[] = "ACGT";
inline constexpr char GAP_CHAR = '-';

// Complementing a code complements the bases it stands for: R (A or G) pairs
// with Y (T or C), B (not A) with V (not T), and so on.
inline constexpr char COMPLEMENT_FROM[] = "ACGTRYSWKMBDHVN";
inline constexpr char COMPLEMENT_TO[]   = "TGCAYRSWMKVHDBN";

// Band sizes of a standard 100 bp DNA ladder, for the simulated gel.
inline constexpr int DNA_LADDER_100BP[] = {
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1200, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000
};
inline constexpr int DNA_LADDER_100BP_COUNT = 19;

}  // namespace enzyme_data

#endif  // ENZYME_DATA_HPP
