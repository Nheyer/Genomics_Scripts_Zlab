// SPDX-License-Identifier: GPL-3.0-only
//
// Static reference data for the PCR protocol generator.
//
// Lookup tables only - no logic, and nothing here depends on the tool, the same
// split as enzyme_data.hpp. Every number carries its source. Thermodynamic
// parameters that are wrong do not crash anything, they just move a Tm, so none
// of these were typed from memory.
//
// The polymerase table is deliberately NOT here: Data/polymerases.csv is the
// source of truth for it, embedded verbatim at configure time (see
// polymerase_csv_embedded.hpp.in).
#ifndef PCR_DATA_HPP
#define PCR_DATA_HPP

namespace pcr_data {

// Unified Watson-Crick nearest-neighbour parameters in 1 M NaCl.
// SantaLucia J Jr. (1998) PNAS 95:1460-1465, Table 2,
// https://doi.org/10.1073/pnas.95.4.1460 (PMC19045).
//
// Each stack is written as 5'-XY-3' on the top strand; the paper writes it as
// XY/X'Y', so AA/TT is "AA", CA/GT is "CA". The other six stacks are the same
// duplexes read from the bottom strand (TT = AA, AC = GT, and so on), which
// nn_stack() in the tool resolves by reverse complementing.
struct nn_row {
    const char *stack;
    double dH;  // kcal/mol
    double dS;  // cal/(K mol)
};
inline constexpr nn_row NN_UNIFIED[] = {
    {"AA",  -7.9, -22.2},
    {"AT",  -7.2, -20.4},
    {"TA",  -7.2, -21.3},
    {"CA",  -8.5, -22.7},
    {"GT",  -8.4, -22.4},
    {"CT",  -7.8, -21.0},
    {"GA",  -8.2, -22.2},
    {"CG", -10.6, -27.2},
    {"GC",  -9.8, -24.4},
    {"GG",  -8.0, -19.9},
};
inline constexpr int NN_UNIFIED_COUNT = 10;

// Same paper, same table: duplex initiation, charged once per terminal pair,
// and the entropic penalty for a self-complementary duplex.
inline constexpr double INIT_GC_DH = 0.1;
inline constexpr double INIT_GC_DS = -2.8;
inline constexpr double INIT_AT_DH = 2.3;
inline constexpr double INIT_AT_DS = 4.1;
inline constexpr double SYMMETRY_DS = -1.4;

// Salt dependence of the entropy (SantaLucia 1998):
//   dS([Na+]) = dS(1 M) + 0.368 * N * ln[Na+]
// with N the number of phosphates in the duplex divided by 2, i.e. the number
// of stacks. primer3's oligotm.c implements it as 0.368 * (len - 1).
inline constexpr double SALT_DS_PER_STACK = 0.368;

// Magnesium as sodium equivalents, [Na+]eq = [Mono] + 120 * sqrt([Mg] - [dNTP])
// in mM, exactly as primer3 does it (oligotm.c, divalent_to_monovalent). dNTPs
// chelate Mg one to one, so only the free magnesium counts; with more dNTP
// than Mg none is left.
inline constexpr double MG_TO_NA_FACTOR = 120.0;

// Gas constant in cal/(K mol), the value primer3 uses, and 0 C in kelvin.
inline constexpr double GAS_CONSTANT = 1.987;
inline constexpr double KELVIN = 273.15;

// Hairpins and dimers are not scored here: they come from primer3's thal
// (External_tools/primer3), which handles mismatches, bulges, internal loops
// and dangling ends. It reports the melting temperature of the most stable
// structure, and primer3 flags one above 47 C by default for all five
// of PRIMER_MAX_SELF_ANY_TH, PRIMER_MAX_SELF_END_TH, PRIMER_PAIR_MAX_COMPL_ANY_TH,
// PRIMER_PAIR_MAX_COMPL_END_TH and PRIMER_MAX_HAIRPIN_TH (primer3 manual,
// External_tools/primer3/src/primer3_manual.htm).
inline constexpr double STRUCTURE_TM_LIMIT = 47.0;

// Primer design limits. Premier Biosoft, "PCR Primer Design Guidelines",
// http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html:
//   - a Tm difference of 5 C or more between the pair can lead to no product
//   - no more than 3 G or C in the last 5 bases
//   - runs of one base up to 4, dinucleotide repeats up to 4
// GC content of 40-60% and length of 20-40 nt are from the NEB datasheets for
// all three built-in polymerases (Premier gives the same GC range).
inline constexpr double TM_DIFFERENCE_LIMIT = 5.0;
inline constexpr int GC_CLAMP_WINDOW = 5;
inline constexpr int GC_CLAMP_MAX = 3;
inline constexpr int MAX_RUN = 4;
inline constexpr int MAX_DINUCLEOTIDE_REPEATS = 4;
inline constexpr double GC_MIN_PERCENT = 40.0;
inline constexpr double GC_MAX_PERCENT = 60.0;
inline constexpr int PRIMER_MIN_NT = 20;
inline constexpr int PRIMER_MAX_NT = 40;

}  // namespace pcr_data

#endif  // PCR_DATA_HPP
