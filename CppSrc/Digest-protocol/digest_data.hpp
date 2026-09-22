// SPDX-License-Identifier: GPL-3.0-only
//
// Static reference data for the restriction digest protocol generator.
//
// Lookup constants only - no logic, and nothing here depends on the tool, the
// same split as enzyme_data.hpp and pcr_data.hpp. Every number carries its
// source, and none was typed from memory: a wrong number here does not crash
// anything, it just quietly ruins someone's digest.
//
// The source for all of it is one document, read in full:
//
//   New England Biolabs, "Restriction Digest" (version 2), protocols.io,
//   DOI 10.17504/protocols.io.isycefw. Its Guidelines tab holds the 'typical'
//   digest, the table for smaller volumes and the rules quoted below; its
//   protocol tab holds the steps. Open access, CC BY.
//
// Nothing per enzyme is here. Buffer, incubation temperature and heat
// inactivation depend on the enzyme, and NEB's own protocol sends the reader to
// its Activity/Performance Chart for them. They live in Data/digest_conditions.csv
// instead, one cited row per enzyme, embedded at configure time.
#ifndef DIGEST_DATA_HPP
#define DIGEST_DATA_HPP

namespace digest_data {

// 'A "Typical" Restriction Digest': enzyme 10 units ("generally 1 ul is used"),
// DNA 1 ug, 10X NEBuffer 5 ul (1X), total 50 ul, 1 hour. The table for smaller
// volumes gives 1 unit / 0.1 ug / 1 ul at 10 ul, 5 units / 0.5 ug / 2.5 ul at
// 25 ul, and 10 units / 1 ug / 5 ul at 50 ul. Read down the rows that is a
// buffer that is one tenth of the volume, and 10 units for every ug of DNA.
inline constexpr double BUFFER_X = 10.0;
inline constexpr double UNITS_PER_UG = 10.0;

// The largest DNA mass NEB pairs with each volume it tabulates. Not linear: the
// 10 ul row is 0.1 ug, half of what the other two rows' ratio would give,
// because "additives in the restriction enzyme storage buffer ... as well as
// contaminants ... can be problematic in smaller reaction volumes".
struct dna_limit {
    double volume_ul;
    double max_ug;
};
inline constexpr dna_limit DNA_LIMITS[] = {
    {10.0, 0.1},
    {25.0, 0.5},
    {50.0, 1.0},
};
inline constexpr int DNA_LIMIT_COUNT = 3;
// Above the last row, "a 50 ul reaction volume is recommended for digestion of
// 1 ug of substrate", read as 1 ug per 50 ul: "this enzyme : DNA : reaction
// volume ratio can be used as a guide when designing reactions".
inline constexpr double DNA_UG_PER_UL_ABOVE_TABLE = 1.0 / 50.0;

// "In general, we recommend 5-10 units of enzyme per ug DNA, and 10-20 units for
// genomic DNA in a 1 hour digest." The tool takes the upper end of each range,
// the way the polymerase table takes the more conservative end of a datasheet's.
inline constexpr double UNITS_PER_UG_MIN = 5.0;
inline constexpr double GENOMIC_UNITS_PER_UG = 20.0;

// "Enzyme volume should not exceed 10% of the total reaction volume to prevent
// star activity due to excess glycerol."
inline constexpr double ENZYME_MAX_FRACTION = 0.10;

// "10 units is sufficient, generally 1 ul is used": that is 10 units per ul, and
// it is what the tool assumes for an enzyme whose stock it has not been told.
// It is the protocol's own ratio, not a fact about any enzyme, so the output
// says when it is being used.
inline constexpr double DEFAULT_STOCK_U_PER_UL = 10.0;

// "Incubation time is typically 1 hour." Time-Saver Qualified enzymes: "can be
// decreased to 5-15 minutes". "It is possible, with many enzymes, to use fewer
// units and digest for up to 16 hours."
inline constexpr int DEFAULT_MINUTES = 60;
inline constexpr int TIME_SAVER_MIN_MINUTES = 5;
inline constexpr int TIME_SAVER_MAX_MINUTES = 15;
inline constexpr int EXTENDED_DIGEST_HOURS = 16;

// "10 ul rxns should not be incubated for longer than 1 hour to avoid
// evaporation."
inline constexpr double SMALL_REACTION_UL = 10.0;

// "Terminate with a stop solution (10 ul per 50 ul rxn)", e.g. NEB #B7024.
inline constexpr double STOP_UL_PER_50UL = 10.0;

}  // namespace digest_data

#endif  // DIGEST_DATA_HPP
