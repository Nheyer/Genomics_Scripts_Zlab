// Restriction digest protocol generator - turns one or more enzymes and an amount
// of DNA into a reaction setup, an incubation and a way to stop it.
//
// Copyright (c) 2026 The Zabel Lab at Colorado State University.
//
// This file is part of Genomics_Scripts_Zlab and is distributed under the
// GPL v3 - see the LICENSE file at the root of the repository.
//
// SPDX-License-Identifier: GPL-3.0-only
//
// What it computes, and where each number comes from:
//
//   Setup       NEB's 'typical' digest and its table for smaller volumes, from
//               the Restriction Digest protocol (see digest_data.hpp for the
//               citation and for every constant): buffer at one tenth of the
//               volume, 10 units of each enzyme per ug of DNA (20 for genomic
//               DNA), enzyme never more than 10% of the volume. The amounts
//               reproduce NEB's 10, 25 and 50 uL rows exactly (pinned by a test).
//   Per enzyme  Buffer, incubation temperature, heat inactivation and stock
//               concentration come from Data/digest_conditions.csv, one cited
//               row per enzyme, or from the command line. NEB's own protocol
//               says only "enzyme dependent" and points at its chart, so the
//               tool has no default for these. Where it does not know, it says
//               so and does not guess: a wrong temperature does not crash, it
//               just gives an undigested plasmid.
//   Enzymes     Names are checked against the same enzyme table Enzyme-digest
//               uses (Data/restriction_enzymes.csv), so a typo is an error and
//               not a protocol for an enzyme that does not exist.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <argparse/argparse.hpp>

#include "digest_csv_embedded.hpp"
#include "digest_data.hpp"
#include "enzyme_csv_embedded.hpp"

using namespace digest_data;

// Single source of truth for the version.
static const char *const VERSION = "0.1.0";

// ---------------------------------------------------------------- strings ---

static std::string upper_of(const std::string &s) {
    std::string out = s;
    for (char &c : out) { c = (char) toupper((unsigned char) c); }
    return out;
}

static std::string trimmed(const std::string &s) {
    size_t first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) { return ""; }
    size_t last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

// Plain split on one character. The CSVs have no quoting, and notes use '-' and
// ';' instead of commas so that stays true.
static std::vector<std::string> split_on(const std::string &line, char sep) {
    std::vector<std::string> fields;
    std::string current;
    for (char c : line) {
        if (c == sep) { fields.push_back(current); current.clear(); }
        else if (c != '\r') { current += c; }
    }
    fields.push_back(current);
    return fields;
}

static double parse_number(const std::string &text, const std::string &column,
                           const std::string &where) {
    char *end = NULL;
    double value = std::strtod(text.c_str(), &end);
    if (text.empty() || *end != '\0') {
        throw std::invalid_argument(where + "column " + column + " is not a number: '" +
                                    text + "'");
    }
    return value;
}

static int parse_integer(const std::string &text, const std::string &column,
                         const std::string &where) {
    double value = parse_number(text, column, where);
    if (value != std::floor(value)) {
        throw std::invalid_argument(where + "column " + column + " must be a whole number: '" +
                                    text + "'");
    }
    return (int) value;
}

// ------------------------------------------------------------- conditions ---

// One row of the conditions table: what NEB lists for one enzyme. Not every
// chart NEB publishes carries every column - the NEBuffer Performance Chart
// gives a heat-inactivation temperature with no time, and neither stock
// concentration nor Time-Saver status at all - so each of those three fields
// has its own "not known" sentinel rather than forcing a row to either lie or
// be refused. incubation_c and buffer have no such state: a row with neither
// is not a conditions row, it is nothing.
struct Conditions {
    std::string name, buffer, notes, source;
    int incubation_c = 0;
    int inactivate_c = -1;       // -1 not known, 0 NEB lists no heat inactivation, else the temperature
    int inactivate_min = -1;     // -1 not known (temperature may still be); else the minutes
    int time_saver = -1;         // -1 not known, 0 not qualified, 1 Time-Saver Qualified
    double stock_u_per_ul = 0;   // 0 not known, else U/uL
};

static const char *const CONDITIONS_COLUMNS[] = {
    "enzyme", "buffer", "incubation_c", "inactivate_c", "inactivate_min", "time_saver",
    "stock_u_per_ul", "notes", "source",
};

// Read a conditions CSV with the columns above, in any order. Every column is
// required, and so is a source on every row: a row nobody can trace to NEB's
// chart is exactly the number this table exists to keep out. A table with no
// rows is fine - it means nothing is on file, not that the file is broken.
static std::vector<Conditions> load_conditions_csv_string(const std::string &text,
                                                          const std::string &path) {
    std::istringstream stream(text);
    std::string line;
    if (!std::getline(stream, line)) {
        throw std::invalid_argument("Nothing to read in " + path);
    }
    std::vector<std::string> header = split_on(line, ',');
    std::vector<int> col;
    std::string missing;
    for (const char *name : CONDITIONS_COLUMNS) {
        int found = -1;
        for (size_t i = 0; i < header.size(); i++) {
            if (trimmed(header[i]) == name) { found = (int) i; }
        }
        if (found < 0) { missing += (missing.empty() ? "" : ", ") + std::string(name); }
        col.push_back(found);
    }
    if (!missing.empty()) {
        throw std::invalid_argument(path + " is missing column(s) " + missing);
    }

    std::vector<Conditions> rows;
    int lineno = 1;
    while (std::getline(stream, line)) {
        lineno++;
        if (trimmed(line).empty()) { continue; }
        std::vector<std::string> row = split_on(line, ',');
        std::ostringstream where_stream;
        where_stream << path << ", line " << lineno << ": ";
        const std::string where = where_stream.str();
        int c = 0;
        auto text = [&](void) -> std::string {
            int idx = col[c];
            std::string value = idx < (int) row.size() ? trimmed(row[idx]) : "";
            c++;
            return value;
        };
        auto number = [&](void) -> double {
            const char *name = CONDITIONS_COLUMNS[c];
            return parse_number(text(), name, where);
        };
        auto integer = [&](void) -> int {
            const char *name = CONDITIONS_COLUMNS[c];
            return parse_integer(text(), name, where);
        };
        // Same order as CONDITIONS_COLUMNS.
        Conditions r;
        r.name = text();
        r.buffer = text();
        r.incubation_c = integer();
        r.inactivate_c = integer();
        r.inactivate_min = integer();
        r.time_saver = integer();
        r.stock_u_per_ul = number();
        r.notes = text();
        r.source = text();
        if (r.name.empty()) { throw std::invalid_argument(where + "enzyme name is empty"); }
        if (r.buffer.empty()) { throw std::invalid_argument(where + "buffer is empty"); }
        if (r.source.empty()) {
            throw std::invalid_argument(where + "source is empty: every row must cite where NEB "
                                                "lists it");
        }
        if (r.incubation_c <= 0) {
            throw std::invalid_argument(where + "incubation_c must be positive");
        }
        // -1/-1 not known; 0/0 NEB lists no heat inactivation; positive/-1 a
        // temperature with no time on this source; positive/positive both known.
        // Anything else is a half fact, which is worse than an admitted unknown.
        if (r.inactivate_c < -1 || r.inactivate_min < -1) {
            throw std::invalid_argument(where + "inactivate_c and inactivate_min cannot be less "
                                                "than -1 (-1 means not known)");
        }
        if (r.inactivate_c == -1 && r.inactivate_min != -1) {
            throw std::invalid_argument(where + "inactivate_c is -1 (not known), so inactivate_min "
                                                "must be -1 too");
        }
        if (r.inactivate_c == 0 && r.inactivate_min != 0) {
            throw std::invalid_argument(where + "inactivate_c is 0 (NEB lists no heat inactivation), "
                                                "so inactivate_min must be 0 too");
        }
        if (r.inactivate_c > 0 && r.inactivate_min == 0) {
            throw std::invalid_argument(where + "inactivate_c is positive but inactivate_min is 0; "
                                                "write -1 if the time is not known, or the minutes");
        }
        if (r.time_saver != -1 && r.time_saver != 0 && r.time_saver != 1) {
            throw std::invalid_argument(where + "time_saver must be -1 (not known), 0 or 1");
        }
        if (r.stock_u_per_ul < 0) {
            throw std::invalid_argument(where + "stock_u_per_ul cannot be negative; write 0 if not "
                                                "known");
        }
        for (const Conditions &seen : rows) {
            if (upper_of(seen.name) == upper_of(r.name)) {
                throw std::invalid_argument(where + "'" + r.name + "' is listed twice");
            }
        }
        rows.push_back(r);
    }
    return rows;
}

static std::vector<Conditions> load_conditions_csv(const std::string &path) {
    std::ifstream handle(path.c_str());
    if (!handle) {
        throw std::runtime_error("No such file or directory: '" + path + "'");
    }
    std::ostringstream buffer;
    buffer << handle.rdbuf();
    return load_conditions_csv_string(buffer.str(), path);
}

static const std::vector<Conditions> &builtin_conditions() {
    static std::vector<Conditions> table =
        load_conditions_csv_string(EMBEDDED_DIGEST_CSV, "digest_conditions.csv");
    return table;
}

// --------------------------------------------------------------- enzymes ---

// The enzyme table Enzyme-digest reads, for the names and the sites. Only those
// two columns are needed here, so this does not use Enzyme-digest's loader: this
// tool has no use for a cut position.
struct SiteRow {
    std::string name, site;
};

static const std::vector<SiteRow> &enzyme_sites() {
    static std::vector<SiteRow> table = [] {
        std::vector<SiteRow> out;
        std::istringstream stream(EMBEDDED_ENZYME_CSV);
        std::string line;
        std::getline(stream, line);   // the header
        while (std::getline(stream, line)) {
            std::vector<std::string> f = split_on(line, ',');
            if (f.size() >= 2 && !trimmed(f[0]).empty()) {
                out.push_back({trimmed(f[0]), trimmed(f[1])});
            }
        }
        return out;
    }();
    return table;
}

// One enzyme as it goes into this run. Whatever is not known stays at its
// sentinel and the output says so.
struct Enzyme {
    std::string name;
    std::string site;              // NEB notation, "" if the enzyme table lacks it
    std::string buffer;            // "" unknown
    int incubation_c = 0;          // 0 unknown
    int inactivate_c = -1;         // -1 unknown, 0 none listed, else the temperature
    int inactivate_min = -1;       // -1 unknown (the temperature above may still be known)
    int time_saver = -1;           // -1 unknown
    double stock_u_per_ul = 0;     // 0 unknown
    std::string notes, source;
    bool from_table = false;       // a conditions row supplied some of the above
    bool given = false;            // NAME:BUFFER:TEMP_C on the command line supplied conditions
    bool overrides_table = false;  // ... and they differ from the row, so its source no longer applies
};

static double effective_stock(const Enzyme &e) {
    return e.stock_u_per_ul > 0 ? e.stock_u_per_ul : DEFAULT_STOCK_U_PER_UL;
}

// "65/20", "65" alone (temperature known, time not), or "no". Anything else is
// refused rather than read as something.
static void parse_inactivation(const std::string &text, const std::string &spec, Enzyme *e) {
    if (upper_of(text) == "NO") {
        e->inactivate_c = 0;
        e->inactivate_min = 0;
        return;
    }
    std::vector<std::string> f = split_on(text, '/');
    if (f.size() == 1) {
        e->inactivate_c = parse_integer(trimmed(f[0]), "inactivation temperature", "'" + spec + "': ");
        if (e->inactivate_c <= 0) {
            throw std::invalid_argument("in '" + spec + "', heat inactivation temperature must be "
                                        "positive; write 'no' for none");
        }
        e->inactivate_min = -1;   // temperature known, no time given for it
        return;
    }
    if (f.size() != 2) {
        throw std::invalid_argument("in '" + spec + "', heat inactivation is '" + text +
                                    "'; write TEMP_C/MINUTES, e.g. 65/20, TEMP_C alone if the time "
                                    "is not known, or 'no'");
    }
    e->inactivate_c = parse_integer(trimmed(f[0]), "inactivation temperature", "'" + spec + "': ");
    e->inactivate_min = parse_integer(trimmed(f[1]), "inactivation minutes", "'" + spec + "': ");
    if (e->inactivate_c <= 0 || e->inactivate_min <= 0) {
        throw std::invalid_argument("in '" + spec + "', heat inactivation needs a positive "
                                    "temperature and time; write 'no' for none");
    }
}

// NAME, or NAME:BUFFER:TEMP_C, or NAME:BUFFER:TEMP_C:INACTIVATION. Values given on
// the command line win over the table row. Without a table row they are the only
// source, which is what lets an enzyme the table does not list be used at all.
static Enzyme parse_enzyme_spec(const std::string &raw, const std::vector<Conditions> &table) {
    std::vector<std::string> f = split_on(raw, ':');
    for (std::string &part : f) { part = trimmed(part); }
    if (f.size() > 4) {
        throw std::invalid_argument("enzyme '" + raw + "' has too many fields; use "
                                    "NAME[:BUFFER:TEMP_C[:INACTIVATION]]");
    }
    if (f.size() == 2) {
        throw std::invalid_argument("enzyme '" + raw + "': give both a buffer and a temperature "
                                    "(NAME:BUFFER:TEMP_C), or just the name");
    }
    if (f[0].empty()) { throw std::invalid_argument("an enzyme has no name"); }

    Enzyme e;
    e.name = f[0];
    bool known = false;
    const std::string wanted = upper_of(f[0]);
    for (const Conditions &r : table) {
        if (upper_of(r.name) != wanted) { continue; }
        e.name = r.name;
        e.buffer = r.buffer;
        e.incubation_c = r.incubation_c;
        e.inactivate_c = r.inactivate_c;
        e.inactivate_min = r.inactivate_min;
        e.time_saver = r.time_saver;
        e.stock_u_per_ul = r.stock_u_per_ul;
        e.notes = r.notes;
        e.source = r.source;
        e.from_table = true;
        known = true;
    }
    for (const SiteRow &s : enzyme_sites()) {
        if (upper_of(s.name) != wanted) { continue; }
        if (!e.from_table) { e.name = s.name; }
        e.site = s.site;
        known = true;
    }
    if (f.size() >= 3) {
        const int temperature = parse_integer(f[2], "temperature", "enzyme '" + raw + "': ");
        if (f[1].empty() || temperature <= 0) {
            throw std::invalid_argument("enzyme '" + raw + "': the buffer cannot be empty and the "
                                        "temperature must be positive");
        }
        // Only a real difference means the row's source no longer describes what
        // is being printed; typing the table's own values back changes nothing.
        if (e.from_table) {
            e.overrides_table = upper_of(f[1]) != upper_of(e.buffer) || temperature != e.incubation_c;
        }
        e.buffer = f[1];
        e.incubation_c = temperature;
        e.given = true;
        known = true;
        if (f.size() == 4) {
            const int before_c = e.inactivate_c, before_min = e.inactivate_min;
            parse_inactivation(f[3], raw, &e);
            if (e.from_table && (e.inactivate_c != before_c || e.inactivate_min != before_min)) {
                e.overrides_table = true;
            }
        }
    }
    if (!known) {
        throw std::invalid_argument("unknown enzyme '" + f[0] + "': not in the enzyme table or the "
                                    "conditions table (see --list-enzymes). To use one that is "
                                    "not, give its conditions: " + f[0] + ":BUFFER:TEMP_C");
    }
    return e;
}

// ----------------------------------------------------------------- output ---

static std::string fixed(double value, int places) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(places) << value;
    return out.str();
}

// Numbers with at most three decimals and no trailing zeros: 2.5, 0.25, 10.
static std::string volume(double ul) {
    std::string s = fixed(ul, 3);
    s.erase(s.find_last_not_of('0') + 1);
    if (!s.empty() && s.back() == '.') { s.pop_back(); }
    return s;
}

static std::string pad(const std::string &s, size_t width) {
    return s.size() >= width ? s + " " : s + std::string(width - s.size(), ' ');
}

// --------------------------------------------------------------- reaction ---

// Extra made up in a master mix to cover pipetting loss, as a fraction. The
// same figure PCR-protocol uses; it is this repo's convention, not NEB's.
static const double MASTER_MIX_OVERAGE = 0.10;

// One reagent of the reaction: what goes in one tube and what it comes to.
struct Component {
    std::string name;
    double ul = 0;      // per reaction
    std::string final;  // final amount, for display
};

// The largest DNA mass NEB pairs with this volume. At and above its 50 uL row
// that is 1 ug per 50 uL; between rows it is the row below; under the smallest
// row it is 0, because NEB says nothing about it.
static double max_dna_ug(double volume_ul) {
    if (volume_ul >= DNA_LIMITS[DNA_LIMIT_COUNT - 1].volume_ul) {
        return volume_ul * DNA_UG_PER_UL_ABOVE_TABLE;
    }
    double best = 0;
    for (int i = 0; i < DNA_LIMIT_COUNT; i++) {
        if (DNA_LIMITS[i].volume_ul <= volume_ul) { best = DNA_LIMITS[i].max_ug; }
    }
    return best;
}

struct Reaction {
    double volume_ul = 0;
    double dna_ug = 0;
    double dna_ul = 0;            // per tube
    double units = 0;             // of each enzyme
    std::vector<Component> parts; // the 10X buffer, then each enzyme
    double water_ul = 0;
    double enzyme_ul = 0;         // all the enzymes together
};

// Buffer, then every enzyme at units_per_ug for each ug of DNA, with the water
// that makes it up. Refuses what NEB says not to do, and what does not fit.
static Reaction plan_reaction(const std::vector<Enzyme> &enzymes, const std::string &buffer_name,
                              double volume_ul, double dna_ug, double dna_ul,
                              double units_per_ug) {
    if (volume_ul <= 0) { throw std::invalid_argument("--volume must be positive"); }
    if (dna_ug <= 0) { throw std::invalid_argument("the DNA amount must be positive"); }
    if (dna_ul < 0) { throw std::invalid_argument("--dna-volume cannot be negative"); }
    if (units_per_ug <= 0) { throw std::invalid_argument("--units-per-ug must be positive"); }
    if (enzymes.empty()) { throw std::invalid_argument("no enzymes"); }

    Reaction r;
    r.volume_ul = volume_ul;
    r.dna_ug = dna_ug;
    r.dna_ul = dna_ul;
    r.units = units_per_ug * dna_ug;
    r.parts.push_back({std::to_string((int) BUFFER_X) + "X " + buffer_name, volume_ul / BUFFER_X,
                       "1X"});
    for (const Enzyme &e : enzymes) {
        const double stock = effective_stock(e);
        const double ul = r.units / stock;
        r.enzyme_ul += ul;
        r.parts.push_back({e.name + " (" + volume(stock) + " U/uL)", ul, volume(r.units) + " U"});
    }
    const double allowed = ENZYME_MAX_FRACTION * volume_ul;
    if (r.enzyme_ul > allowed + 1e-9) {
        throw std::invalid_argument(
            "the enzymes come to " + volume(r.enzyme_ul) + " uL, more than the " +
            volume(allowed) + " uL (10%) NEB allows of a " + volume(volume_ul) +
            " uL reaction, since the glycerol they are stored in causes star activity; use a "
            "larger --volume, less DNA, or fewer --units-per-ug");
    }
    double used = dna_ul;
    for (const Component &c : r.parts) { used += c.ul; }
    r.water_ul = volume_ul - used;
    if (r.water_ul < -1e-9) {
        throw std::invalid_argument(
            "a " + volume(dna_ul) + " uL DNA sample does not fit: the reagents already take " +
            volume(used - dna_ul) + " of the " + volume(volume_ul) + " uL");
    }
    if (r.water_ul < 0) { r.water_ul = 0; }
    return r;
}

// A master mix for `replicates` reactions plus MASTER_MIX_OVERAGE. The DNA goes
// into each tube on its own, so the mix holds everything else, with water up to
// the volume less the DNA: aliquot that much, then add the DNA.
struct MasterMix {
    std::vector<Component> per_reaction;  // buffer, enzymes, then water last
    double reactions = 0;                 // replicates with the overage
    double aliquot_ul = 0;                // mix per tube
};

static MasterMix master_mix(const Reaction &r, int replicates) {
    if (replicates < 1) { throw std::invalid_argument("--replicates must be at least 1"); }
    MasterMix mix;
    mix.per_reaction = r.parts;
    mix.per_reaction.push_back({"Nuclease-free water", r.water_ul, ""});
    mix.reactions = replicates * (1.0 + MASTER_MIX_OVERAGE);
    mix.aliquot_ul = r.volume_ul - r.dna_ul;
    return mix;
}

// "10 ul per 50 ul rxn", scaled.
static double stop_solution_ul(double volume_ul) {
    return volume_ul * STOP_UL_PER_50UL / 50.0;
}

// --------------------------------------------------------- shared settings ---

// What the enzymes in one tube have in common. A digest with two enzymes has one
// buffer and one temperature, so if they disagree there is no protocol to write,
// only a choice for the person to make.
enum class Agreement { unknown, agree, disagree };

struct Common {
    Agreement how = Agreement::unknown;
    std::string text;   // the shared value when they agree
};

static Common common_buffer(const std::vector<Enzyme> &enzymes) {
    Common c;
    bool all_known = true;
    bool differ = false;
    for (const Enzyme &e : enzymes) {
        if (e.buffer.empty()) { all_known = false; continue; }
        if (c.text.empty()) { c.text = e.buffer; }
        else if (upper_of(c.text) != upper_of(e.buffer)) { differ = true; }
    }
    c.how = differ ? Agreement::disagree : (all_known ? Agreement::agree : Agreement::unknown);
    return c;
}

static Common common_temperature(const std::vector<Enzyme> &enzymes) {
    Common c;
    bool all_known = true;
    bool differ = false;
    int seen = 0;
    for (const Enzyme &e : enzymes) {
        if (e.incubation_c <= 0) { all_known = false; continue; }
        if (seen == 0) { seen = e.incubation_c; }
        else if (seen != e.incubation_c) { differ = true; }
    }
    if (seen > 0) { c.text = std::to_string(seen); }
    c.how = differ ? Agreement::disagree : (all_known ? Agreement::agree : Agreement::unknown);
    return c;
}

static std::string names_of(const std::vector<Enzyme> &enzymes, const std::string &sep,
                            const std::string &value_of_kind) {
    std::string out;
    for (const Enzyme &e : enzymes) {
        std::string value = value_of_kind == "buffer" ? e.buffer
                            : (e.incubation_c > 0 ? std::to_string(e.incubation_c) + " C" : "");
        out += (out.empty() ? "" : sep) + e.name + (value.empty() ? "" : " " + value);
    }
    return out;
}

// What to do about the enzyme once the digest is done. Heat inactivation is per
// enzyme, so this only says "heat" when every enzyme agrees on how. A source
// can also know the temperature but not the time (NEB's NEBuffer Performance
// Chart is exactly this), which gets its own bucket rather than a fabricated
// time or being folded into "unknown".
static std::vector<std::string> inactivation_lines(const std::vector<Enzyme> &enzymes) {
    std::vector<std::string> out;
    std::vector<std::string> unknown, cannot, temp_only;
    std::vector<std::pair<int, int> > conditions;
    std::vector<std::string> conditions_names;
    for (const Enzyme &e : enzymes) {
        if (e.inactivate_c < 0) { unknown.push_back(e.name); }
        else if (e.inactivate_c == 0) { cannot.push_back(e.name); }
        else if (e.inactivate_min <= 0) {
            temp_only.push_back(e.name + " at " + std::to_string(e.inactivate_c) + " C");
        } else {
            conditions.push_back({e.inactivate_c, e.inactivate_min});
            conditions_names.push_back(e.name);
        }
    }
    auto joined = [&](const std::vector<std::string> &names) {
        std::string s;
        for (const std::string &n : names) { s += (s.empty() ? "" : ", ") + n; }
        return s;
    };
    if (!unknown.empty()) {
        out.push_back("Heat inactivation for " + joined(unknown) + " is not on file; NEB lists it "
                      "per enzyme");
    }
    if (!cannot.empty()) {
        out.push_back(joined(cannot) + " cannot be heat inactivated (per NEB); clean up with a "
                      "spin column or phenol/chloroform extraction instead");
    }
    if (!temp_only.empty()) {
        out.push_back("Heat inactivate " + joined(temp_only) + "; the source on file gives the "
                      "temperature but not the time");
    }
    if (!conditions.empty()) {
        bool same = true;
        for (const auto &c : conditions) { same = same && c == conditions[0]; }
        if (same && unknown.empty() && cannot.empty() && temp_only.empty()) {
            out.push_back("Heat inactivate at " + std::to_string(conditions[0].first) + " C for " +
                          std::to_string(conditions[0].second) + " min");
        } else {
            for (size_t i = 0; i < conditions.size(); i++) {
                out.push_back(conditions_names[i] + ": heat inactivate at " +
                              std::to_string(conditions[i].first) + " C for " +
                              std::to_string(conditions[i].second) + " min");
            }
        }
    }
    return out;
}

// ------------------------------------------------------------------- main ---

static void print_enzyme_list(const std::vector<Conditions> &table) {
    std::cout << "Enzymes (* = conditions on file):" << std::endl;
    std::vector<std::string> shown;
    auto on_file = [&](const std::string &name) {
        for (const Conditions &c : table) {
            if (upper_of(c.name) == upper_of(name)) { return true; }
        }
        return false;
    };
    for (const SiteRow &s : enzyme_sites()) {
        std::cout << "  " << (on_file(s.name) ? "* " : "  ") << pad(s.name, 14) << s.site
                  << std::endl;
        shown.push_back(upper_of(s.name));
    }
    // Rows whose enzyme the site table lacks, such as the -HF versions.
    for (const Conditions &c : table) {
        if (std::find(shown.begin(), shown.end(), upper_of(c.name)) == shown.end()) {
            std::cout << "  * " << c.name << std::endl;
        }
    }
    std::cout << enzyme_sites().size() << " in the enzyme table, " << table.size()
              << " with conditions on file." << std::endl;
}

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("Digest-protocol", VERSION);
    program.add_description(
        "Restriction digest protocol generator - a reaction setup, incubation and stop for one or "
        "more enzymes. Amounts from NEB's Restriction Digest protocol; buffer, temperature and "
        "heat inactivation per enzyme from the conditions table or the command line.");

    program.add_argument("--enzymes", "-e")
        .help("Comma separated enzymes: NAME, or NAME:BUFFER:TEMP_C[:INACTIVATION] to give its "
              "conditions (INACTIVATION is TEMP_C/MINUTES, e.g. 65/20, or no)");
    program.add_argument("--dna-ug", "-d")
        .help("DNA in ug per reaction (default: the most NEB pairs with the volume, 1 ug at 50 uL)")
        .scan<'g', double>();
    program.add_argument("--dna-conc")
        .help("DNA concentration in ng/uL; works out the volume of DNA per tube")
        .scan<'g', double>();
    program.add_argument("--dna-volume")
        .help("DNA per tube in uL, left out of the master mix; sets its water (default: 1)")
        .default_value(1.0)
        .scan<'g', double>();
    program.add_argument("--volume")
        .help("Reaction volume in uL (default: 50)")
        .default_value(50.0)
        .scan<'g', double>();
    program.add_argument("--replicates", "-n")
        .help("Number of reactions; above 1, a master mix is made up with 10% extra (default: 1)")
        .default_value(1)
        .scan<'i', int>();
    program.add_argument("--genomic")
        .help("Genomic DNA: 20 units per ug instead of 10")
        .flag();
    program.add_argument("--units-per-ug")
        .help("Units of each enzyme per ug of DNA, instead of NEB's 10 (20 with --genomic)")
        .scan<'g', double>();
    program.add_argument("--stock")
        .help("An enzyme's concentration in U/uL, as NAME=UNITS; repeat for more. Read it off the "
              "tube: without it the tool assumes NEB's '10 units, generally 1 uL'")
        .append();
    program.add_argument("--minutes", "-t")
        .help("Incubation time in minutes (default: 60)")
        .default_value(DEFAULT_MINUTES)
        .scan<'i', int>();
    program.add_argument("--buffer")
        .help("The one 10X buffer for the whole reaction, when the enzymes' own do not settle it");
    program.add_argument("--incubate-c")
        .help("The one incubation temperature for the whole reaction")
        .scan<'i', int>();
    program.add_argument("--conditions-db")
        .metavar("CSV")
        .help("Read per-enzyme conditions from a CSV in the format of Data/digest_conditions.csv "
              "instead of the built-in table");
    program.add_argument("--list-enzymes")
        .help("List the known enzymes and which have conditions on file, and exit")
        .flag();

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 2;
    }

    std::vector<Conditions> loaded;
    const std::vector<Conditions> *table = NULL;
    try {
        if (program.is_used("--conditions-db")) {
            loaded = load_conditions_csv(program.get<std::string>("--conditions-db"));
            table = &loaded;
        } else {
            table = &builtin_conditions();
        }
    } catch (const std::exception &exc) {
        std::cerr << "Error reading conditions database: " << exc.what() << std::endl;
        return 1;
    }

    if (program.get<bool>("--list-enzymes")) {
        print_enzyme_list(*table);
        return 0;
    }
    if (!program.is_used("--enzymes")) {
        std::cerr << program;
        std::cerr << "Digest-protocol: error: the following arguments are required: "
                  << "--enzymes/-e" << std::endl;
        return 2;
    }

    // The enzymes, with whatever is known about each.
    std::vector<Enzyme> enzymes;
    try {
        for (const std::string &spec : split_on(program.get<std::string>("--enzymes"), ',')) {
            if (trimmed(spec).empty()) { continue; }
            Enzyme e = parse_enzyme_spec(spec, *table);
            for (const Enzyme &seen : enzymes) {
                if (upper_of(seen.name) == upper_of(e.name)) {
                    throw std::invalid_argument(e.name + " is listed twice");
                }
            }
            enzymes.push_back(e);
        }
        if (enzymes.empty()) { throw std::invalid_argument("no enzymes given"); }
        if (program.is_used("--stock")) {
            for (const std::string &item : program.get<std::vector<std::string> >("--stock")) {
                size_t eq = item.find('=');
                if (eq == std::string::npos) {
                    throw std::invalid_argument("--stock '" + item + "': write NAME=UNITS_PER_UL");
                }
                const std::string name = upper_of(trimmed(item.substr(0, eq)));
                const double stock = parse_number(trimmed(item.substr(eq + 1)), "stock",
                                                  "--stock '" + item + "': ");
                if (stock <= 0) { throw std::invalid_argument("--stock '" + item + "': must be positive"); }
                bool matched = false;
                for (Enzyme &e : enzymes) {
                    if (upper_of(e.name) == name) {
                        e.stock_u_per_ul = stock;
                        matched = true;
                    }
                }
                if (!matched) {
                    throw std::invalid_argument("--stock '" + item + "': that enzyme is not in "
                                                "--enzymes");
                }
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    const double volume_ul = program.get<double>("--volume");
    const int replicates = program.get<int>("--replicates");
    const int minutes = program.get<int>("--minutes");
    const bool genomic = program.get<bool>("--genomic");
    if (minutes <= 0) {
        std::cerr << "Error: --minutes must be positive" << std::endl;
        return 1;
    }

    // One buffer and one temperature for the tube. Enzymes that disagree, and no
    // override, is not something to guess a way out of.
    std::string buffer_name;
    std::string temperature_note;
    int incubation_c = 0;
    std::vector<std::string> warnings;
    {
        Common b = common_buffer(enzymes);
        Common t = common_temperature(enzymes);
        if (program.is_used("--buffer")) {
            buffer_name = program.get<std::string>("--buffer");
        } else if (b.how == Agreement::agree) {
            buffer_name = b.text;
        } else if (b.how == Agreement::disagree) {
            std::cerr << "Error: the enzymes want different buffers (" << names_of(enzymes, ", ", "buffer")
                      << "). A digest with them together needs one buffer they all work in, so "
                      << "either digest them one after the other, or give --buffer if you have "
                      << "checked one. NEB's double digest guidance covers choosing it."
                      << std::endl;
            return 1;
        }
        if (program.is_used("--incubate-c")) {
            incubation_c = program.get<int>("--incubate-c");
            if (incubation_c <= 0) {
                std::cerr << "Error: --incubate-c must be positive" << std::endl;
                return 1;
            }
        } else if (t.how == Agreement::agree) {
            incubation_c = std::atoi(t.text.c_str());
        } else if (t.how == Agreement::disagree) {
            std::cerr << "Error: the enzymes want different temperatures ("
                      << names_of(enzymes, ", ", "temperature") << "). Digest them one after the "
                      << "other, or give --incubate-c if you have checked that they work at one "
                      << "temperature." << std::endl;
            return 1;
        }
        if (buffer_name.empty()) {
            warnings.push_back("no buffer on file for " +
                               std::string(enzymes.size() > 1 ? "every enzyme" : enzymes[0].name) +
                               ": use the 10X NEBuffer NEB's chart lists, or give --buffer or "
                               "NAME:BUFFER:TEMP_C");
        }
        if (incubation_c == 0) {
            warnings.push_back("no incubation temperature on file for " +
                               std::string(enzymes.size() > 1 ? "every enzyme" : enzymes[0].name) +
                               ": it depends on the enzyme, and NEB's chart lists it; give "
                               "--incubate-c or NAME:BUFFER:TEMP_C");
        }
    }

    // The DNA: how much, and how much of the tube it takes.
    double dna_ug = 0;
    double units_per_ug = genomic ? GENOMIC_UNITS_PER_UG : UNITS_PER_UG;
    double dna_ul = program.get<double>("--dna-volume");
    const bool dna_volume_known = program.is_used("--dna-conc");
    try {
        if (program.is_used("--dna-ug")) {
            dna_ug = program.get<double>("--dna-ug");
        } else {
            dna_ug = max_dna_ug(volume_ul);
            if (dna_ug <= 0) {
                throw std::invalid_argument(
                    "a " + volume(volume_ul) + " uL reaction is smaller than the " +
                    volume(DNA_LIMITS[0].volume_ul) + " uL NEB tabulates; give --dna-ug");
            }
        }
        if (dna_ug <= 0) { throw std::invalid_argument("--dna-ug must be positive"); }
        if (program.is_used("--units-per-ug")) { units_per_ug = program.get<double>("--units-per-ug"); }
        if (dna_volume_known) {
            const double conc = program.get<double>("--dna-conc");
            if (conc <= 0) { throw std::invalid_argument("--dna-conc must be positive"); }
            dna_ul = dna_ug * 1000.0 / conc;
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    Reaction rxn;
    MasterMix mix;
    try {
        rxn = plan_reaction(enzymes, buffer_name.empty() ? "NEBuffer" : buffer_name, volume_ul,
                            dna_ug, dna_ul, units_per_ug);
        mix = master_mix(rxn, replicates);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    // Everything worth saying that is not an error.
    const double dna_max = max_dna_ug(volume_ul);
    if (dna_max > 0 && dna_ug > dna_max + 1e-9) {
        warnings.push_back(fixed(dna_ug, 3) + " ug is more than the " + fixed(dna_max, 3) +
                           " ug NEB pairs with " + volume(volume_ul) + " uL; NEB recommends "
                           "50 uL per ug of DNA");
    }
    if (dna_max == 0 && program.is_used("--dna-ug")) {
        warnings.push_back(volume(volume_ul) + " uL is below the " +
                           volume(DNA_LIMITS[0].volume_ul) + " uL NEB tabulates, so there is no "
                           "NEB figure to check the DNA against");
    }
    if (units_per_ug < UNITS_PER_UG_MIN - 1e-9) {
        warnings.push_back(volume(units_per_ug) + " units per ug is under NEB's recommended " +
                           volume(UNITS_PER_UG_MIN) + "-" + volume(UNITS_PER_UG) +
                           "; NEB says fewer units can work for up to " +
                           std::to_string(EXTENDED_DIGEST_HOURS) + " hours");
    }
    if (volume_ul <= SMALL_REACTION_UL + 1e-9 && minutes > DEFAULT_MINUTES) {
        warnings.push_back("NEB says not to incubate a " + volume(SMALL_REACTION_UL) +
                           " uL reaction over " + std::to_string(DEFAULT_MINUTES) +
                           " min, to avoid evaporation");
    }
    if (minutes > EXTENDED_DIGEST_HOURS * 60) {
        warnings.push_back(std::to_string(minutes) + " min is longer than the " +
                           std::to_string(EXTENDED_DIGEST_HOURS) + " hours NEB mentions");
    }
    for (const Enzyme &e : enzymes) {
        if (e.stock_u_per_ul <= 0) {
            warnings.push_back("stock concentration of " + e.name + " not on file, so " +
                               volume(DEFAULT_STOCK_U_PER_UL) + " U/uL is assumed (NEB's \"10 "
                               "units, generally 1 uL\"); if your tube says otherwise, give "
                               "--stock " + e.name + "=UNITS");
        }
    }
    for (const Enzyme &e : enzymes) {
        if (e.overrides_table) {
            warnings.push_back(e.name + "'s conditions on the command line differ from the "
                               "table's, so the table's source is not printed for it");
        }
    }

    std::string title;
    for (const Enzyme &e : enzymes) { title += (title.empty() ? "" : " + ") + e.name; }
    std::cout << "Restriction digest protocol: " << title << std::endl;

    std::cout << "\nEnzymes" << std::endl;
    for (const Enzyme &e : enzymes) {
        std::string what;
        if (e.buffer.empty() && e.incubation_c == 0) {
            what = "conditions not on file";
        } else {
            what = (e.buffer.empty() ? "buffer not on file" : "10X " + e.buffer) + ", " +
                   (e.incubation_c > 0 ? std::to_string(e.incubation_c) + " C" : "temperature not on file");
        }
        std::cout << "  " << pad(e.name, 14) << pad(e.site.empty() ? "-" : e.site, 18) << what;
        if (e.given) { std::cout << " (as given)"; }
        std::cout << std::endl;
    }

    // Reaction setup. The enzyme goes in last, and the water goes in first.
    auto row = [](const Component &c) {
        std::cout << "  " << pad(c.name, 36) << pad(volume(c.ul) + " uL", 12) << c.final
                  << std::endl;
    };
    const std::string dna_final = volume(dna_ug) + " ug (" + volume(dna_ug * 1000.0 / volume_ul) +
                                  " ng/uL)";
    std::cout << "\nReaction setup (" << volume(volume_ul) << " uL, " << volume(dna_ug)
              << " ug DNA, " << volume(rxn.units) << " units of each enzyme)" << std::endl;
    std::cout << "  " << pad("Component", 36) << pad("Volume", 12) << "Final" << std::endl;
    if (dna_volume_known) {
        std::cout << "  " << pad("Nuclease-free water", 36) << volume(rxn.water_ul) << " uL"
                  << std::endl;
    } else {
        std::cout << "  " << pad("Nuclease-free water", 36) << "to " << volume(volume_ul) << " uL"
                  << std::endl;
    }
    row(rxn.parts[0]);   // the buffer
    std::cout << "  " << pad("DNA", 36)
              << (dna_volume_known ? pad(volume(dna_ul) + " uL", 12) : pad("variable", 12))
              << dna_final << std::endl;
    for (size_t k = 1; k < rxn.parts.size(); k++) { row(rxn.parts[k]); }   // the enzymes, last
    std::cout << "  Enzyme(s) " << volume(rxn.enzyme_ul) << " uL of the " << volume(volume_ul)
              << " uL, " << fixed(100.0 * rxn.enzyme_ul / volume_ul, 1) << "% (NEB: at most "
              << fixed(ENZYME_MAX_FRACTION * 100, 0) << "%)." << std::endl;

    if (replicates > 1) {
        std::cout << "\nMaster mix (" << replicates << " reactions + "
                  << fixed(MASTER_MIX_OVERAGE * 100, 0) << "% = " << volume(mix.reactions) << ", "
                  << volume(mix.aliquot_ul * mix.reactions) << " uL)" << std::endl;
        std::cout << "  " << pad("Component", 36) << pad("Per rxn", 12) << "Mix" << std::endl;
        for (const Component &c : mix.per_reaction) {
            std::cout << "  " << pad(c.name, 36) << pad(volume(c.ul) + " uL", 12)
                      << volume(c.ul * mix.reactions) << " uL" << std::endl;
        }
        std::cout << "  Aliquot " << volume(mix.aliquot_ul) << " uL per tube, then add "
                  << volume(dna_ul) << " uL DNA." << std::endl;
    }

    std::cout << "\nProcedure" << std::endl;
    std::cout << "  1. Add the water, buffer and DNA, and the enzyme last. Keep the enzyme on ice "
              << "when it is not in the freezer." << std::endl;
    std::cout << "  2. Mix by pipetting up and down, or by flicking the tube. Give it a quick "
              << "spin in a microcentrifuge. Do not vortex." << std::endl;
    std::cout << "  3. Incubate " << (incubation_c > 0 ? "at " + std::to_string(incubation_c) + " C"
                                                        : std::string("at the enzyme's own temperature"))
              << " for " << minutes << " min." << std::endl;
    std::cout << "  4. Stop the reaction: with " << volume(stop_solution_ul(volume_ul))
              << " uL of stop solution (NEB uses 10 uL per 50 uL reaction) if the DNA needs no "
              << "further handling; if it does, heat inactivate or remove the enzyme with a spin "
              << "column or phenol/chloroform extraction." << std::endl;
    for (const std::string &line : inactivation_lines(enzymes)) {
        std::cout << "     " << line << "." << std::endl;
    }

    std::cout << "\nNotes" << std::endl;
    std::cout << "  Amounts: NEB's Restriction Digest protocol - buffer at 1X, " << volume(UNITS_PER_UG)
              << " units per ug of DNA (" << volume(GENOMIC_UNITS_PER_UG) << " for genomic), enzyme "
              << "at most 10% of the volume. Source: New England Biolabs, Restriction Digest v2, "
              << "protocols.io, DOI 10.17504/protocols.io.isycefw." << std::endl;
    std::cout << "  DNA should be free of phenol, chloroform, alcohol, EDTA, detergents and excess "
              << "salt; methylation can block some enzymes." << std::endl;
    std::string time_savers;
    for (const Enzyme &e : enzymes) {
        if (e.time_saver == 1) { time_savers += (time_savers.empty() ? "" : ", ") + e.name; }
    }
    if (!time_savers.empty()) {
        std::cout << "  Time-Saver Qualified (" << time_savers << "): NEB says "
                  << TIME_SAVER_MIN_MINUTES << "-" << TIME_SAVER_MAX_MINUTES
                  << " min is enough; the incubation here is " << minutes << " min." << std::endl;
    }
    if (rxn.enzyme_ul < 1.0) {
        std::cout << "  An enzyme volume under 1 uL is hard to pipette; NEB says enzymes can be "
                  << "diluted in their recommended diluent buffer when smaller amounts are needed."
                  << std::endl;
    }
    for (const Enzyme &e : enzymes) {
        if (!e.notes.empty()) { std::cout << "  " << e.name << ": " << e.notes << "." << std::endl; }
        if (!e.source.empty() && !e.overrides_table) {
            std::cout << "  " << e.name << " source: " << e.source << std::endl;
        }
    }

    if (!warnings.empty()) {
        std::cout << "\nWarnings" << std::endl;
        for (const std::string &w : warnings) { std::cout << "  ! " << w << std::endl; }
    }
    return 0;
}
