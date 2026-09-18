// PCR protocol generator - turns a primer pair and a polymerase into a reaction
// setup and a thermocycler program.
//
// Copyright (c) 2026 The Zabel Lab at Colorado State University.
//
// This file is part of Genomics_Scripts_Zlab and is distributed under the
// GPL v3 - see the LICENSE file at the root of the repository.
//
// What it computes, and where each number comes from:
//
//   Primer Tm   SantaLucia 1998 unified nearest-neighbour model with its salt
//               correction, magnesium folded in as sodium equivalents. The
//               formulas are the ones primer3 uses by default, so the two agree
//               under the same conditions (pinned by a test).
//   Annealing   The polymerase's own rule from its datasheet (Data/polymerases.csv),
//               applied to the lower primer Tm. The datasheets pitch their rules
//               against NEB's Tm calculator, whose model is not published, and
//               the two disagree: for NEB's own Phusion control primers NEB
//               quotes Tms 0.6-6.2 C above ours (see the README). Treat the
//               annealing temperature as the middle of a gradient, not gospel.
//   Extension   Amplicon length times the datasheet rate, rounded up to a whole
//               second. The length comes from locating both primers on a
//               template (-f), or is given directly (-l).
//   Primer QC   Hairpins and dimers from primer3's thal, vendored in
//               External_tools/primer3, under the same reaction conditions as
//               the Tm, with primer3's default limits. GC, length, runs and
//               the GC clamp against the datasheet and Premier Biosoft limits.
//               See pcr_data.hpp for every source.

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

#include "pcr_data.hpp"
#include "polymerase_csv_embedded.hpp"

// primer3 is C, and thal.h (unlike its other headers) has no extern "C" guard.
extern "C" {
#include "thal.h"
}
#include "thal_parameters.h"

using namespace pcr_data;

// Single source of truth for the version.
static const char *const VERSION = "0.1.0";

// ------------------------------------------------------------ sequences ---

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

static char complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default:  return 'N';
    }
}

static std::string reverse_complement(const std::string &seq) {
    std::string out(seq.rbegin(), seq.rend());
    for (char &c : out) { c = complement(c); }
    return out;
}

// A primer as typed: upper-cased, whitespace dropped, A/C/G/T only. Degenerate
// primers are refused rather than guessed at - the nearest-neighbour model has
// no parameters for an ambiguous stack.
static std::string validate_primer(const std::string &raw, const std::string &label) {
    std::string seq;
    for (char c : raw) {
        if (!isspace((unsigned char) c)) { seq += (char) toupper((unsigned char) c); }
    }
    if (seq.size() < 2) {
        throw std::invalid_argument(label + " primer is too short to have a Tm");
    }
    for (char c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            throw std::invalid_argument(
                label + " primer contains '" + std::string(1, c) +
                "'; only A, C, G and T are supported (no degenerate bases)");
        }
    }
    return seq;
}

static double gc_percent(const std::string &seq) {
    int gc = 0;
    for (char c : seq) {
        if (c == 'G' || c == 'C') { gc++; }
    }
    return 100.0 * gc / (double) seq.size();
}

// ---------------------------------------------------------- thermodynamics ---

struct Thermo {
    double dH = 0.0;  // kcal/mol
    double dS = 0.0;  // cal/(K mol)
};

// The stack 5'-ab-3' paired with its complement. The table holds each duplex
// once, so a stack not listed is looked up as the same duplex read from the
// other strand: 5'-TT-3' is 5'-AA-3'/3'-TT-5' upside down.
static Thermo nn_stack(char a, char b) {
    const std::string forward{a, b};
    const std::string flipped = reverse_complement(forward);
    for (const nn_row &row : NN_UNIFIED) {
        if (forward == row.stack || flipped == row.stack) {
            Thermo t;
            t.dH = row.dH;
            t.dS = row.dS;
            return t;
        }
    }
    throw std::invalid_argument("no nearest-neighbour parameters for " + forward);
}

static Thermo add(Thermo a, const Thermo &b) {
    a.dH += b.dH;
    a.dS += b.dS;
    return a;
}

// Initiation is charged once per terminal pair, by its type.
static Thermo initiation(char terminal) {
    Thermo t;
    bool gc = (terminal == 'G' || terminal == 'C');
    t.dH = gc ? INIT_GC_DH : INIT_AT_DH;
    t.dS = gc ? INIT_GC_DS : INIT_AT_DS;
    return t;
}

// Stacks only, no initiation: seq[from .. from + len) against its complement.
static Thermo stacks_of(const std::string &seq, size_t from, size_t len) {
    Thermo t;
    for (size_t k = from; k + 1 < from + len; k++) {
        t = add(t, nn_stack(seq[k], seq[k + 1]));
    }
    return t;
}

static bool is_self_complementary(const std::string &seq) {
    return seq == reverse_complement(seq);
}

// Full duplex of seq with its perfect complement, in 1 M NaCl.
static Thermo duplex_thermo(const std::string &seq) {
    Thermo t = stacks_of(seq, 0, seq.size());
    t = add(t, initiation(seq.front()));
    t = add(t, initiation(seq.back()));
    if (is_self_complementary(seq)) { t.dS += SYMMETRY_DS; }
    return t;
}

// Monovalent plus free magnesium as sodium equivalents, all in mM.
static double sodium_equivalent_mm(double monovalent_mm, double mg_mm, double dntp_total_mm) {
    double free_mg = mg_mm - dntp_total_mm;
    if (free_mg < 0.0) { free_mg = 0.0; }
    return monovalent_mm + MG_TO_NA_FACTOR * std::sqrt(free_mg);
}

static double salt_corrected_ds(double dS, int stacks, double sodium_mm) {
    return dS + SALT_DS_PER_STACK * stacks * std::log(sodium_mm / 1000.0);
}

// Tm in C of the primer against its complement. primer_um is the primer
// concentration in the reaction; like primer3, the strand concentration term is
// C/4 for a non self-complementary duplex and C for a self-complementary one.
static double melting_temp(const std::string &seq, double sodium_mm, double primer_um) {
    Thermo t = duplex_thermo(seq);
    double dS = salt_corrected_ds(t.dS, (int) seq.size() - 1, sodium_mm);
    double molar = primer_um * 1e-6;
    double strands = is_self_complementary(seq) ? molar : molar / 4.0;
    return 1000.0 * t.dH / (dS + GAS_CONSTANT * std::log(strands)) - KELVIN;
}

// ----------------------------------------------------- secondary structure ---

// thal keeps its parameters in static tables, loaded once from the copies
// compiled into thal_parameters.c - no parameter files at run time.
static void load_thal_once() {
    static bool loaded = false;
    if (loaded) { return; }
    thal_parameters params;
    thal_results result;
    thal_set_null_parameters(&params);
    set_default_thal_parameters(&params);
    int failed = get_thermodynamic_values(&params, &result);
    thal_free_parameters(&params);
    if (failed) {
        throw std::runtime_error(std::string("could not load primer3 thal parameters: ") +
                                 result.msg);
    }
    loaded = true;
}

// Melting temperature in C of the most stable structure thal finds, 0 when
// there is none. Both oligos 5'->3' as ordered, the way primer3 passes them.
// Concentrations are the reaction's, the same ones the primer Tm uses.
static double structure_tm(const std::string &a, const std::string &b, thal_alignment_type type,
                           double monovalent_mm, double mg_mm, double dntp_total_mm,
                           double primer_um) {
    load_thal_once();
    thal_args args;
    set_thal_default_args(&args);
    args.type = type;
    args.mv = monovalent_mm;
    args.dv = mg_mm;
    args.dntp = dntp_total_mm;
    args.dna_conc = primer_um * 1000.0;  // nM
    if (type == thal_hairpin) { args.dimer = 0; }
    thal_results result;
    thal((const unsigned char *) a.c_str(), (const unsigned char *) b.c_str(), &args, THL_FAST,
         &result);
    if (result.temp == THAL_ERROR_SCORE) {
        throw std::runtime_error(std::string("primer3 thal failed: ") + result.msg);
    }
    return result.temp;
}

// -------------------------------------------------------------- primer QC ---

// Longest run of one base, and most consecutive copies of one dinucleotide.
static int longest_run(const std::string &s) {
    int best = s.empty() ? 0 : 1, run = 1;
    for (size_t k = 1; k < s.size(); k++) {
        run = (s[k] == s[k - 1]) ? run + 1 : 1;
        best = std::max(best, run);
    }
    return best;
}

static int most_dinucleotide_repeats(const std::string &s) {
    int best = 0;
    for (size_t start = 0; start + 1 < s.size(); start++) {
        if (s[start] == s[start + 1]) { continue; }  // a run, not a dinucleotide
        int copies = 1;
        size_t k = start + 2;
        while (k + 1 < s.size() && s[k] == s[start] && s[k + 1] == s[start + 1]) {
            copies++;
            k += 2;
        }
        best = std::max(best, copies);
    }
    return best;
}

static int gc_in_last(const std::string &s, int window) {
    int gc = 0;
    size_t from = s.size() > (size_t) window ? s.size() - window : 0;
    for (size_t k = from; k < s.size(); k++) {
        if (s[k] == 'G' || s[k] == 'C') { gc++; }
    }
    return gc;
}

// -------------------------------------------------------------- polymerases ---

struct Polymerase {
    std::string name, full_name, catalog, buffer, notes, source;
    double buffer_x = 0, mg_mm = 0, monovalent_mm = 0, dntp_each_um = 0, primer_um = 0;
    double enzyme_units_per_50ul = 0, enzyme_stock_u_per_ul = 0;
    int initial_denature_c = 0, initial_denature_s = 0, denature_c = 0, denature_s = 0;
    int anneal_offset_c = 0, anneal_offset_min_primer_nt = 0, anneal_min_c = 0, anneal_s = 0;
    int two_step_ta_c = 0, extend_c = 0;
    int extend_s_per_kb = 0, extend_s_per_kb_simple = 0, long_amplicon_bp = 0;
    int extend_s_per_kb_long = 0, extend_min_s = 0;
    int final_extend_c = 0, final_extend_s = 0, hold_c = 0, cycles = 0;
};

// Plain comma split: the CSV has no quoting, and notes use '-' and ';' instead
// of commas so that stays true.
static std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> fields;
    std::string current;
    for (char c : line) {
        if (c == ',') { fields.push_back(current); current.clear(); }
        else if (c != '\r') { current += c; }
    }
    fields.push_back(current);
    return fields;
}

static const char *const POLYMERASE_COLUMNS[] = {
    "polymerase", "full_name", "catalog", "buffer", "buffer_x", "mg_mm", "monovalent_mm",
    "dntp_each_um", "primer_um", "enzyme_units_per_50ul", "enzyme_stock_u_per_ul",
    "initial_denature_c", "initial_denature_s", "denature_c", "denature_s",
    "anneal_offset_c", "anneal_offset_min_primer_nt", "anneal_min_c", "anneal_s",
    "two_step_ta_c", "extend_c", "extend_s_per_kb", "extend_s_per_kb_simple",
    "long_amplicon_bp", "extend_s_per_kb_long", "extend_min_s", "final_extend_c",
    "final_extend_s", "hold_c", "cycles", "notes", "source",
};

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

// Read a polymerase CSV with the columns above, in any order. Every column is
// required: a missing one would otherwise quietly become a zero-second step.
static std::vector<Polymerase> load_polymerase_csv_string(const std::string &text,
                                                          const std::string &path) {
    std::istringstream stream(text);
    std::string line;
    if (!std::getline(stream, line)) {
        throw std::invalid_argument("No polymerases found in " + path);
    }
    std::vector<std::string> header = split_csv_line(line);
    std::vector<int> col;
    std::string missing;
    for (const char *name : POLYMERASE_COLUMNS) {
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

    std::vector<Polymerase> polymerases;
    int lineno = 1;
    while (std::getline(stream, line)) {
        lineno++;
        if (trimmed(line).empty()) { continue; }
        std::vector<std::string> row = split_csv_line(line);
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
            const char *name = POLYMERASE_COLUMNS[c];
            return parse_number(text(), name, where);
        };
        auto integer = [&](void) -> int {
            const char *name = POLYMERASE_COLUMNS[c];
            return parse_integer(text(), name, where);
        };
        // Same order as POLYMERASE_COLUMNS.
        Polymerase p;
        p.name = text();
        p.full_name = text();
        p.catalog = text();
        p.buffer = text();
        p.buffer_x = number();
        p.mg_mm = number();
        p.monovalent_mm = number();
        p.dntp_each_um = number();
        p.primer_um = number();
        p.enzyme_units_per_50ul = number();
        p.enzyme_stock_u_per_ul = number();
        p.initial_denature_c = integer();
        p.initial_denature_s = integer();
        p.denature_c = integer();
        p.denature_s = integer();
        p.anneal_offset_c = integer();
        p.anneal_offset_min_primer_nt = integer();
        p.anneal_min_c = integer();
        p.anneal_s = integer();
        p.two_step_ta_c = integer();
        p.extend_c = integer();
        p.extend_s_per_kb = integer();
        p.extend_s_per_kb_simple = integer();
        p.long_amplicon_bp = integer();
        p.extend_s_per_kb_long = integer();
        p.extend_min_s = integer();
        p.final_extend_c = integer();
        p.final_extend_s = integer();
        p.hold_c = integer();
        p.cycles = integer();
        p.notes = text();
        p.source = text();
        if (p.name.empty()) {
            throw std::invalid_argument(where + "polymerase name is empty");
        }
        if (p.buffer_x <= 0 || p.enzyme_stock_u_per_ul <= 0 || p.primer_um <= 0) {
            throw std::invalid_argument(
                where + "buffer_x, enzyme_stock_u_per_ul and primer_um must be positive");
        }
        polymerases.push_back(p);
    }
    if (polymerases.empty()) {
        throw std::invalid_argument("No polymerases found in " + path);
    }
    return polymerases;
}

static std::vector<Polymerase> load_polymerase_csv(const std::string &path) {
    std::ifstream handle(path.c_str());
    if (!handle) {
        throw std::runtime_error("No such file or directory: '" + path + "'");
    }
    std::ostringstream buffer;
    buffer << handle.rdbuf();
    return load_polymerase_csv_string(buffer.str(), path);
}

static const std::vector<Polymerase> &builtin_polymerases() {
    static std::vector<Polymerase> table =
        load_polymerase_csv_string(EMBEDDED_POLYMERASE_CSV, "polymerases.csv");
    return table;
}

// Case-insensitive, by the short name.
static const Polymerase &find_polymerase(const std::string &name,
                                         const std::vector<Polymerase> &table) {
    const std::string wanted = upper_of(trimmed(name));
    for (const Polymerase &p : table) {
        if (upper_of(p.name) == wanted) { return p; }
    }
    std::string known;
    for (const Polymerase &p : table) { known += (known.empty() ? "" : ", ") + p.name; }
    throw std::invalid_argument("unknown polymerase '" + name + "'; known: " + known);
}

// Sodium equivalents for this polymerase's reaction buffer. dNTPs are listed
// per nucleotide, so the total that chelates Mg is four times that.
static double reaction_sodium_mm(const Polymerase &p) {
    return sodium_equivalent_mm(p.monovalent_mm, p.mg_mm, 4.0 * p.dntp_each_um / 1000.0);
}

// ---------------------------------------------------------------- cycling ---

struct Annealing {
    int ta = 0;             // C, whole degrees
    bool two_step = false;  // annealing merged into extension
    bool two_step_possible = false;
    bool below_range = false;
};

// The datasheet rule, applied to the lower Tm: offset it (only when the
// shorter primer is long enough, if the polymerase says so), round to a whole
// degree. At or above the extension temperature there is nothing to gain from
// a separate annealing step, so it is folded into extension. Above
// two_step_ta_c the datasheet merely allows a 2-step program (Taq: "above
// 65 C"), which is left as a note.
static Annealing annealing_for(const Polymerase &p, double tm_forward, double tm_reverse,
                               int len_forward, int len_reverse) {
    double tm_low = std::min(tm_forward, tm_reverse);
    int shorter = std::min(len_forward, len_reverse);
    int offset = (shorter >= p.anneal_offset_min_primer_nt) ? p.anneal_offset_c : 0;
    Annealing a;
    a.ta = (int) std::floor(tm_low + offset + 0.5);
    if (a.ta >= p.extend_c) {
        a.ta = p.extend_c;
        a.two_step = true;
    } else if (a.ta > p.two_step_ta_c) {
        a.two_step_possible = true;
    }
    a.below_range = a.ta < p.anneal_min_c;
    return a;
}

// Per-cycle extension in whole seconds, rounded up. Integer arithmetic on
// purpose: 1347 bp at 30 s/kb is 40.41 s and must come out as 41, not as
// whatever 1.347 * 30 happens to round to.
static int extension_seconds(const Polymerase &p, long amplicon_bp, bool simple_template) {
    int rate = p.extend_s_per_kb;
    if (simple_template) {
        rate = p.extend_s_per_kb_simple;
    } else if (p.long_amplicon_bp > 0 && amplicon_bp > p.long_amplicon_bp) {
        rate = p.extend_s_per_kb_long;
    }
    long seconds = (amplicon_bp * rate + 999) / 1000;
    if (seconds < p.extend_min_s) { seconds = p.extend_min_s; }
    if (seconds < 1) { seconds = 1; }
    return (int) seconds;
}

// ---------------------------------------------------------------- template ---

static std::vector<std::pair<std::string, std::string> > parse_fasta(const std::string &path) {
    std::ifstream handle(path.c_str());
    if (!handle) {
        throw std::runtime_error("File not found: " + path);
    }
    std::vector<std::pair<std::string, std::string> > records;
    std::string line;
    while (std::getline(handle, line)) {
        std::string text = trimmed(line);
        if (text.empty()) { continue; }
        if (text[0] == '>') {
            records.push_back(std::make_pair(trimmed(text.substr(1)), std::string()));
        } else {
            if (records.empty()) {
                throw std::runtime_error("FASTA file missing header line: " + path);
            }
            records.back().second += upper_of(text);
        }
    }
    if (records.empty()) {
        throw std::runtime_error("No sequences found in FASTA file: " + path);
    }
    return records;
}

static std::vector<size_t> find_all(const std::string &haystack, const std::string &needle) {
    std::vector<size_t> hits;
    for (size_t at = haystack.find(needle); at != std::string::npos;
         at = haystack.find(needle, at + 1)) {
        hits.push_back(at);
    }
    return hits;
}

struct Product {
    std::string record;
    long start = 0;       // 1-based, on the record as given
    long end = 0;         // 1-based, inclusive
    char strand = '+';    // '+' when the forward primer sits on the record's top strand
    long length() const { return end - start + 1; }
};

// Every product the pair would make on the template. Both primers are 5'->3'
// as ordered, so on the strand the forward primer reads along, the reverse
// primer appears as its reverse complement, downstream. The template may be
// either way round, so the mirror image is searched too.
static std::vector<Product> find_products(
        const std::vector<std::pair<std::string, std::string> > &records,
        const std::string &forward, const std::string &reverse) {
    std::vector<Product> products;
    const std::string forward_rc = reverse_complement(forward);
    const std::string reverse_rc = reverse_complement(reverse);
    for (const auto &record : records) {
        const std::string &seq = record.second;
        // (upstream primer, downstream primer's reverse complement, strand)
        const std::string *orientations[2][2] = {{&forward, &reverse_rc}, {&reverse, &forward_rc}};
        for (int o = 0; o < 2; o++) {
            const std::string &up = *orientations[o][0];
            const std::string &down = *orientations[o][1];
            for (size_t u : find_all(seq, up)) {
                for (size_t d : find_all(seq, down)) {
                    if (d < u) { continue; }
                    Product p;
                    p.record = record.first;
                    p.start = (long) u + 1;
                    p.end = (long) (d + down.size());
                    p.strand = (o == 0) ? '+' : '-';
                    if (p.end < (long) (u + up.size())) { continue; }
                    products.push_back(p);
                }
            }
        }
    }
    return products;
}

// Why a pair made nothing. The common slip is giving the reverse primer as the
// top-strand sequence instead of 5'->3' as ordered, so that is checked for by
// name.
static std::string explain_no_product(
        const std::vector<std::pair<std::string, std::string> > &records,
        const std::string &forward, const std::string &reverse) {
    bool fwd_top = false, fwd_bottom = false, rev_top = false, rev_bottom = false;
    for (const auto &record : records) {
        const std::string &seq = record.second;
        fwd_top |= seq.find(forward) != std::string::npos;
        fwd_bottom |= seq.find(reverse_complement(forward)) != std::string::npos;
        rev_top |= seq.find(reverse) != std::string::npos;
        rev_bottom |= seq.find(reverse_complement(reverse)) != std::string::npos;
    }
    if (!fwd_top && !fwd_bottom) { return "the forward primer does not match the template"; }
    if (!rev_top && !rev_bottom) { return "the reverse primer does not match the template"; }
    if ((fwd_top && rev_top && !rev_bottom) || (fwd_bottom && rev_bottom && !rev_top)) {
        return "both primers match the same strand - give the reverse primer 5'->3' as "
               "ordered, i.e. the reverse complement of the top strand";
    }
    return "the primers bind but face away from each other";
}

// ------------------------------------------------------------------ output ---

static std::string fixed(double value, int places) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(places) << value;
    return out.str();
}

// Volumes with at most three decimals and no trailing zeros: 2.5, 0.25, 10.
static std::string volume(double ul) {
    std::string s = fixed(ul, 3);
    s.erase(s.find_last_not_of('0') + 1);
    if (!s.empty() && s.back() == '.') { s.pop_back(); }
    return s;
}

static std::string duration(int seconds) {
    if (seconds < 60) { return std::to_string(seconds) + " s"; }
    std::string out = std::to_string(seconds / 60) + " min";
    if (seconds % 60) { out += " " + std::to_string(seconds % 60) + " s"; }
    return out;
}

static std::string pad(const std::string &s, size_t width) {
    return s.size() >= width ? s + " " : s + std::string(width - s.size(), ' ');
}

// A structure Tm, flagged when primer3 would reject it.
static std::string describe_structure(double tm, bool *bad) {
    if (tm <= 0.0) { return "none"; }
    bool over = tm > STRUCTURE_TM_LIMIT;
    *bad = *bad || over;
    return fixed(tm, 1) + " C" + (over ? " (!)" : "");
}

// Structure Tms as primer3 reports them. "3' end" is thal's END alignment,
// anchored on a 3' end: the structure the polymerase can extend.
struct PrimerReport {
    std::string label, seq;
    double tm = 0, gc = 0;
    double hairpin = 0, self_any = 0, self_end = 0;
};

static double reaction_structure_tm(const std::string &a, const std::string &b,
                                    thal_alignment_type type, const Polymerase &p) {
    return structure_tm(a, b, type, p.monovalent_mm, p.mg_mm, 4.0 * p.dntp_each_um / 1000.0,
                        p.primer_um);
}

static PrimerReport report_primer(const std::string &label, const std::string &seq,
                                  const Polymerase &p) {
    PrimerReport r;
    r.label = label;
    r.seq = seq;
    r.tm = melting_temp(seq, reaction_sodium_mm(p), p.primer_um);
    r.gc = gc_percent(seq);
    r.hairpin = reaction_structure_tm(seq, seq, thal_hairpin, p);
    r.self_any = reaction_structure_tm(seq, seq, thal_any, p);
    r.self_end = reaction_structure_tm(seq, seq, thal_end1, p);
    return r;
}

// The pair annealed to each other. END1 anchors the forward primer's 3' end
// and END2 the reverse primer's, so the worse of the two covers either one
// being extended.
static void report_pair(const std::string &forward, const std::string &reverse,
                        const Polymerase &p, double *any_out, double *end_out) {
    *any_out = reaction_structure_tm(forward, reverse, thal_any, p);
    *end_out = std::max(reaction_structure_tm(forward, reverse, thal_end1, p),
                        reaction_structure_tm(forward, reverse, thal_end2, p));
}

static void primer_warnings(const PrimerReport &r, std::vector<std::string> *warnings) {
    const int len = (int) r.seq.size();
    if (len < PRIMER_MIN_NT || len > PRIMER_MAX_NT) {
        warnings->push_back(r.label + " primer is " + std::to_string(len) + " nt; primers are "
                            "generally " + std::to_string(PRIMER_MIN_NT) + "-" +
                            std::to_string(PRIMER_MAX_NT) + " nt");
    }
    if (r.gc < GC_MIN_PERCENT || r.gc > GC_MAX_PERCENT) {
        warnings->push_back(r.label + " primer GC is " + fixed(r.gc, 1) + "%, outside " +
                            fixed(GC_MIN_PERCENT, 0) + "-" + fixed(GC_MAX_PERCENT, 0) + "%");
    }
    int clamp = gc_in_last(r.seq, GC_CLAMP_WINDOW);
    if (clamp > GC_CLAMP_MAX) {
        warnings->push_back(r.label + " primer has " + std::to_string(clamp) + " G/C in its last " +
                            std::to_string(GC_CLAMP_WINDOW) + " bases (more than " +
                            std::to_string(GC_CLAMP_MAX) + " invites mispriming)");
    }
    int run = longest_run(r.seq);
    if (run > MAX_RUN) {
        warnings->push_back(r.label + " primer has a run of " + std::to_string(run) +
                            " identical bases");
    }
    int repeats = most_dinucleotide_repeats(r.seq);
    if (repeats > MAX_DINUCLEOTIDE_REPEATS) {
        warnings->push_back(r.label + " primer has " + std::to_string(repeats) +
                            " dinucleotide repeats in a row");
    }
}

static void print_polymerase_list(const std::vector<Polymerase> &table) {
    std::cout << "Polymerases:" << std::endl;
    for (const Polymerase &p : table) {
        std::cout << "  " << pad(p.name, 10) << p.full_name << " (" << p.catalog << ")" << std::endl;
    }
}

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("PCR-protocol", VERSION);
    program.add_description(
        "PCR protocol generator - a reaction setup and thermocycler program for a primer pair "
        "and polymerase. Tm by SantaLucia 1998 nearest neighbour; annealing, extension and "
        "reaction setup from the polymerase's datasheet.");

    program.add_argument("--forward", "-F").help("Forward primer, 5'->3'");
    program.add_argument("--reverse", "-R").help("Reverse primer, 5'->3' as ordered");
    program.add_argument("--polymerase", "-p")
        .help("Polymerase by name, see --list-polymerases (case-insensitive)");
    program.add_argument("--fasta", "-f")
        .help("Template FASTA: both primers are located on it to get the amplicon length");
    program.add_argument("--amplicon-length", "-l")
        .help("Amplicon length in bp, instead of --fasta")
        .scan<'i', int>();
    program.add_argument("--cycles", "-c")
        .help("Number of cycles (default: the polymerase's, 30 for all built-ins)")
        .scan<'i', int>();
    program.add_argument("--volume")
        .help("Reaction volume in uL (default: 50)")
        .default_value(50.0)
        .scan<'g', double>();
    program.add_argument("--simple-template")
        .help("Plasmid, lambda or E. coli template: use the datasheet's faster extension rate")
        .flag();
    program.add_argument("--polymerase-db")
        .metavar("CSV")
        .help("Read polymerases from a CSV in the format of Data/polymerases.csv instead of "
              "the built-in table");
    program.add_argument("--list-polymerases")
        .help("List the known polymerases and exit")
        .flag();

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 2;
    }

    std::vector<Polymerase> loaded;
    const std::vector<Polymerase> *table = &builtin_polymerases();
    if (program.is_used("--polymerase-db")) {
        try {
            loaded = load_polymerase_csv(program.get<std::string>("--polymerase-db"));
        } catch (const std::exception &exc) {
            std::cerr << "Error reading polymerase database: " << exc.what() << std::endl;
            return 1;
        }
        table = &loaded;
    }

    if (program.get<bool>("--list-polymerases")) {
        print_polymerase_list(*table);
        return 0;
    }

    if (!program.is_used("--forward") || !program.is_used("--reverse") ||
        !program.is_used("--polymerase")) {
        std::cerr << program;
        std::cerr << "PCR-protocol: error: the following arguments are required: "
                  << "--forward/-F, --reverse/-R, --polymerase/-p" << std::endl;
        return 2;
    }
    if (program.is_used("--fasta") == program.is_used("--amplicon-length")) {
        std::cerr << "PCR-protocol: error: give exactly one of --fasta/-f or "
                  << "--amplicon-length/-l" << std::endl;
        return 2;
    }

    std::string forward, reverse;
    const Polymerase *pol = NULL;
    try {
        forward = validate_primer(program.get<std::string>("--forward"), "Forward");
        reverse = validate_primer(program.get<std::string>("--reverse"), "Reverse");
        pol = &find_polymerase(program.get<std::string>("--polymerase"), *table);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    // The amplicon, from the template or as given.
    long amplicon_bp = 0;
    std::string amplicon_note;
    if (program.is_used("--fasta")) {
        const std::string path = program.get<std::string>("--fasta");
        std::vector<std::pair<std::string, std::string> > records;
        try {
            records = parse_fasta(path);
        } catch (const std::exception &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
        std::vector<Product> products = find_products(records, forward, reverse);
        if (products.empty()) {
            std::cerr << "Error: no product on " << path << ": "
                      << explain_no_product(records, forward, reverse)
                      << ". Primers must match the template exactly; for tailed primers "
                      << "pass --amplicon-length instead." << std::endl;
            return 1;
        }
        if (products.size() > 1) {
            std::cerr << "Error: the pair makes " << products.size() << " products on " << path
                      << ", so the amplicon length is ambiguous:" << std::endl;
            for (const Product &p : products) {
                std::cerr << "  " << p.record << "  " << p.start << ".." << p.end << " ("
                          << p.strand << ")  " << p.length() << " bp" << std::endl;
            }
            return 1;
        }
        const Product &p = products.front();
        amplicon_bp = p.length();
        // The record's ID, not its whole description line.
        const std::string id = p.record.substr(0, p.record.find_first_of(" \t"));
        amplicon_note = id + ", bases " + std::to_string(p.start) + ".." +
                        std::to_string(p.end) +
                        (p.strand == '-' ? ", on the reverse strand" : "");
    } else {
        amplicon_bp = program.get<int>("--amplicon-length");
        if (amplicon_bp < (long) std::max(forward.size(), reverse.size())) {
            std::cerr << "Error: an amplicon of " << amplicon_bp << " bp is shorter than "
                      << "the primers" << std::endl;
            return 1;
        }
        amplicon_note = "as given";
    }

    const int cycles = program.is_used("--cycles") ? program.get<int>("--cycles") : pol->cycles;
    const double volume_ul = program.get<double>("--volume");
    const bool simple_template = program.get<bool>("--simple-template");
    if (cycles < 1 || volume_ul <= 0) {
        std::cerr << "Error: --cycles and --volume must be positive" << std::endl;
        return 1;
    }

    // Primers.
    const double sodium = reaction_sodium_mm(*pol);
    PrimerReport fwd, rev;
    double cross_any = 0, cross_end = 0;
    try {
        fwd = report_primer("Forward", forward, *pol);
        rev = report_primer("Reverse", reverse, *pol);
        report_pair(forward, reverse, *pol, &cross_any, &cross_end);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    Annealing anneal = annealing_for(*pol, fwd.tm, rev.tm, (int) forward.size(),
                                     (int) reverse.size());
    const int extend_s = extension_seconds(*pol, amplicon_bp, simple_template);

    std::vector<std::string> warnings;
    primer_warnings(fwd, &warnings);
    primer_warnings(rev, &warnings);
    if (std::fabs(fwd.tm - rev.tm) >= TM_DIFFERENCE_LIMIT) {
        warnings.push_back("primer Tms differ by " + fixed(std::fabs(fwd.tm - rev.tm), 1) +
                           " C; " + fixed(TM_DIFFERENCE_LIMIT, 0) +
                           " C or more can give no product");
    }
    if (anneal.below_range) {
        warnings.push_back("annealing at " + std::to_string(anneal.ta) + " C is below the " +
                           std::to_string(pol->anneal_min_c) + " C the " + pol->name +
                           " datasheet goes down to");
    }

    std::cout << "PCR protocol: " << pol->full_name << " (" << pol->catalog << ")" << std::endl;

    std::cout << "\nPrimers" << std::endl;
    for (const PrimerReport *r : {&fwd, &rev}) {
        std::cout << "  " << pad(r->label, 8) << "5'-" << r->seq << "-3'" << std::endl;
        std::cout << "  " << pad("", 8) << r->seq.size() << " nt, GC " << fixed(r->gc, 1)
                  << "%, Tm " << fixed(r->tm, 1) << " C" << std::endl;
    }
    std::cout << "  Amplicon " << amplicon_bp << " bp (" << amplicon_note << ")" << std::endl;

    bool structure_bad = false;
    std::cout << "\nPrimer structures (primer3 thal, Tm of the most stable; flagged above "
              << fixed(STRUCTURE_TM_LIMIT, 0) << " C)" << std::endl;
    std::cout << "  " << pad("", 19) << pad("3' end", 12) << "anywhere" << std::endl;
    for (const PrimerReport *r : {&fwd, &rev}) {
        std::cout << "  " << pad(r->label + " hairpin", 19) << pad("", 12)
                  << describe_structure(r->hairpin, &structure_bad) << std::endl;
        std::cout << "  " << pad(r->label + " self-dimer", 19)
                  << pad(describe_structure(r->self_end, &structure_bad), 12)
                  << describe_structure(r->self_any, &structure_bad) << std::endl;
    }
    std::cout << "  " << pad("Cross-dimer", 19) << pad(describe_structure(cross_end, &structure_bad), 12)
              << describe_structure(cross_any, &structure_bad) << std::endl;
    if (structure_bad) {
        warnings.push_back("a primer structure melts above " + fixed(STRUCTURE_TM_LIMIT, 0) +
                           " C, marked (!) above; primer3 would reject it");
    }

    // Reaction setup, scaled from the datasheet's 50 uL.
    const double scale = volume_ul / 50.0;
    const double units = pol->enzyme_units_per_50ul * scale;
    std::cout << "\nReaction setup (" << volume(volume_ul) << " uL)" << std::endl;
    std::cout << "  " << pad("Component", 36) << pad("Volume", 12) << "Final" << std::endl;
    std::cout << "  " << pad(pol->buffer, 36) << pad(volume(volume_ul / pol->buffer_x) + " uL", 12)
              << "1X" << std::endl;
    std::cout << "  " << pad("10 mM dNTPs", 36)
              << pad(volume(volume_ul * pol->dntp_each_um / 10000.0) + " uL", 12)
              << volume(pol->dntp_each_um) << " uM each" << std::endl;
    for (const char *which : {"Forward", "Reverse"}) {
        std::cout << "  " << pad(std::string("10 uM ") + which + " primer", 36)
                  << pad(volume(volume_ul * pol->primer_um / 10.0) + " uL", 12)
                  << volume(pol->primer_um) << " uM" << std::endl;
    }
    std::cout << "  " << pad("Template DNA", 36) << "variable" << std::endl;
    std::cout << "  "
              << pad(pol->name + " polymerase (" + volume(pol->enzyme_stock_u_per_ul) + " U/uL)", 36)
              << pad(volume(units / pol->enzyme_stock_u_per_ul) + " uL", 12) << volume(units)
              << " U" << std::endl;
    std::cout << "  " << pad("Nuclease-free water", 36) << "to " << volume(volume_ul) << " uL"
              << std::endl;

    // Thermocycler program.
    auto step = [](const std::string &name, int temp_c, const std::string &time) {
        std::cout << "  " << pad(name, 24) << pad(std::to_string(temp_c) + " C", 8) << time
                  << std::endl;
    };
    std::cout << "\nThermocycler program" << (anneal.two_step ? " (2-step)" : "") << std::endl;
    step("Initial denaturation", pol->initial_denature_c, duration(pol->initial_denature_s));
    std::cout << "  " << cycles << " cycles of:" << std::endl;
    step("  Denaturation", pol->denature_c, duration(pol->denature_s));
    if (anneal.two_step) {
        step("  Anneal + extension", pol->extend_c, duration(extend_s));
    } else {
        step("  Annealing", anneal.ta, duration(pol->anneal_s));
        step("  Extension", pol->extend_c, duration(extend_s));
    }
    step("Final extension", pol->final_extend_c, duration(pol->final_extend_s));
    step("Hold", pol->hold_c, "forever");

    std::cout << "\nNotes" << std::endl;
    std::cout << "  Annealing: lower Tm " << fixed(std::min(fwd.tm, rev.tm), 1) << " C";
    int shorter = (int) std::min(forward.size(), reverse.size());
    if (shorter >= pol->anneal_offset_min_primer_nt && pol->anneal_offset_c != 0) {
        std::cout << (pol->anneal_offset_c > 0 ? " + " : " - ") << std::abs(pol->anneal_offset_c)
                  << " C";
    }
    std::cout << " per the " << pol->name << " datasheet";
    if (anneal.two_step) {
        std::cout << ", which reaches the " << pol->extend_c
                  << " C extension temperature, so annealing is folded into extension";
    }
    std::cout << ". A gradient around it is the surest check." << std::endl;
    if (anneal.two_step_possible) {
        std::cout << "  At this annealing temperature the datasheet also allows a 2-step "
                  << "program (annealing folded into extension)." << std::endl;
    }
    std::cout << "  Tm: SantaLucia 1998 nearest neighbour, " << volume(pol->monovalent_mm)
              << " mM monovalent + " << volume(pol->mg_mm) << " mM Mg2+ ("
              << fixed(sodium, 0) << " mM Na+ equivalent), " << volume(pol->primer_um)
              << " uM primer." << std::endl;
    if (!pol->notes.empty()) { std::cout << "  Data notes: " << pol->notes << "." << std::endl; }
    std::cout << "  Source: " << pol->source << std::endl;

    if (!warnings.empty()) {
        std::cout << "\nWarnings" << std::endl;
        for (const std::string &w : warnings) { std::cout << "  ! " << w << std::endl; }
    }
    return 0;
}
