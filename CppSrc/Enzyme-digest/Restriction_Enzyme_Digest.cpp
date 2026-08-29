// Restriction Enzyme Digest Simulator - in silico restriction digestion of DNA.
//
// C++ port of the Python tool at
// https://github.com/Nheyer/restriction-enzyme-digest-simulator, ported from
// commit 82994c5 on main. The CLI surface, the output text and the cut
// coordinate convention are all kept identical to the Python original so the
// two can be differentially tested against each other on the same FASTA files.
//
// Copyright (c) 2026 The Zabel Lab at Colorado State University.
//
// This file is part of Genomics_Scripts_Zlab and is distributed under the
// GPL v3 - see the LICENSE file at the root of the repository. That is the
// licence you receive it under.
//
// It is also a derivative work of MIT licensed code, and that notice is
// reproduced below because keeping it is a condition of the permission that
// allows the code to be used here. MIT is GPL compatible, so this does not put
// the file under two licences: the GPL v3 governs, and the notice is
// attribution for the portion it covers. It stays regardless of how far the
// port diverges - a translation is a derivative work, and the extent of the
// changes does not retire the obligation.
//
//   MIT License
//   Copyright (c) 2026 Collins Amatu Gorgerat
//
//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the "Software"),
//   to deal in the Software without restriction, including without limitation
//   the rights to use, copy, modify, merge, publish, distribute, sublicense,
//   and/or sell copies of the Software, and to permit persons to whom the
//   Software is furnished to do so, subject to the following conditions:
//
//   The above copyright notice and this permission notice shall be included in
//   all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//   DEALINGS IN THE SOFTWARE.
//
// The full notice is also kept in LICENSE.restriction-digest at the repo root.
//
// Where this differs from the Python: matching does not go through a regex
// engine. Each IUPAC code is a 4 bit mask over {A,C,G,T} and a site is matched
// by walking the window and testing masks, so there is no pattern compilation
// and no backtracking. The semantics are unchanged - see match_base().

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <argparse/argparse.hpp>

#include "enzyme_data.hpp"
#include "enzyme_csv_embedded.hpp"

using namespace enzyme_data;

// Single source of truth for the version.
static const char *const VERSION = "0.1.0";

// How to treat IUPAC ambiguity codes in the *target sequence*.
//   definite: the site is cut however the ambiguous bases resolve - the bases
//             the sequence allows are a subset of those the enzyme accepts.
//             Answers "where will this definitely cut?".
//   possible: the site is cut for at least one resolution - the two base sets
//             overlap. Answers "could this cut at all?", the question behind
//             checking that an enzyme leaves a sequence intact.
// Every definite site is also a possible site, never the other way round.
static const char *const AMBIGUITY_DEFINITE = "definite";
static const char *const AMBIGUITY_POSSIBLE = "possible";

// ------------------------------------------------------------ IUPAC codes --

// Bit per base, so a code is a set and the set operations are bitwise.
enum : uint8_t { BASE_A = 1, BASE_C = 2, BASE_G = 4, BASE_T = 8 };

// code -> mask, built once from enzyme_data::IUPAC_BASES. Anything that is not
// an IUPAC code maps to 0, which never matches: that is what keeps a stray
// character in a sequence from silently matching everything.
static const uint8_t *iupac_mask_table() {
    static uint8_t table[256] = {0};
    static bool built = false;
    if (!built) {
        for (int i = 0; i < IUPAC_COUNT; i++) {
            uint8_t mask = 0;
            for (const char *b = IUPAC_BASES[i].bases; *b; b++) {
                if (*b == 'A') { mask |= BASE_A; }
                if (*b == 'C') { mask |= BASE_C; }
                if (*b == 'G') { mask |= BASE_G; }
                if (*b == 'T') { mask |= BASE_T; }
            }
            table[(unsigned char) IUPAC_BASES[i].code] = mask;
        }
        built = true;
    }
    return table;
}

static uint8_t mask_of(char code) { return iupac_mask_table()[(unsigned char) code]; }

// Does a target base satisfy a site code? This is the whole of the ambiguity
// policy, and the direction is asymmetric on purpose.
//
//   definite: every base the target could be is one the site accepts, so the
//             cut happens no matter how the target resolves. Site N matches
//             target A and target N; site A matches target A but NOT target N,
//             because that N might be a C.
//   possible: they share at least one base, so the cut happens for at least
//             one resolution. Site T matches target Y (Y may be a T); site A
//             does not (Y is C or T, neither is A).
//
// A target that is not an IUPAC code has mask 0 and matches under neither.
static bool match_base(uint8_t site_mask, uint8_t target_mask, bool definite) {
    if (target_mask == 0) { return false; }
    if (definite) { return (target_mask & ~site_mask) == 0; }
    return (target_mask & site_mask) != 0;
}

static char complement_code(char code) {
    for (int i = 0; COMPLEMENT_FROM[i]; i++) {
        if (COMPLEMENT_FROM[i] == code) { return COMPLEMENT_TO[i]; }
    }
    return code;
}

static std::string reverse_complement(const std::string &seq) {
    std::string out;
    out.reserve(seq.size());
    for (std::string::const_reverse_iterator it = seq.rbegin(); it != seq.rend(); ++it) {
        out += complement_code(*it);
    }
    return out;
}

// Upper-case a recognition site and reject non-IUPAC characters.
static std::string validate_site(const std::string &raw) {
    std::string site;
    for (char c : raw) {
        if (!isspace((unsigned char) c)) { site += (char) toupper((unsigned char) c); }
    }
    if (site.empty()) { throw std::invalid_argument("Recognition site is empty"); }
    std::set<char> bad;
    for (char c : site) {
        if (mask_of(c) == 0) { bad.insert(c); }
    }
    if (!bad.empty()) {
        std::string listed;
        for (char c : bad) {
            if (!listed.empty()) { listed += ", "; }
            listed += c;
        }
        throw std::invalid_argument("Invalid IUPAC character(s) in recognition site: " + listed);
    }
    return site;
}

// Expand IUPAC ambiguity codes into every concrete ACGT sequence. HincII's
// GTYRAC becomes GTCAAC, GTCGAC, GTTAAC, GTTGAC. Matching does not use this -
// it would explode on XcmI's nine Ns - but --convert reports it.
static std::vector<std::string> expand_ambiguity(const std::string &raw,
                                                 unsigned long max_expansions = 65536) {
    std::string site = validate_site(raw);
    unsigned long total = 1;
    for (char code : site) {
        uint8_t mask = mask_of(code);
        int degeneracy = 0;
        for (int bit = 0; bit < 4; bit++) {
            if (mask & (1u << bit)) { degeneracy++; }
        }
        total *= (unsigned long) degeneracy;
        if (total > max_expansions) {
            std::ostringstream oss;
            oss << "Expanding '" << site << "' would give more than " << max_expansions
                << " sequences; raise max_expansions or match with site_to_regex() instead.";
            throw std::invalid_argument(oss.str());
        }
    }
    std::vector<std::string> results(1, std::string());
    // Ascending bit order is ACGT order, which is the order the Python
    // IUPAC_BASES strings are written in, so the variants come out the same.
    static const char BASES[] = "ACGT";
    for (char code : site) {
        uint8_t mask = mask_of(code);
        std::vector<std::string> next;
        for (const std::string &prefix : results) {
            for (int bit = 0; bit < 4; bit++) {
                if (mask & (1u << bit)) { next.push_back(prefix + BASES[bit]); }
            }
        }
        results.swap(next);
    }
    return results;
}

// ---------------------------------------------------------------- enzymes --

// cut_top and cut_bottom are 0-based offsets from the first base of the
// recognition site, both measured along the top strand: the cut falls
// immediately before the base at that offset. cut_bottom is where the enzyme
// nicks the complementary strand, expressed in the same top-strand frame, so
// cut_bottom - cut_top is the overhang. Offsets may lie outside the site -
// type IIS enzymes such as BsaI cut downstream of it.
struct Enzyme {
    std::string name;
    std::string site;
    int cut_top = 0;
    int cut_bottom = 0;
};

// True if the site reads the same on both strands (EcoRI, not BsaI).
static bool is_palindromic(const Enzyme &e) { return e.site == reverse_complement(e.site); }

// True if binding in either orientation gives the same pair of cuts. Note this
// is strictly stronger than palindromic: it also needs the cuts to mirror, and
// it is what decides whether the reverse orientation is searched at all.
static bool is_symmetric(const Enzyme &e) {
    return is_palindromic(e) && e.cut_bottom == (int) e.site.size() - e.cut_top;
}

// Positive for a 5' overhang, 0 for a blunt cut, negative for a 3' overhang.
static int overhang(const Enzyme &e) { return e.cut_bottom - e.cut_top; }

static std::string describe_ends(const Enzyme &e) {
    int hang = overhang(e);
    std::ostringstream oss;
    if (hang > 0) { oss << "5' overhang, " << hang << " nt"; return oss.str(); }
    if (hang < 0) { oss << "3' overhang, " << -hang << " nt"; return oss.str(); }
    return "blunt";
}

// Render an enzyme as a NEB/REBASE recognition specificity string.
static std::string to_neb_notation(const Enzyme &e) {
    int length = (int) e.site.size();
    if (e.cut_top >= 0 && e.cut_top <= length && e.cut_bottom == length - e.cut_top) {
        return e.site.substr(0, e.cut_top) + "^" + e.site.substr(e.cut_top);
    }
    std::ostringstream oss;
    oss << e.site << "(" << (e.cut_top - length) << "/" << (e.cut_bottom - length) << ")";
    return oss.str();
}

// Render an enzyme as the NAME:SEQ:OFFSET spec the tool started with. That form
// records only the top-strand cut and assumes the bottom-strand cut is its
// mirror, so it cannot represent an asymmetric cutter. Rather than silently
// dropping cut_bottom it refuses.
static std::string to_legacy_spec(const Enzyme &e) {
    int length = (int) e.site.size();
    if (!(e.cut_top >= 0 && e.cut_top <= length)) {
        throw std::invalid_argument(
            e.name + " cuts outside its recognition site (" + to_neb_notation(e) +
            "); NAME:SEQ:OFFSET cannot express that. Use NAME:" + to_neb_notation(e) +
            " instead.");
    }
    if (e.cut_bottom != length - e.cut_top) {
        throw std::invalid_argument(
            e.name + " cuts the two strands asymmetrically (" + to_neb_notation(e) +
            "); NAME:SEQ:OFFSET would lose the bottom-strand cut. Use NAME:" +
            to_neb_notation(e) + " instead.");
    }
    std::ostringstream oss;
    oss << e.name << ":" << e.site << ":" << e.cut_top;
    return oss.str();
}

// Parse a NEB/REBASE recognition specificity into site, cut_top, cut_bottom.
// Accepts the two notations NEB publishes: a caret inside the site (G^AATTC,
// GGTAC^C, ^GATC) or cut offsets after it (GGTCTC(1/5), GAATTC(-5/-1)).
// Enzymes that cut on both sides of their site are rejected: two cuts per site
// do not fit the single-cut fragment model.
static void parse_neb_notation(const std::string &raw, std::string *site_out,
                               int *top_out, int *bottom_out) {
    std::string text;
    for (char c : raw) {
        if (!isspace((unsigned char) c)) { text += (char) toupper((unsigned char) c); }
    }
    if (text.empty()) { throw std::invalid_argument("Empty recognition specificity"); }
    if (text[0] == '(') {
        throw std::invalid_argument(
            "'" + text + "' cuts on both sides of its recognition site; dual-cut enzymes "
            "(BaeI, BcgI, CspCI ...) are not supported.");
    }

    static const std::regex OFFSETS("([^()]+)\\((-?[0-9]+)/(-?[0-9]+)\\)");
    std::smatch m;
    if (std::regex_match(text, m, OFFSETS)) {
        if (m[1].str().find('^') != std::string::npos) {
            throw std::invalid_argument(
                "'" + text + "' mixes caret and offset notation; use one or the other");
        }
        std::string site = validate_site(m[1].str());
        int length = (int) site.size();
        *site_out = site;
        *top_out = length + std::stoi(m[2].str());
        *bottom_out = length + std::stoi(m[3].str());
        return;
    }

    if (text.find('^') != std::string::npos) {
        if (std::count(text.begin(), text.end(), '^') > 1) {
            throw std::invalid_argument(
                "'" + text + "' has more than one cut position marker '^'");
        }
        int cut_top = (int) text.find('^');
        std::string bare = text;
        bare.erase(std::remove(bare.begin(), bare.end(), '^'), bare.end());
        std::string site = validate_site(bare);
        // A caret marks only the top-strand cut; the bottom cut is its mirror.
        *site_out = site;
        *top_out = cut_top;
        *bottom_out = (int) site.size() - cut_top;
        return;
    }

    validate_site(text);  // surface bad characters before the generic message
    throw std::invalid_argument(
        "No cut position in '" + text + "'. Mark it with '^' (G^AATTC), append offsets "
        "(GGTCTC(1/5)), or use NAME:SEQ:OFFSET.");
}

// Split a ^/_ marked specificity into site, top cuts, bottom cuts.
//
// This is the notation NEB publishes for its commercially available enzymes:
// '^' marks the top-strand cut, '_' the bottom-strand cut, and Ns pad the gap
// for enzymes that cut away from what they recognise:
//
//   G^AATTC          EcoRI, both cuts inside the site
//   G_ACGT^C         AatII, a 3' overhang
//   GGTCTCN^NNNN_    BsaI, cutting downstream
//
// Offsets are in site coordinates, so the markers take up no position of their
// own. Leading and trailing Ns are spacing rather than recognition, so they are
// trimmed and the offsets shifted to match - BsaI comes back as
// ("GGTCTC", 7, 11), which renders as GGTCTC(1/5).
static void split_cut_markers(const std::string &text, std::string *site_out,
                              std::vector<int> *top_out, std::vector<int> *bottom_out) {
    std::string site;
    std::vector<int> top, bottom;
    for (char raw : text) {
        if (isspace((unsigned char) raw)) { continue; }
        char c = (char) toupper((unsigned char) raw);
        if (c == '^') {
            top.push_back((int) site.size());
        } else if (c == '_') {
            bottom.push_back((int) site.size());
        } else {
            site += c;
        }
    }
    size_t first = site.find_first_not_of('N');
    if (first != std::string::npos) {
        size_t last = site.find_last_not_of('N');
        int lead = (int) first;
        for (int &value : top) { value -= lead; }
        for (int &value : bottom) { value -= lead; }
        site = site.substr(first, last - first + 1);
    }
    *site_out = site;
    *top_out = top;
    *bottom_out = bottom;
}

// Enzymes from a source file, split by what this tool can simulate.
struct EnzymeTables {
    // Insertion ordered, like the Python dict this mirrors: --list-enzymes sorts
    // case-insensitively and ties fall back on this order.
    std::vector<Enzyme> enzymes;
    std::vector<std::pair<std::string, std::string> > nicking;
    std::vector<std::pair<std::string, std::string> > dual_cut;
};

static std::string upper_of(const std::string &s) {
    std::string out = s;
    for (char &c : out) { c = (char) toupper((unsigned char) c); }
    return out;
}

// Split one CSV line. The enzyme CSV has no quoting or embedded newlines, so a
// plain comma split matches what Python's csv.DictReader does with it.
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

static std::string trimmed(const std::string &s) {
    size_t first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) { return ""; }
    size_t last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

// Read an "enzyme,recognition_sequence" CSV of NEB specificities.
//
// Nicking enzymes (one strand only) and type IIB enzymes (a cut on each side of
// the site) are kept aside rather than loaded: neither fits the one-cut fragment
// model, and naming them gives a better error than "unknown enzyme".
static EnzymeTables load_enzyme_csv_string(const std::string &text, const std::string &path) {
    EnzymeTables tables;
    std::istringstream stream(text);
    std::string line;
    int lineno = 1;
    if (!std::getline(stream, line)) {
        throw std::invalid_argument("No usable enzymes found in " + path);
    }
    std::vector<std::string> header = split_csv_line(line);
    int name_col = -1, spec_col = -1;
    for (size_t i = 0; i < header.size(); i++) {
        std::string field = trimmed(header[i]);
        if (field == "enzyme") { name_col = (int) i; }
        if (field == "recognition_sequence") { spec_col = (int) i; }
    }
    if (name_col < 0 || spec_col < 0) {
        // Names listed sorted, to match the Python message.
        std::string missing;
        if (name_col < 0) { missing = "enzyme"; }
        if (spec_col < 0) {
            missing += (missing.empty() ? "" : ", ");
            missing += "recognition_sequence";
        }
        throw std::invalid_argument(
            path + " is missing column(s) " + missing +
            "; expected a header 'enzyme,recognition_sequence'.");
    }

    while (std::getline(stream, line)) {
        lineno++;
        std::vector<std::string> row = split_csv_line(line);
        std::string name = (name_col < (int) row.size()) ? trimmed(row[name_col]) : "";
        std::string spec = (spec_col < (int) row.size()) ? trimmed(row[spec_col]) : "";
        if (name.empty() || spec.empty()) { continue; }

        std::string site;
        std::vector<int> top, bottom;
        split_cut_markers(spec, &site, &top, &bottom);

        std::ostringstream where;
        where << path << ", line " << lineno << ": ";

        if (top.size() > 1 || bottom.size() > 1) {
            tables.dual_cut.push_back(std::make_pair(name, spec));   // type IIB
        } else if (top.empty() && bottom.empty()) {
            throw std::invalid_argument(
                where.str() + "no cut marker in '" + spec + "'; mark the top-strand cut "
                "with '^' and the bottom-strand cut with '_'.");
        } else if ((top.empty() || bottom.empty()) &&
                   (name.rfind("Nb.", 0) == 0 || name.rfind("Nt.", 0) == 0)) {
            // REBASE names nicking enzymes Nb.* (bottom) and Nt.* (top), and
            // writes them with the one marker for the strand they cut.
            tables.nicking.push_back(std::make_pair(name, spec));
        } else {
            // One marker is the usual NEB shorthand for a cut whose partner
            // mirrors it: G^AATTC means GAATTC with cuts at 1 and 5.
            if (bottom.empty()) { bottom.push_back((int) site.size() - top[0]); }
            else if (top.empty()) { top.push_back((int) site.size() - bottom[0]); }
            try {
                site = validate_site(site);
            } catch (const std::invalid_argument &exc) {
                throw std::invalid_argument(where.str() + exc.what());
            }
            Enzyme e;
            e.name = name;
            e.site = site;
            e.cut_top = top[0];
            e.cut_bottom = bottom[0];
            tables.enzymes.push_back(e);
        }
    }
    if (tables.enzymes.empty()) {
        throw std::invalid_argument("No usable enzymes found in " + path);
    }
    return tables;
}

static EnzymeTables load_enzyme_csv(const std::string &path) {
    std::ifstream handle(path.c_str());
    if (!handle) {
        throw std::runtime_error("[Errno 2] No such file or directory: '" + path + "'");
    }
    std::ostringstream buffer;
    buffer << handle.rdbuf();
    return load_enzyme_csv_string(buffer.str(), path);
}

// The built-in table, parsed once from the CSV embedded at build time. Keeping
// the CSV as the only copy of the data is deliberate: see ROADMAP.md upstream.
static const EnzymeTables &builtin_tables() {
    static EnzymeTables tables =
        load_enzyme_csv_string(EMBEDDED_ENZYME_CSV, "restriction_enzymes.csv");
    return tables;
}

// Parse an enzyme specification into an Enzyme. Accepted forms:
//   EcoRI                - a built-in, matched case-insensitively
//   NAME:G^AATTC         - NEB notation, caret form
//   NAME:GGTCTC(1/5)     - NEB notation, offset form (non-palindromic)
//   NAME:SEQ:OFFSET      - top-strand cut, bottom strand mirrored
//   NAME:SEQ:TOP:BOTTOM  - both cuts, as offsets from the site start
static Enzyme parse_enzyme_spec(const std::string &raw, const EnzymeTables *tables = NULL) {
    std::string spec = trimmed(raw);
    const EnzymeTables &db = tables ? *tables : builtin_tables();
    std::string wanted = upper_of(spec);

    for (const Enzyme &e : db.enzymes) {
        if (upper_of(e.name) == wanted) { return e; }
    }

    // A name we know but cannot simulate deserves a better error than "unknown".
    const char *NICK_REASON = "nicks one strand and leaves the duplex intact, so it "
                              "produces no fragments to simulate";
    const char *DUAL_REASON = "cuts on both sides of its recognition site; two cuts "
                              "per site do not fit the single-cut fragment model";
    for (int table = 0; table < 2; table++) {
        const std::vector<std::pair<std::string, std::string> > &rows =
            table == 0 ? db.nicking : db.dual_cut;
        const char *reason = table == 0 ? NICK_REASON : DUAL_REASON;
        for (const std::pair<std::string, std::string> &row : rows) {
            if (upper_of(row.first) == wanted) {
                throw std::invalid_argument(
                    row.first + " (" + row.second + ") " + reason + ".");
            }
        }
    }

    std::vector<std::string> parts;
    {
        std::string current;
        for (char c : spec) {
            if (c == ':') { parts.push_back(current); current.clear(); }
            else { current += c; }
        }
        parts.push_back(current);
    }

    if (parts.size() == 1) {
        if (spec.find('^') != std::string::npos || spec.find('(') != std::string::npos) {
            // A bare recognition specificity, e.g. G^AATTC or GGTCTC(1/5).
            Enzyme e;
            e.name = "custom";
            parse_neb_notation(spec, &e.site, &e.cut_top, &e.cut_bottom);
            return e;
        }
        throw std::invalid_argument(
            "Unknown enzyme '" + spec + "'. Provide custom as NAME:SEQ:OFFSET or "
            "NAME:RECOGNITION_SITE (e.g. MyI:G^AATTC, MyI:GGTCTC(1/5)).");
    }

    std::string name = trimmed(parts[0]);
    if (name.empty()) {
        throw std::invalid_argument("Malformed enzyme spec: '" + spec + "' (missing enzyme name)");
    }

    if (parts.size() == 2) {
        std::string recog = trimmed(parts[1]);
        if (recog.find('^') == std::string::npos && recog.find('(') == std::string::npos) {
            throw std::invalid_argument(
                "Cut offset required for custom enzyme '" + spec + "'. Use NAME:SEQ:OFFSET or "
                "mark the cut in the site itself (NAME:" + recog.substr(0, 1) + "^" +
                (recog.size() > 1 ? recog.substr(1) : std::string()) + ").");
        }
        Enzyme e;
        e.name = name;
        parse_neb_notation(recog, &e.site, &e.cut_top, &e.cut_bottom);
        return e;
    }

    if (parts.size() == 3 || parts.size() == 4) {
        std::string site = validate_site(parts[1]);
        std::vector<int> cuts;
        for (size_t i = 2; i < parts.size(); i++) {
            std::string value = trimmed(parts[i]);
            size_t consumed = 0;
            int parsed = 0;
            bool ok = !value.empty();
            if (ok) {
                try { parsed = std::stoi(value, &consumed); }
                catch (const std::exception &) { ok = false; }
            }
            if (!ok || consumed != value.size()) {
                throw std::invalid_argument("Offset must be integer, got '" + value + "'");
            }
            cuts.push_back(parsed);
        }
        Enzyme e;
        e.name = name;
        e.site = site;
        e.cut_top = cuts[0];
        e.cut_bottom = cuts.size() == 2 ? cuts[1] : (int) site.size() - cuts[0];
        if (cuts.size() == 1 && !(e.cut_top >= 0 && e.cut_top <= (int) site.size())) {
            std::ostringstream oss;
            oss << "Offset " << e.cut_top << " out of range for recognition seq length "
                << site.size();
            throw std::invalid_argument(oss.str());
        }
        return e;
    }

    throw std::invalid_argument("Malformed enzyme spec: '" + spec + "'");
}

// ----------------------------------------------------------- site finding --

// One recognition site found in a sequence. start is where the site begins on
// the top strand, strand is '+' when the site reads in the given sequence's
// direction and '-' when it was found as the reverse complement, and cut is the
// top-strand cut position in sequence coordinates, before normalisation.
struct SiteHit {
    int start = 0;
    char strand = '+';
    int cut = 0;
};

static bool site_hit_less(const SiteHit &a, const SiteHit &b) {
    if (a.start != b.start) { return a.start < b.start; }
    if (a.strand != b.strand) { return a.strand < b.strand; }
    return a.cut < b.cut;
}

// Precomputed masks for one orientation of a site.
static std::vector<uint8_t> site_masks(const std::string &site) {
    std::vector<uint8_t> masks;
    masks.reserve(site.size());
    for (char c : site) { masks.push_back(mask_of(c)); }
    return masks;
}

// Find every recognition site for an enzyme, on both strands.
//
// A non-palindromic enzyme such as BsaI recognises its site on either strand,
// and a site on the bottom strand cuts the top strand on the other side of the
// site, so both orientations are searched. For a fully symmetric enzyme the two
// searches coincide, so each site is reported once.
//
// ambiguity selects how IUPAC codes in the *sequence* are judged. When circular
// is set, sites that straddle the origin are found too.
static std::vector<SiteHit> find_sites(const std::string &sequence, const Enzyme &enzyme,
                                       const std::string &ambiguity = AMBIGUITY_DEFINITE,
                                       bool circular = false) {
    const bool definite = (ambiguity != AMBIGUITY_POSSIBLE);
    const int length = (int) enzyme.site.size();
    const int seq_len = (int) sequence.size();

    // A circular molecule has no ends: extend the scan window by length-1 bases
    // so a site reading through the join is contiguous. Start indices stay
    // within 0..seq_len-1, so no site is counted twice.
    std::string scan = sequence;
    if (circular && length > 1 && length <= seq_len) {
        scan += sequence.substr(0, length - 1);
    }
    const int scan_len = (int) scan.size();

    // Target masks for the whole window, computed once for both orientations.
    std::vector<uint8_t> target(scan_len);
    for (int i = 0; i < scan_len; i++) { target[i] = mask_of(scan[i]); }

    struct Orientation { char strand; std::string pattern; };
    std::vector<Orientation> orientations;
    orientations.push_back(Orientation{'+', enzyme.site});
    if (!is_symmetric(enzyme)) {
        orientations.push_back(Orientation{'-', reverse_complement(enzyme.site)});
    }

    // Deduplicated on (start, cut), '+' inserted first so a collision keeps it.
    std::vector<SiteHit> hits;
    std::set<std::pair<int, int> > seen;
    for (const Orientation &orientation : orientations) {
        std::vector<uint8_t> masks = site_masks(orientation.pattern);
        for (int start = 0; start + length <= scan_len; start++) {
            bool matched = true;
            for (int i = 0; i < length; i++) {
                if (!match_base(masks[i], target[start + i], definite)) { matched = false; break; }
            }
            if (!matched) { continue; }
            SiteHit hit;
            hit.start = start;
            hit.strand = orientation.strand;
            if (orientation.strand == '+') {
                hit.cut = start + enzyme.cut_top;
            } else {
                // Mirror the enzyme frame: an offset p sits at start + length - p,
                // so the enzyme's bottom-strand cut lands on our top strand.
                hit.cut = start + length - enzyme.cut_bottom;
            }
            if (seen.insert(std::make_pair(hit.start, hit.cut)).second) { hits.push_back(hit); }
        }
    }
    std::sort(hits.begin(), hits.end(), site_hit_less);
    return hits;
}

// Bring cut positions into sequence coordinates.
//
// Type IIS enzymes cut a fixed distance away from their site, which can fall
// past the end - or before the start - of the sequence. On circular DNA that
// wraps around; on linear DNA there is nothing there to cut, so the cut is
// dropped. Returns the usable cut positions and how many were dropped.
static std::vector<int> normalise_cuts(const std::vector<int> &cuts, int seq_len,
                                       bool circular, int *dropped_out) {
    std::set<int> kept;
    int dropped = 0;
    for (int cut : cuts) {
        if (circular) {
            // Python's % is a floor mod and always non-negative here; C++'s is
            // not, and a bottom-strand type IIS hit near the origin goes negative.
            if (seq_len) { kept.insert(((cut % seq_len) + seq_len) % seq_len); }
        } else if (cut >= 1 && cut <= seq_len - 1) {
            kept.insert(cut);
        } else {
            dropped++;
        }
    }
    if (dropped_out) { *dropped_out = dropped; }
    return std::vector<int>(kept.begin(), kept.end());
}

// ------------------------------------------------------------- digestion ---

// Compute fragment lengths for linear DNA.
static std::vector<int> digest_linear(int seq_len, const std::vector<int> &cut_positions,
                                      int min_fragment) {
    std::vector<int> all_cuts = cut_positions;
    std::sort(all_cuts.begin(), all_cuts.end());
    std::vector<int> fragments;
    int prev = 0;
    for (int cut : all_cuts) {
        int frag_len = cut - prev;
        if (frag_len >= min_fragment) { fragments.push_back(frag_len); }
        prev = cut;
    }
    int final_frag = seq_len - prev;
    if (final_frag >= min_fragment) { fragments.push_back(final_frag); }
    return fragments;
}

// Compute fragment lengths for circular DNA.
static std::vector<int> digest_circular(int seq_len, const std::vector<int> &cut_positions,
                                        int min_fragment) {
    std::vector<int> fragments;
    if (cut_positions.empty()) {
        if (seq_len >= min_fragment) { fragments.push_back(seq_len); }
        return fragments;
    }
    std::vector<int> sorted_cuts = cut_positions;
    std::sort(sorted_cuts.begin(), sorted_cuts.end());
    int n = (int) sorted_cuts.size();
    for (int i = 0; i < n; i++) {
        int start = sorted_cuts[i];
        int end = sorted_cuts[(i + 1) % n];
        int dist = (i == n - 1) ? (seq_len - start + end) : (end - start);
        if (dist >= min_fragment) { fragments.push_back(dist); }
    }
    return fragments;
}

// ------------------------------------------------------------ FASTA input --

// The set of characters accepted in an input sequence.
static std::set<char> sequence_alphabet(bool allow_ambiguity) {
    std::set<char> allowed;
    if (allow_ambiguity) {
        for (int i = 0; IUPAC_ALPHABET[i]; i++) { allowed.insert(IUPAC_ALPHABET[i]); }
        allowed.insert(GAP_CHAR);
    } else {
        for (int i = 0; CONCRETE_BASES[i]; i++) { allowed.insert(CONCRETE_BASES[i]); }
    }
    return allowed;
}

static std::string describe_alphabet(bool allow_ambiguity) {
    if (allow_ambiguity) {
        return std::string(IUPAC_DISPLAY_ORDER) + " and '" + std::string(1, GAP_CHAR) + "' (gap)";
    }
    return std::string(CONCRETE_BASES);  // already in sorted ACGT order
}

// Python renders these with repr(), which for a single character is the
// character in single quotes unless it is itself a quote.
static std::string py_repr(char c) {
    if (c == '\'') { return "\"'\""; }
    if (c == '\\') { return "'\\\\'"; }
    return std::string("'") + c + "'";
}

static std::string join_reprs(const std::vector<char> &chars) {
    std::string out;
    for (size_t i = 0; i < chars.size(); i++) {
        if (i) { out += ", "; }
        out += py_repr(chars[i]);
    }
    return out;
}

// Spell out every offending character and how to proceed.
static std::string bad_alphabet_message(const std::vector<char> &bad, const std::string &header,
                                        int lineno, bool allow_ambiguity) {
    std::string listed = join_reprs(bad);
    std::ostringstream where;
    where << "sequence '" << header << "', line " << lineno;
    if (allow_ambiguity) {
        return "REJECTED - invalid character(s) " + listed + " in " + where.str() +
               ". Allowed: " + describe_alphabet(true) + ".";
    }

    std::set<char> permissive = sequence_alphabet(true);
    std::vector<char> codes, unknown;
    for (char c : bad) {
        if (permissive.count(c)) { codes.push_back(c); } else { unknown.push_back(c); }
    }
    std::string msg = "non-ACGT character(s) " + listed + " in " + where.str() + ".";
    if (!unknown.empty()) {
        msg += " " + join_reprs(unknown) + " " + (unknown.size() > 1 ? "are" : "is") +
               " not valid DNA under any setting.";
    }
    if (!codes.empty()) {
        msg += " " + join_reprs(codes) + " " +
               (codes.size() > 1 ? "are IUPAC ambiguity/gap codes" : "is an IUPAC ambiguity/gap code") +
               ", which this tool refuses to guess at: rerun with --ambiguity definite to report"
               " only cuts that are certain, or --ambiguity possible to report every cut that"
               " could happen. Without --ambiguity only A/C/G/T is accepted, so an ambiguous"
               " base is never silently treated as a non-match.";
    }
    return msg;
}

// Drop alignment gaps, returning the sequence and how many were removed.
//
// A gap is not a base - the molecule reads straight through it - so 'AAT-ATT'
// really is an SspI site. Removing them keeps cut positions and fragment
// lengths in real (ungapped) coordinates.
static std::string strip_gaps(const std::string &sequence, int *removed_out) {
    std::string stripped;
    stripped.reserve(sequence.size());
    for (char c : sequence) {
        if (c != GAP_CHAR) { stripped += c; }
    }
    if (removed_out) { *removed_out = (int) (sequence.size() - stripped.size()); }
    return stripped;
}

// Parse a FASTA file into (header, sequence) pairs. Accepts A/C/G/T only,
// unless allow_ambiguity also permits every IUPAC code and the '-' gap.
static std::vector<std::pair<std::string, std::string> > parse_fasta(const std::string &filepath,
                                                                     bool allow_ambiguity) {
    std::set<char> allowed = sequence_alphabet(allow_ambiguity);
    std::vector<std::pair<std::string, std::string> > sequences;
    std::ifstream handle(filepath.c_str());
    if (!handle) {
        std::cerr << "Error: File not found: " << filepath << std::endl;
        std::exit(1);
    }
    std::string line;
    bool have_header = false;
    std::string current_header;
    std::string current_seq;
    int lineno = 0;
    while (std::getline(handle, line)) {
        lineno++;
        std::string text = trimmed(line);
        if (text.empty()) { continue; }
        if (text[0] == '>') {
            if (have_header) { sequences.push_back(std::make_pair(current_header, current_seq)); }
            current_header = trimmed(text.substr(1));
            current_seq.clear();
            have_header = true;
        } else {
            if (!have_header) {
                std::cerr << "Error parsing FASTA: FASTA file missing header line" << std::endl;
                std::exit(1);
            }
            std::string upper = upper_of(text);
            std::set<char> bad_set;
            for (char c : upper) {
                if (!allowed.count(c)) { bad_set.insert(c); }
            }
            if (!bad_set.empty()) {
                std::vector<char> bad(bad_set.begin(), bad_set.end());
                std::cerr << "Error parsing FASTA: "
                          << bad_alphabet_message(bad, current_header, lineno, allow_ambiguity)
                          << std::endl;
                std::exit(1);
            }
            current_seq += upper;
        }
    }
    if (have_header) {
        sequences.push_back(std::make_pair(current_header, current_seq));
    } else if (sequences.empty()) {
        std::cerr << "Error parsing FASTA: No sequences found in FASTA file" << std::endl;
        std::exit(1);
    }
    return sequences;
}

// ---------------------------------------------------------------- output ---

// Python's str.format alignment, which the expected output is written against.
static std::string pad_left(const std::string &s, size_t width) {
    return s.size() >= width ? s : std::string(width - s.size(), ' ') + s;
}
static std::string pad_right(const std::string &s, size_t width) {
    return s.size() >= width ? s : s + std::string(width - s.size(), ' ');
}
static std::string pad_center(const std::string &s, size_t width) {
    if (s.size() >= width) { return s; }
    size_t total = width - s.size();
    size_t left = total / 2;                 // Python puts the odd space on the right
    return std::string(left, ' ') + s + std::string(total - left, ' ');
}

static std::string to_string_int(int value) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

// Python's repr of a list of ints: [1, 2, 3].
static std::string render_int_list(const std::vector<int> &values) {
    std::string out = "[";
    for (size_t i = 0; i < values.size(); i++) {
        if (i) { out += ", "; }
        out += to_string_int(values[i]);
    }
    return out + "]";
}

static void print_fragment_table(const std::vector<int> &fragments, const std::string &note) {
    std::cout << "\n--- Fragment Summary ---" << std::endl;
    if (!note.empty()) { std::cout << note << std::endl; }
    std::vector<int> sorted_fragments = fragments;
    std::sort(sorted_fragments.begin(), sorted_fragments.end(), std::greater<int>());
    long long total_len = 0;
    for (int size : sorted_fragments) { total_len += size; }
    std::cout << pad_left("#", 3) << " " << pad_left("Size (bp)", 9) << std::endl;
    std::cout << std::string(14, '-') << std::endl;
    for (size_t i = 0; i < sorted_fragments.size(); i++) {
        std::cout << pad_left(to_string_int((int) i + 1), 3) << " "
                  << pad_left(to_string_int(sorted_fragments[i]), 9) << std::endl;
    }
    std::cout << std::string(14, '-') << std::endl;
    std::ostringstream total;
    total << total_len;
    std::cout << pad_left("Total", 3) << " " << pad_left(total.str(), 9) << std::endl;
    std::cout << std::endl;
}

// Draw an ASCII gel with a ladder lane and a sample lane.
static void draw_ascii_gel(const std::vector<int> &fragments, const int *ladder_sizes,
                           int ladder_count, int height) {
    std::vector<int> all_sizes(ladder_sizes, ladder_sizes + ladder_count);
    all_sizes.insert(all_sizes.end(), fragments.begin(), fragments.end());
    if (all_sizes.empty()) {
        std::cout << "No fragments to display." << std::endl;
        return;
    }
    int max_s = *std::max_element(all_sizes.begin(), all_sizes.end());
    int min_s = *std::min_element(all_sizes.begin(), all_sizes.end());
    double max_log = std::log10((double) std::max(max_s, 1));
    double min_log = std::log10((double) std::max(min_s, 1));
    if (max_log == min_log) { max_log += 0.1; }  // avoid division by zero

    // Python's round() is half-to-even; std::nearbyint matches it under the
    // default rounding mode, where std::round (half-away-from-zero) would not.
    struct RowFor {
        double max_log, min_log;
        int height;
        int operator()(int size) const {
            double log_s = size > 0 ? std::log10((double) size) : min_log;
            double frac = (max_log - log_s) / (max_log - min_log);  // larger -> top (row 0)
            int row = (int) std::nearbyint(frac * (height - 1));
            return std::max(0, std::min(height - 1, row));
        }
    };
    RowFor row_for{max_log, min_log, height};

    std::map<int, std::vector<int> > ladder_rows;
    for (int i = 0; i < ladder_count; i++) { ladder_rows[row_for(ladder_sizes[i])].push_back(ladder_sizes[i]); }

    std::set<int> sample_bands;
    for (int size : fragments) { sample_bands.insert(row_for(size)); }

    std::cout << pad_left("Ladder (bp)", 12) << " | " << pad_center("Sample", 5) << std::endl;
    std::cout << std::string(25, '-') << std::endl;
    for (int r = 0; r < height; r++) {
        std::string labels;
        std::map<int, std::vector<int> >::iterator found = ladder_rows.find(r);
        if (found != ladder_rows.end()) {
            std::vector<int> sizes = found->second;
            std::sort(sizes.begin(), sizes.end(), std::greater<int>());
            for (size_t i = 0; i < sizes.size(); i++) {
                if (i) { labels += ","; }
                labels += to_string_int(sizes[i]);
            }
        }
        std::string band_char = sample_bands.count(r) ? "X" : "|";
        std::cout << pad_left(labels, 12) << " | " << pad_center(band_char, 5) << std::endl;
    }
}

// Case-insensitive name order, ties falling back on table order like Python's
// stable sorted(..., key=str.lower) over an insertion ordered dict.
static std::vector<int> name_sorted_indices(const std::vector<std::string> &names) {
    std::vector<int> order;
    for (size_t i = 0; i < names.size(); i++) { order.push_back((int) i); }
    std::stable_sort(order.begin(), order.end(), [&names](int a, int b) {
        std::string la = names[a], lb = names[b];
        for (char &c : la) { c = (char) tolower((unsigned char) c); }
        for (char &c : lb) { c = (char) tolower((unsigned char) c); }
        return la < lb;
    });
    return order;
}

static void print_enzyme_list(const EnzymeTables *tables) {
    const EnzymeTables &db = tables ? *tables : builtin_tables();
    std::cout << pad_right("Enzyme", 11) << " " << pad_right("Recognition site", 24) << " "
              << pad_right("Ends", 18) << " Type" << std::endl;
    std::cout << std::string(68, '-') << std::endl;

    std::vector<std::string> names;
    for (const Enzyme &e : db.enzymes) { names.push_back(e.name); }
    for (int index : name_sorted_indices(names)) {
        const Enzyme &e = db.enzymes[index];
        std::string kind = is_palindromic(e) ? "palindromic" : "non-palindromic";
        if (!(e.cut_top >= 0 && e.cut_top <= (int) e.site.size())) { kind += ", cuts outside site"; }
        std::cout << pad_right(e.name, 11) << " " << pad_right(to_neb_notation(e), 24) << " "
                  << pad_right(describe_ends(e), 18) << " " << kind << std::endl;
    }
    std::cout << "\n" << db.enzymes.size()
              << " enzymes. Specificities from REBASE (rebase.neb.com)." << std::endl;

    for (int table = 0; table < 2; table++) {
        const std::vector<std::pair<std::string, std::string> > &rows =
            table == 0 ? db.nicking : db.dual_cut;
        if (rows.empty()) { continue; }
        std::vector<std::string> row_names;
        for (const std::pair<std::string, std::string> &row : rows) { row_names.push_back(row.first); }
        std::string listed;
        for (int index : name_sorted_indices(row_names)) {
            if (!listed.empty()) { listed += ", "; }
            listed += row_names[index];
        }
        if (table == 0) {
            std::cout << rows.size() << " nicking enzyme(s) not listed - they cut one strand only "
                      << "and leave no fragments: " << listed << std::endl;
        } else {
            std::cout << rows.size() << " type IIB enzyme(s) not listed - they cut on both sides of "
                      << "their site: " << listed << std::endl;
        }
    }
}

// Show an enzyme spec in both notations.
static void print_conversion(const std::string &spec) {
    Enzyme enzyme = parse_enzyme_spec(spec);
    std::cout << "Enzyme:            " << enzyme.name << std::endl;
    std::cout << "Recognition site:  " << to_neb_notation(enzyme) << std::endl;
    // Build the value before printing the label: to_legacy_spec refuses on
    // enzymes the legacy form would garble, and a half-written line would show
    // the label twice.
    std::string legacy;
    try {
        legacy = to_legacy_spec(enzyme);
    } catch (const std::invalid_argument &e) {
        legacy = std::string("not representable - ") + e.what();
    }
    std::cout << "NAME:SEQ:OFFSET:   " << legacy << std::endl;
    std::cout << "Cut (top/bottom):  " << enzyme.cut_top << "/" << enzyme.cut_bottom
              << " from the start of the site" << std::endl;
    std::cout << "Ends:              " << describe_ends(enzyme) << std::endl;
    std::cout << "Palindromic:       " << (is_palindromic(enzyme) ? "yes" : "no") << std::endl;

    bool has_codes = false;
    for (char c : enzyme.site) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { has_codes = true; break; }
    }
    if (has_codes) {
        try {
            std::vector<std::string> variants = expand_ambiguity(enzyme.site);
            std::string shown;
            for (size_t i = 0; i < variants.size() && i < 8; i++) {
                if (i) { shown += ", "; }
                shown += variants[i];
            }
            if (variants.size() > 8) { shown += ", ..."; }
            std::cout << "Ambiguity codes:   " << variants.size() << " concrete sites ("
                      << shown << ")" << std::endl;
        } catch (const std::invalid_argument &e) {
            std::cout << "Ambiguity codes:   too degenerate to expand - " << e.what() << std::endl;
        }
    }
}

// ------------------------------------------------------------------ main ---

int main(int argc, char *argv[]) {
    // Named for the binary we install, not for the Python original's
    // enzyme_digest.py. This is the one place the CLI surface deliberately
    // differs from upstream: it shows up in the usage line and in the
    // required-argument error, so those two do not compare byte for byte.
    argparse::ArgumentParser program("Enzyme-digest", VERSION);
    program.add_description(
        "Restriction Enzyme Digest Simulator - in silico restriction digestion of DNA.");

    program.add_argument("--fasta", "-f").help("Path to FASTA file");
    program.add_argument("--enzymes", "-e")
        .help("Comma-separated enzyme names or custom specs "
              "(NAME:SEQ:OFFSET, NAME:G^AATTC, NAME:GGTCTC(1/5))");
    program.add_argument("--circular", "-c")
        .help("Treat DNA as circular (default: linear)")
        .flag();
    program.add_argument("--enzyme-db")
        .metavar("CSV")
        .help("Read enzymes from an \"enzyme,recognition_sequence\" CSV of NEB specificities "
              "instead of the built-in table. The built-in table is generated from "
              "restriction_enzymes.csv in the same format");
    program.add_argument("--ambiguity", "-a")
        .help("Permit IUPAC ambiguity codes and \"-\" gaps in the input sequence, and say how "
              "to treat them: \"definite\" counts only sites cut however the ambiguous bases "
              "resolve; \"possible\" also counts sites that could be cut, for checking that an "
              "enzyme cannot cut at all. Without this flag only A/C/G/T is accepted")
        .choices(AMBIGUITY_DEFINITE, AMBIGUITY_POSSIBLE);
    program.add_argument("--min-fragment", "-m")
        .help("Minimum fragment length to report (default: 1)")
        .default_value(1)
        .scan<'i', int>();
    program.add_argument("--output", "-o")
        .help("Output format: table (default), gel, or both")
        .default_value(std::string("table"))
        .choices("table", "gel", "both");
    program.add_argument("--gel-height")
        .help("Height of ASCII gel (default: 30)")
        .default_value(30)
        .scan<'i', int>();
    program.add_argument("--list-enzymes")
        .help("List the built-in enzymes and exit")
        .flag();
    program.add_argument("--convert")
        .metavar("SPEC")
        .help("Show an enzyme in both notations and exit");

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 2;
    }

    // --enzyme-db swaps the whole table out, so it has to be read before
    // anything that resolves an enzyme name.
    EnzymeTables loaded;
    const EnzymeTables *tables = NULL;
    if (program.is_used("--enzyme-db")) {
        std::string path = program.get<std::string>("--enzyme-db");
        try {
            loaded = load_enzyme_csv(path);
        } catch (const std::exception &exc) {
            std::cerr << "Error reading enzyme database: " << exc.what() << std::endl;
            return 1;
        }
        tables = &loaded;
        std::cerr << "Loaded " << loaded.enzymes.size() << " enzymes from " << path
                  << " (" << loaded.nicking.size() << " nicking, " << loaded.dual_cut.size()
                  << " dual-cut not simulated)" << std::endl;
    }

    if (program.get<bool>("--list-enzymes")) {
        print_enzyme_list(tables);
        return 0;
    }

    if (program.is_used("--convert")) {
        try {
            print_conversion(program.get<std::string>("--convert"));
        } catch (const std::exception &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
        return 0;
    }

    // Both are "required" only at this point, so --list-enzymes and --convert
    // still work without them. An empty string counts as missing, the way it
    // does in Python, where the check is a plain truthiness test.
    const std::string fasta_path =
        program.is_used("--fasta") ? program.get<std::string>("--fasta") : std::string();
    const std::string enzymes_arg =
        program.is_used("--enzymes") ? program.get<std::string>("--enzymes") : std::string();
    if (fasta_path.empty() || enzymes_arg.empty()) {
        std::cerr << program;
        std::cerr << "Enzyme-digest: error: the following arguments are required: "
                  << "--fasta/-f, --enzymes/-e" << std::endl;
        return 2;
    }

    // Passing --ambiguity at all is what opens the input alphabet beyond ACGT;
    // its value then selects how ambiguous bases are judged.
    const bool allow_ambiguity = program.is_used("--ambiguity");
    const std::string ambiguity =
        allow_ambiguity ? program.get<std::string>("--ambiguity") : AMBIGUITY_DEFINITE;
    const bool circular = program.get<bool>("--circular");
    const int min_fragment = program.get<int>("--min-fragment");
    const std::string output = program.get<std::string>("--output");
    const int gel_height = program.get<int>("--gel-height");

    std::vector<std::pair<std::string, std::string> > sequences =
        parse_fasta(fasta_path, allow_ambiguity);

    std::vector<std::string> enzyme_specs;
    {
        const std::string &raw = enzymes_arg;
        std::string current;
        for (char c : raw) {
            if (c == ',') { enzyme_specs.push_back(current); current.clear(); }
            else { current += c; }
        }
        enzyme_specs.push_back(current);
        std::vector<std::string> kept;
        for (const std::string &spec : enzyme_specs) {
            std::string t = trimmed(spec);
            if (!t.empty()) { kept.push_back(t); }
        }
        enzyme_specs.swap(kept);
    }
    if (enzyme_specs.empty()) {
        std::cerr << "Error: No enzymes specified." << std::endl;
        return 1;
    }

    std::vector<Enzyme> enzymes;
    std::vector<std::string> enzyme_names;
    for (const std::string &spec : enzyme_specs) {
        try {
            Enzyme enzyme = parse_enzyme_spec(spec, tables);
            enzymes.push_back(enzyme);
            enzyme_names.push_back(enzyme.name);
        } catch (const std::exception &e) {
            std::cerr << "Error parsing enzyme spec '" << spec << "': " << e.what() << std::endl;
            return 1;
        }
    }
    std::string joined_names;
    for (size_t i = 0; i < enzyme_names.size(); i++) {
        if (i) { joined_names += ", "; }
        joined_names += enzyme_names[i];
    }

    for (const std::pair<std::string, std::string> &entry : sequences) {
        const std::string &header = entry.first;
        int gaps_removed = 0;
        std::string seq = strip_gaps(entry.second, &gaps_removed);
        const int seq_len = (int) seq.size();

        std::cout << "\n=== Digest of " << header << " ===" << std::endl;
        std::cout << "Sequence length: " << seq_len << " bp, "
                  << (circular ? "circular" : "linear");
        if (gaps_removed) {
            std::cout << " (" << gaps_removed << " gap character(s) removed)";
        }
        std::cout << std::endl;
        std::cout << "Enzymes: " << joined_names << std::endl;

        // Only worth a second scan when the sequence actually carries codes.
        bool seq_is_ambiguous = false;
        for (char c : seq) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { seq_is_ambiguous = true; break; }
        }

        std::vector<int> all_cuts;
        std::vector<int> uncertain_cuts;
        for (const Enzyme &enzyme : enzymes) {
            int dropped = 0;
            std::vector<SiteHit> def_hits = find_sites(seq, enzyme, AMBIGUITY_DEFINITE, circular);
            std::vector<int> def_raw;
            for (const SiteHit &h : def_hits) { def_raw.push_back(h.cut); }
            std::vector<int> def_cuts = normalise_cuts(def_raw, seq_len, circular, &dropped);

            std::vector<SiteHit> pos_hits;
            std::vector<int> pos_cuts;
            int pos_dropped = 0;
            if (seq_is_ambiguous) {
                pos_hits = find_sites(seq, enzyme, AMBIGUITY_POSSIBLE, circular);
                std::vector<int> pos_raw;
                for (const SiteHit &h : pos_hits) { pos_raw.push_back(h.cut); }
                pos_cuts = normalise_cuts(pos_raw, seq_len, circular, &pos_dropped);
            } else {
                pos_hits = def_hits;
                pos_cuts = def_cuts;
                pos_dropped = dropped;
            }
            std::set<int> definite_set(def_cuts.begin(), def_cuts.end());
            std::vector<int> possible_only;
            for (int c : pos_cuts) {
                if (!definite_set.count(c)) { possible_only.push_back(c); }
            }

            std::string label = enzyme.name + " (" + to_neb_notation(enzyme) + ")";
            const std::vector<SiteHit> *hits = NULL;
            std::vector<int> cuts;

            if (ambiguity == AMBIGUITY_POSSIBLE) {
                hits = &pos_hits;
                cuts = pos_cuts;
                dropped = pos_dropped;
                std::set<int> speculative(possible_only.begin(), possible_only.end());
                if (!cuts.empty()) {
                    std::string rendered;
                    for (size_t i = 0; i < cuts.size(); i++) {
                        if (i) { rendered += ", "; }
                        rendered += to_string_int(cuts[i]);
                        if (speculative.count(cuts[i])) { rendered += "?"; }
                    }
                    std::cout << "  " << label << " cut sites at: [" << rendered << "]" << std::endl;
                } else if (!hits->empty()) {
                    std::cout << "  " << label << ": " << hits->size()
                              << " recognition site(s) found, none cut" << std::endl;
                } else {
                    // The origin wrap only runs under --circular, so a linear
                    // scan cannot promise anything about a circular molecule.
                    std::string scope = circular ? "" :
                        " (linear scan; pass --circular if this molecule is circular)";
                    std::cout << "  " << label
                              << ": no cut possible under any resolution of ambiguous bases"
                              << scope << std::endl;
                }
                uncertain_cuts.insert(uncertain_cuts.end(), possible_only.begin(), possible_only.end());
            } else {
                hits = &def_hits;
                cuts = def_cuts;
                if (!cuts.empty()) {
                    std::cout << "  " << label << " cut sites at: " << render_int_list(cuts)
                              << std::endl;
                } else if (!hits->empty()) {
                    std::cout << "  " << label << ": " << hits->size()
                              << " recognition site(s) found, none cut" << std::endl;
                } else {
                    std::cout << "  " << label << ": no cut sites found" << std::endl;
                }
                if (!possible_only.empty()) {
                    std::cout << "      note: " << possible_only.size()
                              << " further site(s) could cut depending on how ambiguous bases "
                              << "resolve - rerun with --ambiguity possible to include them"
                              << std::endl;
                }
            }

            if (!hits->empty() && !is_symmetric(enzyme)) {
                size_t bottom = 0;
                for (const SiteHit &h : *hits) {
                    if (h.strand == '-') { bottom++; }
                }
                std::cout << "    " << (hits->size() - bottom) << " site(s) on the top strand, "
                          << bottom << " on the bottom strand" << std::endl;
            }
            if (dropped) {
                std::cout << "    " << dropped << " site(s) recognised but cutting beyond the end "
                          << "of the linear sequence - not cut" << std::endl;
            }
            all_cuts.insert(all_cuts.end(), cuts.begin(), cuts.end());
        }

        if (!uncertain_cuts.empty()) {
            std::cout << "  (? = possible only: cut depends on how ambiguous bases resolve)"
                      << std::endl;
        }

        std::set<int> unique_cuts(all_cuts.begin(), all_cuts.end());
        std::vector<int> final_cuts(unique_cuts.begin(), unique_cuts.end());

        std::string frag_note;
        if (!uncertain_cuts.empty()) {
            frag_note = "(maximal-cut scenario: includes cuts that depend on ambiguous bases, "
                        "so these sizes may not match any single real sequence)";
        }

        std::vector<int> fragments = circular
            ? digest_circular(seq_len, final_cuts, min_fragment)
            : digest_linear(seq_len, final_cuts, min_fragment);

        if (output == "table" || output == "both") {
            print_fragment_table(fragments, frag_note);
        }
        if (output == "gel" || output == "both") {
            if (!frag_note.empty()) { std::cout << frag_note << std::endl; }
            std::cout << "\n--- Simulated Gel Electrophoresis (100 bp ladder) ---" << std::endl;
            draw_ascii_gel(fragments, DNA_LADDER_100BP, DNA_LADDER_100BP_COUNT, gel_height);
        }
    }
    return 0;
}
