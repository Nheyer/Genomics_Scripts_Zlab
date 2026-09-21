// SPDX-License-Identifier: GPL-3.0-only
//
// FASTA reader shared by all three tools.
//
// It grew out of the reader in MSA-to-consensus and replaces the separate ones
// Enzyme-digest and PCR-protocol used to carry. Nothing here was taken from
// Enzyme-digest's parse_fasta, so the MIT notice in that file does not apply to
// this one - see the License section of the README.
//
// Header only, and deliberately so: each tool is a single translation unit that
// owns its main() and the tests #include it whole, so a shared piece cannot be a
// second .cpp. Keep this free of anything tool specific.
//
// What it does is the part every caller agreed on: lines are trimmed (which is
// what makes CRLF files work), blank lines are skipped, a '>' line opens a
// record, and every other line is appended to it. What it does NOT do is decide
// what a valid sequence is. MSA-to-consensus takes protein, Enzyme-digest takes
// ACGT or IUPAC depending on a flag, PCR-protocol checks its own primers, so
// each of them supplies that itself, through options::check_line if it wants
// the offending line number.
#ifndef FASTA_HPP
#define FASTA_HPP

#include <cctype>
#include <fstream>
#include <functional>
#include <istream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fasta {

struct entry {
    std::string name;         // the header up to its first space
    std::string description;  // everything after that space; empty if there is none
    std::string seq;          // every sequence line, joined
    // The whole header, without the '>' or the whitespace around it. This is
    // name + ' ' + description whenever there is a description, and exactly what
    // the callers that want one string (Enzyme-digest, PCR-protocol) print. Last
    // so that anything brace-initialising name, description, seq still works.
    std::string header;
    void clear() {
        name.clear();
        description.clear();
        seq.clear();
        header.clear();
    }
};

// The ways a read can fail. Callers map these onto their own messages, since
// each tool reports errors in its own voice and exit convention.
enum class problem {
    cannot_open,             // the file could not be opened
    sequence_before_header,  // a sequence line with no '>' line above it
    no_records               // no '>' line at all - covers empty and blank-only input
};

class error : public std::runtime_error {
public:
    error(problem kind, const std::string &what) : std::runtime_error(what), kind_(kind) {}
    problem kind() const { return kind_; }
private:
    problem kind_;
};

struct options {
    // Upper-case every sequence line before it is checked or stored.
    bool uppercase = false;
    // Called for each sequence line once it is trimmed (and upper-cased, if
    // asked), with the record it is about to join and its 1-based line number in
    // the file, blank lines counted. Runs before the line is appended, so
    // whatever it does (throw, print and exit) leaves the record untouched.
    std::function<void(const entry &current, const std::string &line, int lineno)> check_line;
};

// Whitespace at either end. A plain space, tab, CR or LF: the CR is what a
// Windows line ending leaves behind, and it is the only reason this exists.
inline std::string trim(const std::string &s) {
    const char *ws = " \t\r\n";
    size_t first = s.find_first_not_of(ws);
    if (first == std::string::npos) { return ""; }
    size_t last = s.find_last_not_of(ws);
    return s.substr(first, last - first + 1);
}

// Read every record from a stream. *lines_out, if given, receives the number of
// lines read, blanks included; it is only set on success.
inline std::vector<entry> read(std::istream &in, const options &opt = options(),
                               int *lines_out = nullptr) {
    std::vector<entry> entries;
    entry current;
    bool open = false;   // a '>' line has been seen, so `current` is being filled
    std::string raw;
    int lineno = 0;
    while (std::getline(in, raw)) {
        ++lineno;
        std::string text = trim(raw);
        if (text.empty()) { continue; }
        if (text[0] == '>') {
            if (open) { entries.push_back(std::move(current)); }
            current.clear();   // also resets what the move left behind
            current.header = trim(text.substr(1));
            size_t space = current.header.find(' ');
            if (space == std::string::npos) {
                current.name = current.header;
            } else {
                current.name = current.header.substr(0, space);
                current.description = current.header.substr(space + 1);
            }
            open = true;
        } else {
            if (!open) {
                throw error(problem::sequence_before_header,
                            "line " + std::to_string(lineno) +
                            ": sequence data before any '>' header line");
            }
            if (opt.uppercase) {
                for (char &c : text) { c = (char) std::toupper((unsigned char) c); }
            }
            if (opt.check_line) { opt.check_line(current, text, lineno); }
            current.seq += text;
        }
    }
    // Records are only pushed when the next one starts, so the last one is still
    // waiting. Nothing open means no '>' line was ever seen.
    if (!open) {
        throw error(problem::no_records, "no '>' header line found");
    }
    entries.push_back(std::move(current));
    if (lines_out) { *lines_out = lineno; }
    return entries;
}

inline std::vector<entry> read_file(const std::string &path, const options &opt = options(),
                                    int *lines_out = nullptr) {
    std::ifstream in(path.c_str());
    if (!in) { throw error(problem::cannot_open, "cannot open " + path); }
    return read(in, opt, lines_out);
}

// (header, sequence) pairs, for the callers that only want those two strings.
inline std::vector<std::pair<std::string, std::string> > header_and_seq(
        const std::vector<entry> &entries) {
    std::vector<std::pair<std::string, std::string> > out;
    out.reserve(entries.size());
    for (const entry &e : entries) { out.push_back(std::make_pair(e.header, e.seq)); }
    return out;
}

}  // namespace fasta

#endif  // FASTA_HPP
