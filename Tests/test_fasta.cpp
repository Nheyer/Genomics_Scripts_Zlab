// SPDX-License-Identifier: GPL-3.0-only
//
// Unit tests for the shared FASTA reader, CppSrc/common/fasta.hpp.
//
// Unlike the other suites this one has no main() to rename away: the reader is a
// header, so it is included the way a tool includes it.
//
// Every expectation is written out by hand from the FASTA rules, not by running
// the reader and copying what it said. The two real files at the end are pinned
// to numbers that come from outside the code: a grep count of the records in one,
// and for the other the README's statement of which bases of phage lambda it holds.
#include "fasta.hpp"

#include <CUnit/Basic.h>

#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// ---------------------------------------------------------------- helpers --

static std::vector<fasta::entry> parse(const std::string &text,
                                       const fasta::options &opt = fasta::options()) {
    std::istringstream in(text);
    return fasta::read(in, opt);
}

// True if reading `text` fails with exactly this problem.
static bool fails_with(const std::string &text, fasta::problem expected,
                       const fasta::options &opt = fasta::options()) {
    try {
        parse(text, opt);
    } catch (const fasta::error &e) {
        return e.kind() == expected;
    }
    return false;
}

// --------------------------------------------------------------- records ---

static void test_two_records(void) {
    std::vector<fasta::entry> r = parse(">a\nAC\n>b\nGT\n");
    CU_ASSERT_EQUAL(r.size(), 2u);
    CU_ASSERT_EQUAL(r[0].name, std::string("a"));
    CU_ASSERT_EQUAL(r[0].seq, std::string("AC"));
    CU_ASSERT_EQUAL(r[1].name, std::string("b"));
    CU_ASSERT_EQUAL(r[1].seq, std::string("GT"));
}

static void test_sequence_lines_are_joined(void) {
    std::vector<fasta::entry> r = parse(">a\nAC\nGT\nTT\n");
    CU_ASSERT_EQUAL(r.size(), 1u);
    CU_ASSERT_EQUAL(r[0].seq, std::string("ACGTTT"));
}

static void test_name_stops_at_the_first_space(void) {
    std::vector<fasta::entry> r = parse(">seq1 first second\nA\n");
    CU_ASSERT_EQUAL(r[0].name, std::string("seq1"));
    CU_ASSERT_EQUAL(r[0].description, std::string("first second"));
    CU_ASSERT_EQUAL(r[0].header, std::string("seq1 first second"));
}

// The reader this grew out of asked substr for `space` characters where it meant
// "up to the space", and so kept the space on the end of every name that had a
// description.
static void test_name_has_no_trailing_space(void) {
    std::vector<fasta::entry> r = parse(">seq1 described\nA\n");
    CU_ASSERT_EQUAL(r[0].name.size(), 4u);
    CU_ASSERT_EQUAL(r[0].name, std::string("seq1"));
}

static void test_no_description(void) {
    std::vector<fasta::entry> r = parse(">solo\nA\n");
    CU_ASSERT_EQUAL(r[0].name, std::string("solo"));
    CU_ASSERT_EQUAL(r[0].description, std::string(""));
    CU_ASSERT_EQUAL(r[0].header, std::string("solo"));
}

// Only the first space splits, so a second one is part of the description and
// name + ' ' + description still gives the header back exactly.
static void test_extra_spaces_stay_in_the_description(void) {
    std::vector<fasta::entry> r = parse(">a  b\nA\n");
    CU_ASSERT_EQUAL(r[0].name, std::string("a"));
    CU_ASSERT_EQUAL(r[0].description, std::string(" b"));
    CU_ASSERT_EQUAL(r[0].header, std::string("a  b"));
}

// Only a space separates the name; a tab does not. Nobody split on tabs before.
static void test_a_tab_does_not_split_the_name(void) {
    std::vector<fasta::entry> r = parse(">a\tb c\nA\n");
    CU_ASSERT_EQUAL(r[0].name, std::string("a\tb"));
    CU_ASSERT_EQUAL(r[0].description, std::string("c"));
}

static void test_header_is_trimmed(void) {
    std::vector<fasta::entry> r = parse("  >  spaced out  \nA\n");
    CU_ASSERT_EQUAL(r[0].header, std::string("spaced out"));
    CU_ASSERT_EQUAL(r[0].name, std::string("spaced"));
}

static void test_header_only_record_has_an_empty_sequence(void) {
    std::vector<fasta::entry> r = parse(">only\n");
    CU_ASSERT_EQUAL(r.size(), 1u);
    CU_ASSERT_EQUAL(r[0].seq, std::string(""));
}

static void test_empty_record_between_records(void) {
    std::vector<fasta::entry> r = parse(">1\n>2\nA\n");
    CU_ASSERT_EQUAL(r.size(), 2u);
    CU_ASSERT_EQUAL(r[0].seq, std::string(""));
    CU_ASSERT_EQUAL(r[1].seq, std::string("A"));
}

static void test_bare_gt_is_a_record_with_no_name(void) {
    std::vector<fasta::entry> r = parse(">\nA\n");
    CU_ASSERT_EQUAL(r.size(), 1u);
    CU_ASSERT_EQUAL(r[0].header, std::string(""));
    CU_ASSERT_EQUAL(r[0].seq, std::string("A"));
}

// ----------------------------------------------------------------- lines ---

// The reason lines are trimmed: a Windows file leaves a CR on every one.
static void test_crlf_line_endings(void) {
    std::vector<fasta::entry> r = parse(">a x\r\nAC\r\nGT\r\n");
    CU_ASSERT_EQUAL(r[0].header, std::string("a x"));
    CU_ASSERT_EQUAL(r[0].seq, std::string("ACGT"));
}

static void test_blank_lines_are_skipped(void) {
    std::vector<fasta::entry> r = parse("\n>a\n\nAC\n \n\t\nGT\n\n");
    CU_ASSERT_EQUAL(r.size(), 1u);
    CU_ASSERT_EQUAL(r[0].seq, std::string("ACGT"));
}

static void test_sequence_lines_are_trimmed(void) {
    std::vector<fasta::entry> r = parse(">a\n  AC \t\n\tGT\n");
    CU_ASSERT_EQUAL(r[0].seq, std::string("ACGT"));
}

static void test_no_newline_at_the_end(void) {
    std::vector<fasta::entry> r = parse(">a\nAC");
    CU_ASSERT_EQUAL(r[0].seq, std::string("AC"));
}

// Whether a character is valid is the caller's call, not the reader's.
static void test_whitespace_inside_a_line_is_kept(void) {
    std::vector<fasta::entry> r = parse(">a\nAC GT\n");
    CU_ASSERT_EQUAL(r[0].seq, std::string("AC GT"));
}

// Space, tab, CR and LF are trimmed and nothing else. Enzyme-digest's output
// depends on that: a form feed has to reach its alphabet check and be refused.
static void test_only_space_tab_cr_lf_are_trimmed(void) {
    std::vector<fasta::entry> r = parse(">a\n\fAC\v\n");
    CU_ASSERT_EQUAL(r[0].seq, std::string("\fAC\v"));
}

// ------------------------------------------------------------------ case ---

static void test_case_is_kept_by_default(void) {
    std::vector<fasta::entry> r = parse(">a\nacGT\n");
    CU_ASSERT_EQUAL(r[0].seq, std::string("acGT"));
}

static void test_uppercase_option(void) {
    fasta::options opt;
    opt.uppercase = true;
    std::vector<fasta::entry> r = parse(">Mixed Case\nacGt\n", opt);
    CU_ASSERT_EQUAL(r[0].seq, std::string("ACGT"));
    CU_ASSERT_EQUAL(r[0].header, std::string("Mixed Case"));   // headers are never touched
}

// ---------------------------------------------------------------- errors ---

static void test_sequence_before_header_is_refused(void) {
    CU_ASSERT_TRUE(fails_with("AC\n>a\nGT\n", fasta::problem::sequence_before_header));
    CU_ASSERT_TRUE(fails_with("ACGT\n", fasta::problem::sequence_before_header));
}

static void test_empty_input_has_no_records(void) {
    CU_ASSERT_TRUE(fails_with("", fasta::problem::no_records));
}

static void test_blank_only_input_has_no_records(void) {
    CU_ASSERT_TRUE(fails_with("\n \n\t\n", fasta::problem::no_records));
}

static void test_missing_file(void) {
    bool refused = false;
    try {
        fasta::read_file("/no/such/directory/x.fna");
    } catch (const fasta::error &e) {
        refused = (e.kind() == fasta::problem::cannot_open);
    }
    CU_ASSERT_TRUE(refused);
}

// The message names the line, so a person can go and look at it.
static void test_the_message_names_the_line(void) {
    std::string message;
    try {
        parse("\n\nAC\n");
    } catch (const fasta::error &e) {
        message = e.what();
    }
    CU_ASSERT_TRUE(message.find("line 3") != std::string::npos);
}

// ------------------------------------------------------------ check_line ---

struct seen_line {
    int lineno;
    std::string line;
    std::string header;
    bool operator==(const seen_line &o) const {
        return lineno == o.lineno && line == o.line && header == o.header;
    }
};

// Blank lines count towards the number, and headers are not sequence lines.
static void test_check_line_gets_the_trimmed_line_and_its_number(void) {
    std::vector<seen_line> seen;
    fasta::options opt;
    opt.check_line = [&](const fasta::entry &current, const std::string &line, int lineno) {
        seen.push_back({lineno, line, current.header});
    };
    parse(">a\n\n  AC \nGT\n>b\nTT\n", opt);
    std::vector<seen_line> expected = {{3, "AC", "a"}, {4, "GT", "a"}, {6, "TT", "b"}};
    CU_ASSERT_TRUE(seen == expected);
}

static void test_check_line_sees_the_uppercased_text(void) {
    std::string got;
    fasta::options opt;
    opt.uppercase = true;
    opt.check_line = [&](const fasta::entry &, const std::string &line, int) { got = line; };
    parse(">a\nacgt\n", opt);
    CU_ASSERT_EQUAL(got, std::string("ACGT"));
}

// It runs before the line joins the record, so what it is shown is what came
// before.
static void test_check_line_runs_before_the_line_is_appended(void) {
    std::vector<std::string> before;
    fasta::options opt;
    opt.check_line = [&](const fasta::entry &current, const std::string &, int) {
        before.push_back(current.seq);
    };
    parse(">a\nAC\nGT\nTT\n", opt);
    std::vector<std::string> expected = {"", "AC", "ACGT"};
    CU_ASSERT_TRUE(before == expected);
}

static void test_check_line_can_stop_the_read(void) {
    int calls = 0;
    bool propagated = false;
    fasta::options opt;
    opt.check_line = [&](const fasta::entry &, const std::string &line, int) {
        calls++;
        if (line == "GT") { throw std::invalid_argument("bad line"); }
    };
    try {
        parse(">a\nAC\nGT\nTT\nCC\n", opt);
    } catch (const std::invalid_argument &) {
        propagated = true;
    }
    CU_ASSERT_TRUE(propagated);
    CU_ASSERT_EQUAL(calls, 2);   // AC, then GT; nothing after it was read
}

// It is the reader's own rule that comes first: a sequence line with no header
// is refused before the caller is asked about it.
static void test_check_line_is_not_asked_about_a_headerless_line(void) {
    int calls = 0;
    fasta::options opt;
    opt.check_line = [&](const fasta::entry &, const std::string &, int) { calls++; };
    CU_ASSERT_TRUE(fails_with("AC\n", fasta::problem::sequence_before_header, opt));
    CU_ASSERT_EQUAL(calls, 0);
}

// ---------------------------------------------------------------- counts ---

static void test_line_count_includes_blanks(void) {
    std::istringstream in("\n>a\n\nAC\n");
    int lines = 0;
    fasta::read(in, fasta::options(), &lines);
    CU_ASSERT_EQUAL(lines, 4);
}

static void test_line_count_is_left_alone_on_failure(void) {
    std::istringstream in("AC\n");
    int lines = -7;
    try {
        fasta::read(in, fasta::options(), &lines);
    } catch (const fasta::error &) {
    }
    CU_ASSERT_EQUAL(lines, -7);
}

static void test_header_and_seq_pairs(void) {
    std::vector<std::pair<std::string, std::string> > p =
        fasta::header_and_seq(parse(">a x\nAC\n>b\nGT\n"));
    CU_ASSERT_EQUAL(p.size(), 2u);
    CU_ASSERT_EQUAL(p[0].first, std::string("a x"));
    CU_ASSERT_EQUAL(p[0].second, std::string("AC"));
    CU_ASSERT_EQUAL(p[1].first, std::string("b"));
    CU_ASSERT_EQUAL(p[1].second, std::string("GT"));
}

// ------------------------------------------------------------ real files ---

// 1000 records, counted with grep on the '>' lines rather than by this reader.
static void test_sheep_products_has_1000_records(void) {
    std::vector<fasta::entry> r = fasta::read_file(FASTA_TEST_SHEEP);
    CU_ASSERT_EQUAL(r.size(), 1000u);
}

// README: "bases 29,951-40,100 of phage lambda", which is 10,150 of them, all
// unambiguous. Wrapped over many lines, so this is the joining checked end to end.
static void test_lambda_fragment_is_10150_unambiguous_bases(void) {
    std::vector<fasta::entry> r = fasta::read_file(FASTA_TEST_LAMBDA);
    CU_ASSERT_EQUAL(r.size(), 1u);
    CU_ASSERT_EQUAL(r[0].seq.size(), 10150u);
    CU_ASSERT_EQUAL(r[0].seq.find_first_not_of("ACGT"), std::string::npos);
}

// -------------------------------------------------------------- registry ---

struct test_case {
    const char *name;
    void (*fn)(void);
};

struct test_suite {
    const char *name;
    std::initializer_list<test_case> cases;
};

int main(void) {
    if (CU_initialize_registry() != CUE_SUCCESS) { return CU_get_error(); }

    const std::initializer_list<test_suite> suites = {
        {"records", {
            {"two records", test_two_records},
            {"sequence lines are joined", test_sequence_lines_are_joined},
            {"name stops at the first space", test_name_stops_at_the_first_space},
            {"name has no trailing space", test_name_has_no_trailing_space},
            {"no description", test_no_description},
            {"extra spaces stay in the description", test_extra_spaces_stay_in_the_description},
            {"a tab does not split the name", test_a_tab_does_not_split_the_name},
            {"header is trimmed", test_header_is_trimmed},
            {"header only record", test_header_only_record_has_an_empty_sequence},
            {"empty record between records", test_empty_record_between_records},
            {"bare > is a record", test_bare_gt_is_a_record_with_no_name},
        }},
        {"lines", {
            {"CRLF line endings", test_crlf_line_endings},
            {"blank lines skipped", test_blank_lines_are_skipped},
            {"sequence lines trimmed", test_sequence_lines_are_trimmed},
            {"no newline at the end", test_no_newline_at_the_end},
            {"inner whitespace kept", test_whitespace_inside_a_line_is_kept},
            {"only space tab CR LF trimmed", test_only_space_tab_cr_lf_are_trimmed},
        }},
        {"case", {
            {"kept by default", test_case_is_kept_by_default},
            {"uppercase option", test_uppercase_option},
        }},
        {"errors", {
            {"sequence before header", test_sequence_before_header_is_refused},
            {"empty input", test_empty_input_has_no_records},
            {"blank only input", test_blank_only_input_has_no_records},
            {"missing file", test_missing_file},
            {"message names the line", test_the_message_names_the_line},
        }},
        {"check_line", {
            {"trimmed line and its number", test_check_line_gets_the_trimmed_line_and_its_number},
            {"sees uppercased text", test_check_line_sees_the_uppercased_text},
            {"runs before the append", test_check_line_runs_before_the_line_is_appended},
            {"can stop the read", test_check_line_can_stop_the_read},
            {"not asked about a headerless line", test_check_line_is_not_asked_about_a_headerless_line},
        }},
        {"counts and helpers", {
            {"line count includes blanks", test_line_count_includes_blanks},
            {"line count untouched on failure", test_line_count_is_left_alone_on_failure},
            {"header and seq pairs", test_header_and_seq_pairs},
        }},
        {"real files", {
            {"Sheep_products has 1000 records", test_sheep_products_has_1000_records},
            {"lambda fragment is 10150 bases", test_lambda_fragment_is_10150_unambiguous_bases},
        }},
    };

    for (const test_suite &s : suites) {
        CU_pSuite suite = CU_add_suite(s.name, NULL, NULL);
        if (suite == NULL) {
            CU_cleanup_registry();
            return CU_get_error();
        }
        for (const test_case &t : s.cases) {
            if (CU_add_test(suite, t.name, t.fn) == NULL) {
                CU_cleanup_registry();
                return CU_get_error();
            }
        }
    }

    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();

    // Exit non zero when anything failed, so the build or CI notices.
    unsigned int failed = CU_get_number_of_tests_failed();
    unsigned int failed_asserts = CU_get_number_of_failures();
    CU_cleanup_registry();
    return (failed > 0 || failed_asserts > 0) ? 1 : 0;
}
