// Unit tests for MSA_fasta_to_Consensus.
//
// The program is one translation unit that owns its own main(), so we rename
// that main out of the way and pull the whole file in. The tests then call the
// real functions directly, and Scripts/MSA_fasta_to_Consensus.cpp needs no
// changes to be testable.
#define main disabled_main
#include "MSA_fasta_to_Consensus.cpp"
#undef main

#include <CUnit/Basic.h>

#include <initializer_list>
#include <string>
#include <vector>

// ---------------------------------------------------------------- helpers --

// Build an alignment out of plain strings, one per sequence.
static std::vector<fasta_entry> alignment(std::initializer_list<const char *> seqs) {
    std::vector<fasta_entry> out;
    int i = 0;
    for (const char *s : seqs) {
        fasta_entry e;
        e.name = "seq" + std::to_string(i++);
        e.seq = s;
        out.push_back(e);
    }
    return out;
}

// An alignment of `count` identical rows, used for the compression tests.
static std::vector<fasta_entry> repeated(const char *seq, int count) {
    std::vector<fasta_entry> out;
    for (int i = 0; i < count; i++) {
        fasta_entry e;
        e.name = "seq" + std::to_string(i);
        e.seq = seq;
        out.push_back(e);
    }
    return out;
}

// Run make_consensus and hand back the consensus string. Returns "" if it
// failed, so a test that expects sequence never silently passes on an error.
static std::string consensus_of(const std::vector<fasta_entry> &in, bool strict) {
    std::vector<fasta_entry> out;
    if (make_consensus(in, &out, strict) != 0) { return ""; }
    if (out.size() != 1) { return ""; }
    return out[0].seq;
}

// Shorthand: the two modes we care about.
static std::string ambiguous(const std::vector<fasta_entry> &in) { return consensus_of(in, false); }
static std::string strict(const std::vector<fasta_entry> &in) { return consensus_of(in, true); }

// ------------------------------------------------------- clean_nucliotide --

static void test_clean_uppercases(void) {
    CU_ASSERT_EQUAL(clean_nucliotide('a'), 'A');
    CU_ASSERT_EQUAL(clean_nucliotide('c'), 'C');
    CU_ASSERT_EQUAL(clean_nucliotide('g'), 'G');
    CU_ASSERT_EQUAL(clean_nucliotide('t'), 'T');
}

static void test_clean_leaves_upper_alone(void) {
    CU_ASSERT_EQUAL(clean_nucliotide('A'), 'A');
    CU_ASSERT_EQUAL(clean_nucliotide('C'), 'C');
    CU_ASSERT_EQUAL(clean_nucliotide('G'), 'G');
    CU_ASSERT_EQUAL(clean_nucliotide('T'), 'T');
}

// RNA in, DNA out.
static void test_clean_maps_u_to_t(void) {
    CU_ASSERT_EQUAL(clean_nucliotide('U'), 'T');
    CU_ASSERT_EQUAL(clean_nucliotide('u'), 'T');
}

static void test_clean_passes_gap_through(void) {
    CU_ASSERT_EQUAL(clean_nucliotide('-'), '-');
}

static void test_clean_passes_ambiguity_codes_through(void) {
    CU_ASSERT_EQUAL(clean_nucliotide('r'), 'R');
    CU_ASSERT_EQUAL(clean_nucliotide('n'), 'N');
    CU_ASSERT_EQUAL(clean_nucliotide('V'), 'V');
}

// ------------------------------------------------------ encode_nucliotide --

static void test_encode_bases(void) {
    CU_ASSERT_EQUAL(encode_nucliotide('A'), 3);
    CU_ASSERT_EQUAL(encode_nucliotide('T'), 5);
    CU_ASSERT_EQUAL(encode_nucliotide('C'), 7);
    CU_ASSERT_EQUAL(encode_nucliotide('G'), 11);
}

// Every ambiguity code is the product of the bases it stands for.
static void test_encode_two_base_codes(void) {
    CU_ASSERT_EQUAL(encode_nucliotide('R'), 3 * 11);  // A or G
    CU_ASSERT_EQUAL(encode_nucliotide('Y'), 7 * 5);   // C or T
    CU_ASSERT_EQUAL(encode_nucliotide('K'), 11 * 5);  // G or T
    CU_ASSERT_EQUAL(encode_nucliotide('M'), 3 * 7);   // A or C
    CU_ASSERT_EQUAL(encode_nucliotide('S'), 7 * 11);  // C or G
    CU_ASSERT_EQUAL(encode_nucliotide('W'), 3 * 5);   // A or T
}

static void test_encode_three_base_codes(void) {
    CU_ASSERT_EQUAL(encode_nucliotide('B'), 7 * 11 * 5);  // not A
    CU_ASSERT_EQUAL(encode_nucliotide('D'), 3 * 11 * 5);  // not C
    CU_ASSERT_EQUAL(encode_nucliotide('H'), 3 * 7 * 5);   // not G
    CU_ASSERT_EQUAL(encode_nucliotide('V'), 3 * 7 * 11);  // not T
}

static void test_encode_n_is_every_base(void) {
    CU_ASSERT_EQUAL(encode_nucliotide('N'), 3 * 7 * 5 * 11);
}

// ------------------------------------------------------ decode_nucliotide --

// Every code the encoder produces must come back as the same letter.
static void test_decode_round_trips_every_code(void) {
    const char *codes = "ATCGRYKMSWBDHVN";
    for (const char *c = codes; *c; c++) {
        CU_ASSERT_EQUAL(decode_nucliotide(encode_nucliotide(*c)), *c);
    }
}

// Decoding a product of single bases gives the code covering that set.
static void test_decode_combines_bases(void) {
    CU_ASSERT_EQUAL(decode_nucliotide(nuc_A * nuc_G), 'R');
    CU_ASSERT_EQUAL(decode_nucliotide(nuc_C * nuc_T), 'Y');
    CU_ASSERT_EQUAL(decode_nucliotide(nuc_A * nuc_C * nuc_G), 'V');
    CU_ASSERT_EQUAL(decode_nucliotide(nuc_A * nuc_C * nuc_G * nuc_T), 'N');
}

// Repeating a base does not widen the code: A*A is still just A.
static void test_decode_ignores_repeats(void) {
    CU_ASSERT_EQUAL(decode_nucliotide(nuc_A * nuc_A), 'A');
    CU_ASSERT_EQUAL(decode_nucliotide(nuc_A * nuc_G * nuc_G), 'R');
}

// --------------------------------------------- make_consensus, agreement ---

static void test_identical_sequences_are_unchanged(void) {
    auto in = alignment({"ACGT", "ACGT", "ACGT"});
    CU_ASSERT_STRING_EQUAL(strict(in).c_str(), "ACGT");
    CU_ASSERT_STRING_EQUAL(ambiguous(in).c_str(), "ACGT");
}

// U in the input is normalised to T on the way out.
static void test_rna_input_becomes_dna(void) {
    auto in = alignment({"ACGU", "ACGU"});
    CU_ASSERT_STRING_EQUAL(strict(in).c_str(), "ACGT");
    CU_ASSERT_STRING_EQUAL(ambiguous(in).c_str(), "ACGT");
}

// ------------------------------------------- make_consensus, ambiguity ----

// One column per IUPAC code, same known truth as the test_files fixtures.
static void test_two_base_ambiguity_codes(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "G"})).c_str(), "R");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"C", "T"})).c_str(), "Y");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"G", "C"})).c_str(), "S");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "T"})).c_str(), "W");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"G", "T"})).c_str(), "K");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "C"})).c_str(), "M");
}

static void test_three_base_ambiguity_codes(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"C", "G", "T"})).c_str(), "B");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "G", "T"})).c_str(), "D");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "C", "T"})).c_str(), "H");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "C", "G"})).c_str(), "V");
}

static void test_all_four_bases_give_n(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "C", "G", "T"})).c_str(), "N");
}

// Row order must not change the answer.
static void test_ambiguity_is_order_independent(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "G"})).c_str(), "R");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"G", "A"})).c_str(), "R");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"T", "A", "C"})).c_str(), "H");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"C", "T", "A"})).c_str(), "H");
}

// An ambiguity code already in the alignment is expanded and merged, not
// treated as an unknown letter: {A,G} + {C} = {A,C,G} = V.
static void test_ambiguity_code_in_input_is_merged(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"R", "C"})).c_str(), "V");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"Y", "A"})).c_str(), "H");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"R", "Y"})).c_str(), "N");
}

// A literal N is already every base, so it swallows whatever it meets.
static void test_literal_n_absorbs_the_column(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "G", "N"})).c_str(), "N");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"N", "A", "G"})).c_str(), "N");
}

// ---------------------------------------------- make_consensus, strict ----

// Strict mode does not reason about ambiguity, any disagreement is masked.
static void test_strict_masks_disagreement(void) {
    CU_ASSERT_STRING_EQUAL(strict(alignment({"A", "G"})).c_str(), "N");
    CU_ASSERT_STRING_EQUAL(strict(alignment({"C", "T"})).c_str(), "N");
    CU_ASSERT_STRING_EQUAL(strict(alignment({"A", "C", "G"})).c_str(), "N");
}

// Agreeing columns survive, only the disagreeing one is masked.
static void test_strict_masks_only_the_bad_column(void) {
    CU_ASSERT_STRING_EQUAL(strict(alignment({"ACGT", "ACTT"})).c_str(), "ACNT");
}

static void test_ambiguity_only_touches_the_bad_column(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"ACGT", "ACTT"})).c_str(), "ACKT");
}

// ------------------------------------------------- make_consensus, gaps ---

// A gap beats everything else in the column, in both modes and wherever it is.
static void test_gap_wins_in_any_row(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "-", "G"})).c_str(), "-");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"-", "A", "G"})).c_str(), "-");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"A", "G", "-"})).c_str(), "-");
}

static void test_gap_wins_in_strict_mode_too(void) {
    CU_ASSERT_STRING_EQUAL(strict(alignment({"A", "-", "G"})).c_str(), "-");
    CU_ASSERT_STRING_EQUAL(strict(alignment({"-", "A", "G"})).c_str(), "-");
}

// The gap must not leak into columns either side of it.
static void test_gap_is_confined_to_its_column(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"AAA", "A-A", "AAA"})).c_str(), "A-A");
    CU_ASSERT_STRING_EQUAL(strict(alignment({"AAA", "A-A", "AAA"})).c_str(), "A-A");
}

// A gap beats an N as well, N does not win the column.
static void test_gap_beats_n(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"N", "-"})).c_str(), "-");
    CU_ASSERT_STRING_EQUAL(ambiguous(alignment({"-", "N"})).c_str(), "-");
}

// ------------------------------------------ make_consensus, compression ---
//
// The accumulator multiplies a prime per sequence, so it eventually runs out of
// headroom in an unsigned long long and the code squarefree-reduces it back to
// the distinct primes seen so far. Measured against this build, an all-G column
// stays under the limit through 16 rows and first compresses at 17.
//
// Whatever the row count, the answer must not drift: every row says G, so the
// consensus is G.

static void test_just_under_compression_threshold(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(repeated("G", 16)).c_str(), "G");
}

static void test_at_compression_threshold(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(repeated("G", 17)).c_str(), "G");
}

static void test_far_past_compression_threshold(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(repeated("G", 100)).c_str(), "G");
    CU_ASSERT_STRING_EQUAL(ambiguous(repeated("G", 1000)).c_str(), "G");
}

// Compression must not lose bases that were already accumulated: these columns
// disagree *and* run long, so the answer is the wide code, not a narrow one.
static void test_compression_keeps_every_base_seen(void) {
    std::vector<fasta_entry> in = repeated("A", 30);
    in.back().seq = "G";  // 29 A's then one G
    CU_ASSERT_STRING_EQUAL(ambiguous(in).c_str(), "R");

    std::vector<fasta_entry> wide = repeated("A", 40);
    wide[10].seq = "C";
    wide[20].seq = "G";
    CU_ASSERT_STRING_EQUAL(ambiguous(wide).c_str(), "V");
}

// A long alignment that still agrees must stay that single base.
static void test_compression_on_agreeing_columns(void) {
    CU_ASSERT_STRING_EQUAL(ambiguous(repeated("ACGT", 50)).c_str(), "ACGT");
}

// Long alignments work in strict mode too (no accumulator involved).
static void test_long_alignment_strict(void) {
    CU_ASSERT_STRING_EQUAL(strict(repeated("ACGT", 50)).c_str(), "ACGT");
}

// ---------------------------------------------- make_consensus, errors ----

// A consensus of one sequence is meaningless, so it is refused.
static void test_single_sequence_is_refused(void) {
    std::vector<fasta_entry> out;
    auto in = alignment({"ACGT"});
    CU_ASSERT_EQUAL(make_consensus(in, &out, true), -1);
}

// The input is supposed to be aligned, so ragged lengths are refused.
static void test_mismatched_lengths_are_refused(void) {
    std::vector<fasta_entry> out;
    auto in = alignment({"ACGT", "ACG"});
    CU_ASSERT_EQUAL(make_consensus(in, &out, true), -1);
    CU_ASSERT_EQUAL(make_consensus(in, &out, false), -1);
}

// ------------------------------------------------------------- registry ---

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
        {"clean_nucliotide", {
            {"uppercases input", test_clean_uppercases},
            {"leaves uppercase alone", test_clean_leaves_upper_alone},
            {"maps U to T", test_clean_maps_u_to_t},
            {"passes gaps through", test_clean_passes_gap_through},
            {"passes ambiguity codes through", test_clean_passes_ambiguity_codes_through},
        }},
        {"encode_nucliotide", {
            {"bases are primes", test_encode_bases},
            {"two base codes", test_encode_two_base_codes},
            {"three base codes", test_encode_three_base_codes},
            {"N covers every base", test_encode_n_is_every_base},
        }},
        {"decode_nucliotide", {
            {"round trips every code", test_decode_round_trips_every_code},
            {"combines bases", test_decode_combines_bases},
            {"ignores repeats", test_decode_ignores_repeats},
        }},
        {"make_consensus: agreement", {
            {"identical sequences unchanged", test_identical_sequences_are_unchanged},
            {"RNA becomes DNA", test_rna_input_becomes_dna},
        }},
        {"make_consensus: ambiguity", {
            {"two base codes", test_two_base_ambiguity_codes},
            {"three base codes", test_three_base_ambiguity_codes},
            {"all four bases give N", test_all_four_bases_give_n},
            {"order independent", test_ambiguity_is_order_independent},
            {"merges codes already in input", test_ambiguity_code_in_input_is_merged},
            {"literal N absorbs the column", test_literal_n_absorbs_the_column},
            {"only the bad column changes", test_ambiguity_only_touches_the_bad_column},
        }},
        {"make_consensus: strict", {
            {"masks disagreement", test_strict_masks_disagreement},
            {"masks only the bad column", test_strict_masks_only_the_bad_column},
        }},
        {"make_consensus: gaps", {
            {"gap wins in any row", test_gap_wins_in_any_row},
            {"gap wins in strict mode", test_gap_wins_in_strict_mode_too},
            {"gap confined to its column", test_gap_is_confined_to_its_column},
            {"gap beats N", test_gap_beats_n},
        }},
        {"make_consensus: compression", {
            {"just under threshold (16 rows)", test_just_under_compression_threshold},
            {"at threshold (17 rows)", test_at_compression_threshold},
            {"far past threshold", test_far_past_compression_threshold},
            {"keeps every base seen", test_compression_keeps_every_base_seen},
            {"agreeing columns survive", test_compression_on_agreeing_columns},
            {"long alignment in strict mode", test_long_alignment_strict},
        }},
        {"make_consensus: errors", {
            {"single sequence refused", test_single_sequence_is_refused},
            {"mismatched lengths refused", test_mismatched_lengths_are_refused},
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
