// Unit tests for PCR_Protocol.
//
// Same arrangement as the other suites: the program is one translation unit
// that owns its own main(), so we rename that main out of the way and pull the
// whole file in.
//
// None of the expected numbers here come from running this code. Each is
// pinned to something outside it:
//   - the SantaLucia 1998 paper's own Table 1, against the Table 2 we compute
//     from (two tables, one paper, checked against each other),
//   - primer3's oligotm(), linked in and called side by side with ours, and
//     values from primer3-py 2.3.1 (a separate build of the same code),
//   - the worked examples in the primer3 manual,
//   - NEB's lambda control primers and the coordinates NEB gives for them,
//   - arithmetic done by hand in the comments.
#define main disabled_main
#include "PCR_Protocol.cpp"
#undef main

#include <CUnit/Basic.h>

#include "oligotm.h"

#include <initializer_list>
#include <string>
#include <vector>

#ifndef PCR_TEST_PRIMERS
#error "PCR_TEST_PRIMERS must point at Data/test_files/TEST_PCR_PRIMERS.fna"
#endif
#ifndef PCR_TEST_FIXTURE
#error "PCR_TEST_FIXTURE must point at Data/test_files/TEST_PCR_LAMBDA.fna"
#endif

// ---------------------------------------------------------------- helpers --

// NEB's control primers for lambda, from the Phusion High Fidelity PCR Kit
// (E0553) manual, with the genome coordinates NEB states for them.
static const char *const NEB_1_3KB_P1 = "GTCACCAGTGCAGTGCTTGATAACAGG";   // 30,006-30,032
static const char *const NEB_1_3KB_P2 = "GATGACGCATCCTCACGATAATATCCGG";  // 31,325-31,352
static const char *const NEB_10KB_P1 = "CAGTGCAGTGCTTGATAACAGG";         // 30,011-30,032
static const char *const NEB_10KB_P2 = "GTAGTGCGCGTTTGATTTCC";           // 40,024-40,043

// The fixture holds lambda bases 29,951..40,100, so genome position g is
// fixture position g - 29,950.
static const long FIXTURE_OFFSET = 29950;

static const std::vector<std::pair<std::string, std::string> > &lambda_fixture() {
    static std::vector<std::pair<std::string, std::string> > records =
        parse_fasta(PCR_TEST_FIXTURE);
    return records;
}

static const Polymerase &builtin(const char *name) {
    return find_polymerase(name, builtin_polymerases());
}

template <typename Fn>
static bool throws(Fn fn) {
    try {
        fn();
    } catch (const std::exception &) {
        return true;
    }
    return false;
}

static bool near(double a, double b, double tolerance) { return std::fabs(a - b) <= tolerance; }

// primer3's own Tm, SantaLucia stacks and salt correction, as it is built here.
static double primer3_tm(const std::string &seq, double mono, double mg, double dntp,
                         double primer_nm) {
    return oligotm(seq.c_str(), primer_nm, mono, mg, dntp, 0.0, 0.6, 0.0, santalucia_auto,
                   santalucia, -10.0).Tm;
}

// ----------------------------------------------------- nearest neighbour ---

// SantaLucia 1998 gives the same parameters twice: dH and dS in Table 2, and
// dG37 in the "Unified" column of Table 1. We compute from Table 2, so check
// every row of it against Table 1, typed in separately below. A transcription
// slip in either shows up as a mismatch. Table 1 rounds to 0.01, and dG from
// dH and dS rounded to 0.1 lands within 0.05 of it.
static void test_nn_table_agrees_with_the_papers_dg_column(void) {
    struct { const char *stack; double dg; } TABLE_1_UNIFIED[] = {
        {"AA", -1.00}, {"AT", -0.88}, {"TA", -0.58}, {"CA", -1.45}, {"GT", -1.44},
        {"CT", -1.28}, {"GA", -1.30}, {"CG", -2.17}, {"GC", -2.24}, {"GG", -1.84},
    };
    CU_ASSERT_EQUAL(NN_UNIFIED_COUNT, 10);
    for (const auto &row : TABLE_1_UNIFIED) {
        Thermo t = nn_stack(row.stack[0], row.stack[1]);
        double dg = t.dH - 310.15 * t.dS / 1000.0;
        CU_ASSERT_TRUE(near(dg, row.dg, 0.05));
    }
    // Initiation: Table 1 gives +0.98 (terminal G.C) and +1.03 (terminal A.T).
    CU_ASSERT_TRUE(near(INIT_GC_DH - 310.15 * INIT_GC_DS / 1000.0, 0.98, 0.05));
    CU_ASSERT_TRUE(near(INIT_AT_DH - 310.15 * INIT_AT_DS / 1000.0, 1.03, 0.05));
}

// The ten listed stacks cover all sixteen by reading the other strand.
static void test_unlisted_stacks_are_the_other_strand(void) {
    // 5'-TT-3' is the AA/TT duplex upside down, 5'-AC-3' is GT/CA.
    Thermo tt = nn_stack('T', 'T'), aa = nn_stack('A', 'A');
    CU_ASSERT_DOUBLE_EQUAL(tt.dH, aa.dH, 1e-12);
    CU_ASSERT_DOUBLE_EQUAL(tt.dS, aa.dS, 1e-12);
    Thermo ac = nn_stack('A', 'C'), gt = nn_stack('G', 'T');
    CU_ASSERT_DOUBLE_EQUAL(ac.dH, gt.dH, 1e-12);
    CU_ASSERT_DOUBLE_EQUAL(ac.dS, gt.dS, 1e-12);
    // And all sixteen resolve.
    const char *bases = "ACGT";
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            CU_ASSERT_FALSE(throws([&] { nn_stack(bases[i], bases[j]); }));
        }
    }
}

// Palindromes take the symmetry penalty and the C (not C/4) strand term;
// getting the test backwards would swap them silently.
static void test_self_complementary_is_pinned(void) {
    CU_ASSERT_TRUE(is_self_complementary("GAATTC"));   // EcoRI site
    CU_ASSERT_TRUE(is_self_complementary("ACGT"));
    CU_ASSERT_FALSE(is_self_complementary("GTCACC"));
    CU_ASSERT_FALSE(is_self_complementary("AAAA"));
}

// Hand arithmetic: 50 + 120 * sqrt(2.0 - 0.8) = 50 + 131.45 = 181.45 mM.
static void test_sodium_equivalents(void) {
    CU_ASSERT_TRUE(near(sodium_equivalent_mm(50, 2.0, 0.8), 181.45, 0.01));
    // All the Mg chelated by dNTP: only the monovalent is left.
    CU_ASSERT_DOUBLE_EQUAL(sodium_equivalent_mm(50, 0.5, 0.8), 50.0, 1e-12);
    CU_ASSERT_DOUBLE_EQUAL(sodium_equivalent_mm(50, 0.0, 0.0), 50.0, 1e-12);
}

// ------------------------------------------------------------------ Tm ---

// Against primer3's oligotm() over a spread of sequences and conditions,
// including the NEB control primers, a palindrome, a negative Tm and runs.
static void test_tm_matches_primer3_oligotm(void) {
    const char *seqs[] = {
        NEB_1_3KB_P1, NEB_1_3KB_P2, NEB_10KB_P1, NEB_10KB_P2, "GAATTC", "GCGCAAAAGCGC",
        "ACGTACGTACGTACGT", "AAAAAAAAAATTTTTTTTTT", "GGGGCCCCAAATTTGGGGCCCC",
        "TCTAGAGGATCCCCGGGTACC",
    };
    struct { double mono, mg, dntp, primer_um; } conditions[] = {
        {50, 2.0, 0.8, 0.5},   // Q5
        {50, 1.5, 0.8, 0.2},   // Taq
        {50, 0.0, 0.0, 0.25},  // no magnesium at all
        {100, 3.0, 0.8, 1.0},
    };
    for (const auto &c : conditions) {
        double sodium = sodium_equivalent_mm(c.mono, c.mg, c.dntp);
        for (const char *s : seqs) {
            double ours = melting_temp(s, sodium, c.primer_um);
            double theirs = primer3_tm(s, c.mono, c.mg, c.dntp, c.primer_um * 1000.0);
            CU_ASSERT_TRUE(near(ours, theirs, 1e-6));
        }
    }
}

// The same, pinned to numbers from primer3-py 2.3.1 so the check does not rest
// only on the copy of primer3 built alongside it.
static void test_tm_matches_primer3_py_for_each_polymerase(void) {
    struct { const char *pol; const char *seq; double tm; } cases[] = {
        {"Q5", NEB_1_3KB_P1, 68.9940},      {"Q5", NEB_1_3KB_P2, 67.9634},
        {"Q5", NEB_10KB_P1, 63.3089},       {"Q5", NEB_10KB_P2, 61.7664},
        {"Taq", NEB_1_3KB_P1, 67.0334},     {"Taq", NEB_1_3KB_P2, 66.0510},
        {"Taq", NEB_10KB_P1, 61.1843},      {"Taq", NEB_10KB_P2, 59.5977},
        {"Phusion", NEB_1_3KB_P1, 68.0176}, {"Phusion", NEB_1_3KB_P2, 66.9930},
        {"Phusion", NEB_10KB_P1, 62.3635},  {"Phusion", NEB_10KB_P2, 60.8546},
    };
    for (const auto &c : cases) {
        const Polymerase &p = builtin(c.pol);
        CU_ASSERT_TRUE(near(melting_temp(c.seq, reaction_sodium_mm(p), p.primer_um), c.tm,
                            0.00005 + 1e-9));
    }
}

// ---------------------------------------------------------- structures ---

// The primer3 manual's own worked examples, which it computes "with default
// Primer3 parameters" of that era: 50 mM monovalent, no Mg, 50 nM oligo.
static void test_thal_reproduces_the_primer3_manual(void) {
    const double hp1 = structure_tm("ACGCTGTGCTGCGA", "ACGCTGTGCTGCGA", thal_hairpin,
                                    50, 0, 0, 0.05);
    const double hp2 = structure_tm("CCGCAGTAAGCTGCGG", "CCGCAGTAAGCTGCGG", thal_hairpin,
                                    50, 0, 0, 0.05);
    CU_ASSERT_TRUE(near(hp1, 53.7263, 0.00005));
    CU_ASSERT_TRUE(near(hp2, 71.0918, 0.00005));
    // ATTAGATAGAGCATC against its complement.
    const std::string s = "ATTAGATAGAGCATC";
    CU_ASSERT_TRUE(near(structure_tm(s, reverse_complement(s), thal_any, 50, 0, 0, 0.05),
                        32.1493, 0.00005));
}

// Under Q5 conditions, against primer3-py 2.3.1's calc_hairpin and
// calc_homodimer. The first primer is NEB's own control primer: no hairpin.
static void test_thal_under_reaction_conditions(void) {
    const Polymerase &q5 = builtin("Q5");
    struct { const char *seq; double hairpin, self_any; } cases[] = {
        {NEB_10KB_P1, 0.0, 24.1606},
        {"GCGCAAAAGCGC", 51.4250, 25.2629},
        {"TCTAGAGGATCCCCGGGTACC", 41.4782, 28.9239},
        {NEB_10KB_P2, 0.0, 16.0978},
    };
    for (const auto &c : cases) {
        CU_ASSERT_TRUE(near(reaction_structure_tm(c.seq, c.seq, thal_hairpin, q5), c.hairpin,
                            0.00005));
        CU_ASSERT_TRUE(near(reaction_structure_tm(c.seq, c.seq, thal_any, q5), c.self_any,
                            0.00005));
    }
}

// A structure that never forms has a Tm at or below zero and reads as none;
// only one above primer3's 47 C limit is flagged.
static void test_structure_description(void) {
    bool bad = false;
    CU_ASSERT_EQUAL(describe_structure(0.0, &bad), std::string("none"));
    CU_ASSERT_EQUAL(describe_structure(-21.7, &bad), std::string("none"));
    CU_ASSERT_FALSE(bad);
    CU_ASSERT_EQUAL(describe_structure(47.0, &bad), std::string("47.0 C"));
    CU_ASSERT_FALSE(bad);
    CU_ASSERT_EQUAL(describe_structure(51.425, &bad), std::string("51.4 C (!)"));
    CU_ASSERT_TRUE(bad);
}

// ------------------------------------------------------------- annealing ---

// Rules, from each datasheet, worked by hand. Tms are given directly so these
// test the rule and not the Tm.
static void test_q5_anneals_three_above_the_lower_tm(void) {
    Annealing a = annealing_for(builtin("Q5"), 60.0, 62.0, 20, 20);
    CU_ASSERT_EQUAL(a.ta, 63);            // 60 + 3
    CU_ASSERT_FALSE(a.two_step);
    // Order does not matter, the lower Tm is used.
    CU_ASSERT_EQUAL(annealing_for(builtin("Q5"), 62.0, 60.0, 20, 20).ta, 63);
    // Rounds half up: 60.5 + 3 = 63.5 -> 64; 60.4 + 3 = 63.4 -> 63.
    CU_ASSERT_EQUAL(annealing_for(builtin("Q5"), 60.5, 70.0, 20, 20).ta, 64);
    CU_ASSERT_EQUAL(annealing_for(builtin("Q5"), 60.4, 70.0, 20, 20).ta, 63);
}

// Phusion: +3 only for primers longer than 20 nt, otherwise the lower Tm.
static void test_phusion_offset_depends_on_primer_length(void) {
    const Polymerase &p = builtin("Phusion");
    CU_ASSERT_EQUAL(annealing_for(p, 60.0, 65.0, 22, 25).ta, 63);
    CU_ASSERT_EQUAL(annealing_for(p, 60.0, 65.0, 20, 25).ta, 60);  // shorter one is 20
    CU_ASSERT_EQUAL(annealing_for(p, 60.0, 65.0, 21, 21).ta, 63);
}

// Taq: the gradient starts 5 C under the Tm, and 45 C is the bottom of range.
static void test_taq_anneals_five_below(void) {
    const Polymerase &p = builtin("Taq");
    Annealing a = annealing_for(p, 60.0, 62.0, 20, 20);
    CU_ASSERT_EQUAL(a.ta, 55);
    CU_ASSERT_FALSE(a.below_range);
    CU_ASSERT_FALSE(a.two_step_possible);
    CU_ASSERT_TRUE(annealing_for(p, 49.0, 55.0, 20, 20).below_range);   // 44 < 45
    CU_ASSERT_FALSE(annealing_for(p, 50.0, 55.0, 20, 20).below_range);  // 45 is in range
    // "Above 65 C" the datasheet allows a 2-step program but does not require
    // it. 65 itself is not above.
    Annealing at_65 = annealing_for(p, 70.0, 71.0, 20, 20);
    CU_ASSERT_EQUAL(at_65.ta, 65);
    CU_ASSERT_FALSE(at_65.two_step_possible);
    Annealing at_66 = annealing_for(p, 71.0, 72.0, 20, 20);
    CU_ASSERT_EQUAL(at_66.ta, 66);
    CU_ASSERT_TRUE(at_66.two_step_possible);
    CU_ASSERT_FALSE(at_66.two_step);
}

// Q5 at 70/71 C wants 73 C, past the 72 C extension: fold into 2-step at 72.
static void test_annealing_at_extension_goes_two_step(void) {
    Annealing a = annealing_for(builtin("Q5"), 70.0, 71.0, 25, 25);
    CU_ASSERT_TRUE(a.two_step);
    CU_ASSERT_EQUAL(a.ta, 72);
    // 69 + 3 = 72 exactly also goes 2-step.
    CU_ASSERT_TRUE(annealing_for(builtin("Q5"), 69.0, 71.0, 25, 25).two_step);
    CU_ASSERT_FALSE(annealing_for(builtin("Q5"), 68.0, 71.0, 25, 25).two_step);
}

// ------------------------------------------------------------- extension ---

// All by hand: ceil(bp * rate / 1000), then the datasheet floor.
static void test_extension_times(void) {
    const Polymerase &q5 = builtin("Q5"), &taq = builtin("Taq"), &phusion = builtin("Phusion");
    CU_ASSERT_EQUAL(extension_seconds(q5, 1347, false), 41);   // 40.41 s at 30 s/kb
    CU_ASSERT_EQUAL(extension_seconds(q5, 1347, true), 14);    // 13.47 s at 10 s/kb
    CU_ASSERT_EQUAL(extension_seconds(q5, 1000, false), 30);   // exact, no rounding up
    CU_ASSERT_EQUAL(extension_seconds(q5, 6000, false), 180);  // 6 kb is not "> 6 kb"
    CU_ASSERT_EQUAL(extension_seconds(q5, 6001, false), 301);  // 300.05 s at 50 s/kb
    CU_ASSERT_EQUAL(extension_seconds(q5, 10, false), 1);      // never zero
    CU_ASSERT_EQUAL(extension_seconds(taq, 1500, false), 90);  // 1 min/kb
    CU_ASSERT_EQUAL(extension_seconds(taq, 500, false), 45);   // 30 s, floored to 45
    CU_ASSERT_EQUAL(extension_seconds(taq, 500, true), 45);    // no faster rate for Taq
    CU_ASSERT_EQUAL(extension_seconds(phusion, 2000, false), 60);
    CU_ASSERT_EQUAL(extension_seconds(phusion, 2000, true), 30);  // 15 s/kb
}

// NEB's own cycling for its lambda controls (E0553 manual): 20 s extension for
// the 1.3 kb product and 2 min 30 s for the 10 kb, at Phusion's 15 s/kb.
// Rounding up to the second puts us 1 s over on both - never under.
static void test_extension_against_nebs_lambda_controls(void) {
    const Polymerase &p = builtin("Phusion");
    int short_amp = extension_seconds(p, 1347, true);
    int long_amp = extension_seconds(p, 10033, true);
    CU_ASSERT_TRUE(short_amp >= 20 && short_amp <= 21);
    CU_ASSERT_TRUE(long_amp >= 150 && long_amp <= 151);
}

// ------------------------------------------------------------- template ---

// NEB's 1.3 kb pair: 30,006..31,352 on the genome is 1,347 bp, at fixture
// positions 56..1,402.
static void test_neb_control_amplicon_found_at_nebs_coordinates(void) {
    std::vector<Product> products = find_products(lambda_fixture(), NEB_1_3KB_P1, NEB_1_3KB_P2);
    CU_ASSERT_EQUAL_FATAL(products.size(), 1u);
    CU_ASSERT_EQUAL(products[0].start + FIXTURE_OFFSET, 30006);
    CU_ASSERT_EQUAL(products[0].end + FIXTURE_OFFSET, 31352);
    CU_ASSERT_EQUAL(products[0].length(), 1347);
    CU_ASSERT_EQUAL(products[0].strand, '+');
}

// NEB's 10 kb pair: 30,011..40,043 is 10,033 bp.
static void test_neb_10kb_amplicon(void) {
    std::vector<Product> products = find_products(lambda_fixture(), NEB_10KB_P1, NEB_10KB_P2);
    CU_ASSERT_EQUAL_FATAL(products.size(), 1u);
    CU_ASSERT_EQUAL(products[0].start + FIXTURE_OFFSET, 30011);
    CU_ASSERT_EQUAL(products[0].end + FIXTURE_OFFSET, 40043);
    CU_ASSERT_EQUAL(products[0].length(), 10033);
}

// The reverse primer is 5'->3' as ordered, so it is found reverse complemented.
// Given as the top-strand sequence instead, there is no product, and the
// explanation says why. Pins the direction of the complement.
static void test_reverse_primer_orientation_is_pinned(void) {
    const std::string top_strand_p2 = reverse_complement(NEB_1_3KB_P2);
    CU_ASSERT_TRUE(find_products(lambda_fixture(), NEB_1_3KB_P1, top_strand_p2).empty());
    std::string why = explain_no_product(lambda_fixture(), NEB_1_3KB_P1, top_strand_p2);
    CU_ASSERT_TRUE(why.find("same strand") != std::string::npos);
}

// Which primer is called forward does not change the product, only which
// strand it reads along.
static void test_swapped_primers_make_the_same_product(void) {
    std::vector<Product> products = find_products(lambda_fixture(), NEB_1_3KB_P2, NEB_1_3KB_P1);
    CU_ASSERT_EQUAL_FATAL(products.size(), 1u);
    CU_ASSERT_EQUAL(products[0].length(), 1347);
    CU_ASSERT_EQUAL(products[0].strand, '-');
}

// A template given the other way round still yields the product.
static void test_reverse_complemented_template(void) {
    std::vector<std::pair<std::string, std::string> > flipped = {
        {"flipped", reverse_complement(lambda_fixture()[0].second)}};
    std::vector<Product> products = find_products(flipped, NEB_1_3KB_P1, NEB_1_3KB_P2);
    CU_ASSERT_EQUAL_FATAL(products.size(), 1u);
    CU_ASSERT_EQUAL(products[0].length(), 1347);
    CU_ASSERT_EQUAL(products[0].strand, '-');
}

static void test_primer_not_on_template(void) {
    CU_ASSERT_TRUE(find_products(lambda_fixture(), NEB_1_3KB_P1, "ACGTACGTACGTACGTACGT").empty());
    std::string why =
        explain_no_product(lambda_fixture(), NEB_1_3KB_P1, "ACGTACGTACGTACGTACGT");
    CU_ASSERT_TRUE(why.find("reverse primer does not match") != std::string::npos);
}

// Two sites for the pair: both reported, so the caller can refuse to guess.
static void test_repeated_site_gives_two_products(void) {
    const std::string unit = std::string("AAAA") + NEB_10KB_P1 + "CCCCCCCCCC" +
                             reverse_complement(NEB_10KB_P2) + "TTTT";
    std::vector<std::pair<std::string, std::string> > twice = {{"twice", unit + unit}};
    CU_ASSERT_TRUE(find_products(twice, NEB_10KB_P1, NEB_10KB_P2).size() >= 2);
}

// ----------------------------------------------------------- primer FASTA ---

// The fixture holds NEB's 1.3 kb pair, the reverse one lower case and wrapped
// over two lines, and must read back as the same two primers in order.
static void test_primer_fasta_reads_forward_then_reverse(void) {
    std::pair<std::string, std::string> pair =
        primers_from_records(parse_fasta(PCR_TEST_PRIMERS), "t.fna");
    CU_ASSERT_EQUAL(pair.first, std::string(NEB_1_3KB_P1));
    CU_ASSERT_EQUAL(pair.second, std::string(NEB_1_3KB_P2));
    // And they make NEB's product.
    std::vector<Product> products = find_products(lambda_fixture(), pair.first, pair.second);
    CU_ASSERT_EQUAL_FATAL(products.size(), 1u);
    CU_ASSERT_EQUAL(products[0].length(), 1347);
}

// Order is all that says which primer is which, so anything but two is refused.
static void test_primer_fasta_needs_exactly_two(void) {
    std::vector<std::pair<std::string, std::string> > one = {{"a", NEB_1_3KB_P1}};
    std::vector<std::pair<std::string, std::string> > three = {
        {"a", NEB_1_3KB_P1}, {"b", NEB_1_3KB_P2}, {"c", NEB_10KB_P1}};
    CU_ASSERT_TRUE(throws([&] { primers_from_records(one, "t.fna"); }));
    CU_ASSERT_TRUE(throws([&] { primers_from_records(three, "t.fna"); }));
    // Same validation as -F/-R.
    std::vector<std::pair<std::string, std::string> > degenerate = {
        {"a", NEB_1_3KB_P1}, {"b", "GATGRCGC"}};
    CU_ASSERT_TRUE(throws([&] { primers_from_records(degenerate, "t.fna"); }));
}

// ------------------------------------------------------------ polymerases ---

static void test_builtin_table(void) {
    const std::vector<Polymerase> &table = builtin_polymerases();
    CU_ASSERT_EQUAL(table.size(), 3u);
    CU_ASSERT_EQUAL(builtin("q5").name, std::string("Q5"));      // case insensitive
    CU_ASSERT_EQUAL(builtin("PHUSION").catalog, std::string("NEB M0530"));
    CU_ASSERT_TRUE(throws([] { builtin("Pfu"); }));
}

// Spot checks against the datasheets, so an edit to the CSV that breaks a
// sourced number is caught.
static void test_builtin_values_match_the_datasheets(void) {
    const Polymerase &q5 = builtin("Q5");
    CU_ASSERT_EQUAL(q5.denature_c, 98);
    CU_ASSERT_EQUAL(q5.extend_c, 72);
    CU_ASSERT_EQUAL(q5.final_extend_s, 120);
    CU_ASSERT_DOUBLE_EQUAL(q5.mg_mm, 2.0, 1e-12);
    const Polymerase &taq = builtin("Taq");
    CU_ASSERT_EQUAL(taq.denature_c, 95);
    CU_ASSERT_EQUAL(taq.extend_c, 68);
    CU_ASSERT_EQUAL(taq.extend_s_per_kb, 60);
    CU_ASSERT_DOUBLE_EQUAL(taq.mg_mm, 1.5, 1e-12);
    CU_ASSERT_DOUBLE_EQUAL(taq.enzyme_units_per_50ul, 1.25, 1e-12);
    const Polymerase &phusion = builtin("Phusion");
    CU_ASSERT_EQUAL(phusion.anneal_offset_min_primer_nt, 21);
    CU_ASSERT_EQUAL(phusion.extend_s_per_kb_simple, 15);
}

static void test_csv_missing_column_rejected(void) {
    CU_ASSERT_TRUE(throws([] { load_polymerase_csv_string("polymerase,buffer\nX,Y\n", "t.csv"); }));
}

static void test_csv_bad_number_rejected(void) {
    std::string text = EMBEDDED_POLYMERASE_CSV;
    size_t at = text.find(",98,30,98,10,");  // Q5's denaturation columns
    CU_ASSERT_FATAL(at != std::string::npos);
    text.replace(at, 4, ",9x,");
    CU_ASSERT_TRUE(throws([&] { load_polymerase_csv_string(text, "t.csv"); }));
}

// ------------------------------------------------------------ primer QC ---

static void test_primer_validation(void) {
    CU_ASSERT_EQUAL(validate_primer("gtc acc\tag", "F"), std::string("GTCACCAG"));
    CU_ASSERT_TRUE(throws([] { validate_primer("GTCRACC", "F"); }));   // degenerate
    CU_ASSERT_TRUE(throws([] { validate_primer("G", "F"); }));
}

static void test_runs_and_repeats(void) {
    CU_ASSERT_EQUAL(longest_run("ACGTTTTTAC"), 5);
    CU_ASSERT_EQUAL(longest_run("ACGT"), 1);
    CU_ASSERT_EQUAL(most_dinucleotide_repeats("GATATATATATC"), 5);   // AT x5 (TA x5 too)
    CU_ASSERT_EQUAL(most_dinucleotide_repeats("AAAAAA"), 0);         // a run, not a repeat
    CU_ASSERT_EQUAL(gc_in_last("AAAAAGGCGC", 5), 5);
    CU_ASSERT_EQUAL(gc_in_last("GGGGGAAAAT", 5), 0);
    CU_ASSERT_TRUE(near(gc_percent("GGCCAATT"), 50.0, 1e-12));
}

// ------------------------------------------------------------ master mix ---

// By hand, Q5 at 50 uL, 8 reactions, 1 uL template: reagents per reaction are
// 10 + 1 + 2.5 + 2.5 + 0.5 = 16.5 uL, so water is 50 - 1 - 16.5 = 32.5 uL.
// With 10% extra that is 8.8 reactions: 88, 8.8, 22, 22, 4.4 and 286 uL.
static void test_master_mix_by_hand(void) {
    MasterMix mix = master_mix(builtin("Q5"), 50.0, 8, 1.0);
    CU_ASSERT_TRUE(near(mix.reactions, 8.8, 1e-9));
    CU_ASSERT_TRUE(near(mix.aliquot_ul, 49.0, 1e-9));
    const double expected_mix[] = {88.0, 8.8, 22.0, 22.0, 4.4, 286.0};
    CU_ASSERT_EQUAL_FATAL(mix.per_reaction.size(), 6u);
    for (size_t k = 0; k < 6; k++) {
        CU_ASSERT_TRUE(near(mix.per_reaction[k].ul * mix.reactions, expected_mix[k], 1e-9));
    }
    CU_ASSERT_EQUAL(mix.per_reaction.back().name, std::string("Nuclease-free water"));
    // Per tube, mix plus template is the reaction volume.
    double per_tube = 0;
    for (const Component &c : mix.per_reaction) { per_tube += c.ul; }
    CU_ASSERT_TRUE(near(per_tube + 1.0, 50.0, 1e-9));
}

// Taq at 25 uL, 3 reactions, 2 uL template: 2.5 + 0.5 + 0.5 + 0.5 + 0.125 =
// 4.125 uL of reagents, water 25 - 2 - 4.125 = 18.875, times 3.3 reactions.
static void test_master_mix_scales_with_volume(void) {
    MasterMix mix = master_mix(builtin("Taq"), 25.0, 3, 2.0);
    CU_ASSERT_TRUE(near(mix.reactions, 3.3, 1e-9));
    CU_ASSERT_TRUE(near(mix.per_reaction.back().ul, 18.875, 1e-9));
    CU_ASSERT_TRUE(near(mix.per_reaction[4].ul * mix.reactions, 0.4125, 1e-9));  // enzyme
    CU_ASSERT_TRUE(near(mix.aliquot_ul, 23.0, 1e-9));
}

static void test_master_mix_rejects_what_does_not_fit(void) {
    // 33.5 uL of template leaves exactly no water; 34 does not fit.
    CU_ASSERT_FALSE(throws([] { master_mix(builtin("Q5"), 50.0, 8, 33.5); }));
    CU_ASSERT_TRUE(throws([] { master_mix(builtin("Q5"), 50.0, 8, 34.0); }));
    CU_ASSERT_TRUE(throws([] { master_mix(builtin("Q5"), 50.0, 0, 1.0); }));
    CU_ASSERT_TRUE(throws([] { master_mix(builtin("Q5"), 50.0, 8, -1.0); }));
}

// ---------------------------------------------------------------- output ---

static void test_formatting(void) {
    CU_ASSERT_EQUAL(volume(2.5), std::string("2.5"));
    CU_ASSERT_EQUAL(volume(0.25), std::string("0.25"));
    CU_ASSERT_EQUAL(volume(10.0), std::string("10"));
    CU_ASSERT_EQUAL(duration(41), std::string("41 s"));
    CU_ASSERT_EQUAL(duration(120), std::string("2 min"));
    CU_ASSERT_EQUAL(duration(150), std::string("2 min 30 s"));
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
        {"nearest neighbour parameters", {
            {"Table 2 agrees with Table 1 dG", test_nn_table_agrees_with_the_papers_dg_column},
            {"unlisted stacks are the other strand", test_unlisted_stacks_are_the_other_strand},
            {"self complementary is pinned", test_self_complementary_is_pinned},
            {"sodium equivalents", test_sodium_equivalents},
        }},
        {"Tm", {
            {"matches primer3 oligotm", test_tm_matches_primer3_oligotm},
            {"matches primer3-py per polymerase", test_tm_matches_primer3_py_for_each_polymerase},
        }},
        {"structures (primer3 thal)", {
            {"reproduces the primer3 manual", test_thal_reproduces_the_primer3_manual},
            {"under reaction conditions", test_thal_under_reaction_conditions},
            {"description and 47 C limit", test_structure_description},
        }},
        {"annealing", {
            {"Q5 is lower Tm + 3", test_q5_anneals_three_above_the_lower_tm},
            {"Phusion offset depends on length", test_phusion_offset_depends_on_primer_length},
            {"Taq is lower Tm - 5", test_taq_anneals_five_below},
            {"at extension temp goes 2-step", test_annealing_at_extension_goes_two_step},
        }},
        {"extension", {
            {"hand worked times", test_extension_times},
            {"against NEB's lambda controls", test_extension_against_nebs_lambda_controls},
        }},
        {"template", {
            {"NEB 1.3 kb control at NEB's coordinates", test_neb_control_amplicon_found_at_nebs_coordinates},
            {"NEB 10 kb control", test_neb_10kb_amplicon},
            {"reverse primer orientation is pinned", test_reverse_primer_orientation_is_pinned},
            {"swapped primers, same product", test_swapped_primers_make_the_same_product},
            {"reverse complemented template", test_reverse_complemented_template},
            {"primer not on template", test_primer_not_on_template},
            {"repeated site gives two products", test_repeated_site_gives_two_products},
        }},
        {"primer FASTA", {
            {"reads forward then reverse", test_primer_fasta_reads_forward_then_reverse},
            {"needs exactly two", test_primer_fasta_needs_exactly_two},
        }},
        {"polymerase table", {
            {"built in table", test_builtin_table},
            {"values match the datasheets", test_builtin_values_match_the_datasheets},
            {"missing column rejected", test_csv_missing_column_rejected},
            {"bad number rejected", test_csv_bad_number_rejected},
        }},
        {"primer QC", {
            {"primer validation", test_primer_validation},
            {"runs, repeats, GC clamp", test_runs_and_repeats},
        }},
        {"master mix", {
            {"by hand", test_master_mix_by_hand},
            {"scales with volume", test_master_mix_scales_with_volume},
            {"rejects what does not fit", test_master_mix_rejects_what_does_not_fit},
        }},
        {"output", {
            {"formatting", test_formatting},
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
