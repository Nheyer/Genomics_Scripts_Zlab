// SPDX-License-Identifier: GPL-3.0-only
//
// Unit tests for Digest_Protocol.
//
// Same arrangement as the other suites: the program is one translation unit
// that owns its own main(), so we rename that main out of the way and pull the
// whole file in.
//
// None of the expected numbers here come from running this code. The amounts are
// pinned to NEB's own table for 10, 25 and 50 uL reactions, typed in here
// separately from digest_data.hpp so a slip in either fails; everything else is
// worked out by hand in the comment above the test.
#define main disabled_main
#include "Digest_Protocol.cpp"
#undef main

#include <CUnit/Basic.h>

#include <initializer_list>
#include <string>
#include <vector>

// ---------------------------------------------------------------- helpers --

template <typename Fn>
static bool throws(Fn fn) {
    try {
        fn();
    } catch (const std::exception &) {
        return true;
    }
    return false;
}

template <typename Fn>
static std::string message_of(Fn fn) {
    try {
        fn();
    } catch (const std::exception &e) {
        return e.what();
    }
    return "";
}

static bool contains(const std::string &text, const std::string &part) {
    return text.find(part) != std::string::npos;
}

// An enzyme with nothing known about it but its name, or with the given stock.
static Enzyme plain(const char *name, double stock = 0) {
    Enzyme e;
    e.name = name;
    e.stock_u_per_ul = stock;
    return e;
}

// An enzyme with its buffer and temperature known.
static Enzyme known(const char *name, const char *buffer, int temperature) {
    Enzyme e = plain(name);
    e.buffer = buffer;
    e.incubation_c = temperature;
    return e;
}

static const double EPS = 1e-9;

static const char *ROW_HEADER =
    "enzyme,buffer,incubation_c,inactivate_c,inactivate_min,time_saver,stock_u_per_ul,notes,source\n";

static std::vector<Conditions> table_of(const std::string &rows) {
    return load_conditions_csv_string(std::string(ROW_HEADER) + rows, "t.csv");
}

// ---------------------------------------------------------- NEB's own table ---

// From NEB's Restriction Digest protocol, the table for smaller volumes:
//   10 uL  1 unit    0.1 ug   1 uL buffer
//   25 uL  5 units   0.5 ug   2.5 uL buffer
//   50 uL  10 units  1 ug     5 uL buffer
// with an enzyme stock of 10 U/uL, "10 units ... generally 1 uL".
static void test_neb_small_volume_table(void) {
    struct row { double volume, units, ug, buffer_ul; };
    const row rows[] = {{10, 1, 0.1, 1}, {25, 5, 0.5, 2.5}, {50, 10, 1, 5}};
    for (const row &r : rows) {
        const double ug = max_dna_ug(r.volume);
        CU_ASSERT_DOUBLE_EQUAL(ug, r.ug, EPS);
        Reaction rxn = plan_reaction({plain("EcoRI")}, "rCutSmart", r.volume, ug, 1.0, UNITS_PER_UG);
        CU_ASSERT_DOUBLE_EQUAL(rxn.units, r.units, EPS);
        CU_ASSERT_DOUBLE_EQUAL(rxn.parts[0].ul, r.buffer_ul, EPS);
        CU_ASSERT_DOUBLE_EQUAL(rxn.parts[1].ul, r.units / 10.0, EPS);
    }
}

// Above the last row it is 1 ug per 50 uL, and between rows it is the row below.
static void test_dna_limit_between_and_beyond_the_rows(void) {
    CU_ASSERT_DOUBLE_EQUAL(max_dna_ug(100), 2.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(max_dna_ug(60), 1.2, EPS);
    CU_ASSERT_DOUBLE_EQUAL(max_dna_ug(30), 0.5, EPS);    // the 25 uL row
    CU_ASSERT_DOUBLE_EQUAL(max_dna_ug(17), 0.1, EPS);    // the 10 uL row
    CU_ASSERT_DOUBLE_EQUAL(max_dna_ug(9), 0.0, EPS);     // under the table: no figure
}

// ---------------------------------------------------------------- setup ---

// 50 uL, 1 ug in 1 uL, one enzyme at 10 U/uL: buffer 5, enzyme 1, DNA 1, so the
// water is 50 - 5 - 1 - 1 = 43.
static void test_typical_setup_by_hand(void) {
    Reaction r = plan_reaction({plain("EcoRI")}, "rCutSmart", 50, 1.0, 1.0, 10.0);
    CU_ASSERT_EQUAL(r.parts.size(), 2u);
    CU_ASSERT_EQUAL(r.parts[0].name, std::string("10X rCutSmart"));
    CU_ASSERT_DOUBLE_EQUAL(r.parts[0].ul, 5.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(r.parts[1].ul, 1.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(r.enzyme_ul, 1.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(r.water_ul, 43.0, EPS);
}

// Genomic DNA takes 20 units per ug: 20 units at 10 U/uL is 2 uL, water 42.
static void test_genomic_takes_twenty_units_per_ug(void) {
    Reaction r = plan_reaction({plain("EcoRI")}, "rCutSmart", 50, 1.0, 1.0, GENOMIC_UNITS_PER_UG);
    CU_ASSERT_DOUBLE_EQUAL(r.units, 20.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(r.enzyme_ul, 2.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(r.water_ul, 42.0, EPS);
}

// A double digest gives each enzyme its whole dose: two enzymes, 10 units each,
// 1 uL each, water 50 - 5 - 1 - 1 - 1 = 42.
static void test_each_enzyme_gets_the_full_dose(void) {
    Reaction r = plan_reaction({plain("EcoRI"), plain("BamHI")}, "rCutSmart", 50, 1.0, 1.0, 10.0);
    CU_ASSERT_EQUAL(r.parts.size(), 3u);
    CU_ASSERT_DOUBLE_EQUAL(r.parts[1].ul, 1.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(r.parts[2].ul, 1.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(r.enzyme_ul, 2.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(r.water_ul, 42.0, EPS);
}

// A 20 U/uL stock gives 10 units in half a microlitre.
static void test_stock_sets_the_enzyme_volume(void) {
    Reaction r = plan_reaction({plain("EcoRI", 20.0)}, "rCutSmart", 50, 1.0, 1.0, 10.0);
    CU_ASSERT_DOUBLE_EQUAL(r.parts[1].ul, 0.5, EPS);
    CU_ASSERT_TRUE(contains(r.parts[1].name, "20 U/uL"));
    CU_ASSERT_EQUAL(r.parts[1].final, std::string("10 U"));
}

// An enzyme that has no stock on file gets the protocol's 10 U/uL, and says so.
static void test_default_stock_is_the_protocols_ratio(void) {
    CU_ASSERT_DOUBLE_EQUAL(effective_stock(plain("X")), 10.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(effective_stock(plain("X", 20.0)), 20.0, EPS);
}

// ---------------------------------------------------- what NEB says not to do ---

// 5 ug in 50 uL is 50 units; at 10 U/uL that is 5 uL, exactly 10% of 50 uL and
// so allowed (water 50 - 5 - 5 - 1 = 39). A second enzyme doubles it to 20%.
static void test_enzyme_volume_is_capped_at_ten_percent(void) {
    Reaction ok = plan_reaction({plain("A")}, "b", 50, 5.0, 1.0, 10.0);
    CU_ASSERT_DOUBLE_EQUAL(ok.enzyme_ul, 5.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(ok.water_ul, 39.0, EPS);
    CU_ASSERT_TRUE(throws([] {
        plan_reaction({plain("A"), plain("B")}, "b", 50, 5.0, 1.0, 10.0);
    }));
    CU_ASSERT_TRUE(contains(message_of([] {
        plan_reaction({plain("A"), plain("B")}, "b", 50, 5.0, 1.0, 10.0);
    }), "10%"));
}

// The reagents take 5 + 1 = 6 uL of 50, so 45 uL of DNA is 1 uL too many.
static void test_dna_that_does_not_fit_is_refused(void) {
    CU_ASSERT_TRUE(throws([] { plan_reaction({plain("A")}, "b", 50, 1.0, 45.0, 10.0); }));
    Reaction ok = plan_reaction({plain("A")}, "b", 50, 1.0, 44.0, 10.0);
    CU_ASSERT_DOUBLE_EQUAL(ok.water_ul, 0.0, EPS);
}

static void test_nonsense_inputs_are_refused(void) {
    CU_ASSERT_TRUE(throws([] { plan_reaction({plain("A")}, "b", 0, 1.0, 1.0, 10.0); }));
    CU_ASSERT_TRUE(throws([] { plan_reaction({plain("A")}, "b", 50, 0.0, 1.0, 10.0); }));
    CU_ASSERT_TRUE(throws([] { plan_reaction({plain("A")}, "b", 50, 1.0, -1.0, 10.0); }));
    CU_ASSERT_TRUE(throws([] { plan_reaction({plain("A")}, "b", 50, 1.0, 1.0, 0.0); }));
    CU_ASSERT_TRUE(throws([] { plan_reaction({}, "b", 50, 1.0, 1.0, 10.0); }));
}

// ------------------------------------------------------------- master mix ---

// 4 reactions + 10% = 4.4. Per reaction: buffer 5, enzyme 1, water 43; the DNA
// is not in the mix, so the aliquot is 50 - 1 = 49. Mix: 22, 4.4 and
// 43 x 4.4 = 189.2.
static void test_master_mix_by_hand(void) {
    Reaction r = plan_reaction({plain("EcoRI")}, "rCutSmart", 50, 1.0, 1.0, 10.0);
    MasterMix mix = master_mix(r, 4);
    CU_ASSERT_DOUBLE_EQUAL(mix.reactions, 4.4, EPS);
    CU_ASSERT_DOUBLE_EQUAL(mix.aliquot_ul, 49.0, EPS);
    CU_ASSERT_EQUAL(mix.per_reaction.size(), 3u);
    CU_ASSERT_DOUBLE_EQUAL(mix.per_reaction[0].ul * mix.reactions, 22.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(mix.per_reaction[1].ul * mix.reactions, 4.4, EPS);
    CU_ASSERT_EQUAL(mix.per_reaction[2].name, std::string("Nuclease-free water"));
    CU_ASSERT_DOUBLE_EQUAL(mix.per_reaction[2].ul * mix.reactions, 189.2, EPS);
}

static void test_master_mix_needs_a_reaction(void) {
    Reaction r = plan_reaction({plain("EcoRI")}, "rCutSmart", 50, 1.0, 1.0, 10.0);
    CU_ASSERT_TRUE(throws([&] { master_mix(r, 0); }));
}

// "10 ul per 50 ul rxn", so a fifth of the volume.
static void test_stop_solution_is_a_fifth_of_the_volume(void) {
    CU_ASSERT_DOUBLE_EQUAL(stop_solution_ul(50), 10.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(stop_solution_ul(25), 5.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(stop_solution_ul(10), 2.0, EPS);
    CU_ASSERT_DOUBLE_EQUAL(stop_solution_ul(100), 20.0, EPS);
}

// ------------------------------------------------------- the conditions CSV ---

static void test_conditions_row_is_read(void) {
    std::vector<Conditions> t = table_of("EcoRI-HF,rCutSmart,37,65,20,1,20,a note,NEB chart\n");
    CU_ASSERT_EQUAL(t.size(), 1u);
    CU_ASSERT_EQUAL(t[0].name, std::string("EcoRI-HF"));
    CU_ASSERT_EQUAL(t[0].buffer, std::string("rCutSmart"));
    CU_ASSERT_EQUAL(t[0].incubation_c, 37);
    CU_ASSERT_EQUAL(t[0].inactivate_c, 65);
    CU_ASSERT_EQUAL(t[0].inactivate_min, 20);
    CU_ASSERT_EQUAL(t[0].time_saver, 1);
    CU_ASSERT_DOUBLE_EQUAL(t[0].stock_u_per_ul, 20.0, EPS);
    CU_ASSERT_EQUAL(t[0].notes, std::string("a note"));
    CU_ASSERT_EQUAL(t[0].source, std::string("NEB chart"));
}

// Nothing on file is a fact about the table, not a broken file.
static void test_header_only_table_is_empty_not_an_error(void) {
    CU_ASSERT_EQUAL(table_of("").size(), 0u);
}

static void test_conditions_columns_can_be_in_any_order(void) {
    std::vector<Conditions> t = load_conditions_csv_string(
        "source,notes,stock_u_per_ul,time_saver,inactivate_min,inactivate_c,incubation_c,buffer,"
        "enzyme\nS,N,10,0,0,0,25,NEBuffer 4,SmaI\n", "t.csv");
    CU_ASSERT_EQUAL(t[0].name, std::string("SmaI"));
    CU_ASSERT_EQUAL(t[0].buffer, std::string("NEBuffer 4"));
    CU_ASSERT_EQUAL(t[0].incubation_c, 25);
}

static void test_conditions_missing_column_rejected(void) {
    CU_ASSERT_TRUE(throws([] {
        load_conditions_csv_string("enzyme,buffer\nEcoRI,rCutSmart\n", "t.csv");
    }));
    CU_ASSERT_TRUE(contains(message_of([] {
        load_conditions_csv_string("enzyme,buffer\nEcoRI,rCutSmart\n", "t.csv");
    }), "incubation_c"));
}

static void test_conditions_bad_number_rejected(void) {
    CU_ASSERT_TRUE(throws([] { table_of("EcoRI,rCutSmart,warm,65,20,0,10,,src\n"); }));
    CU_ASSERT_TRUE(throws([] { table_of("EcoRI,rCutSmart,37.5,65,20,0,10,,src\n"); }));
    CU_ASSERT_TRUE(throws([] { table_of("EcoRI,rCutSmart,37,65,20,0,lots,,src\n"); }));
}

// The whole point of the table: no row without a source.
static void test_conditions_row_without_a_source_rejected(void) {
    CU_ASSERT_TRUE(throws([] { table_of("EcoRI,rCutSmart,37,65,20,0,10,,\n"); }));
    CU_ASSERT_TRUE(contains(message_of([] { table_of("EcoRI,rCutSmart,37,65,20,0,10,,\n"); }),
                            "source"));
}

// A temperature with no time, or a time with no temperature, is half a fact.
static void test_conditions_inactivation_is_both_or_neither(void) {
    CU_ASSERT_TRUE(throws([] { table_of("A,b,37,65,0,0,10,,s\n"); }));
    CU_ASSERT_TRUE(throws([] { table_of("A,b,37,0,20,0,10,,s\n"); }));
    CU_ASSERT_TRUE(throws([] { table_of("A,b,37,-65,20,0,10,,s\n"); }));
    CU_ASSERT_EQUAL(table_of("A,b,37,0,0,0,10,,s\n")[0].inactivate_c, 0);   // "none listed"
    CU_ASSERT_EQUAL(table_of("A,b,37,80,20,0,10,,s\n")[0].inactivate_c, 80);
}

static void test_conditions_other_fields_are_checked(void) {
    CU_ASSERT_TRUE(throws([] { table_of("A,b,0,65,20,0,10,,s\n"); }));     // no temperature
    CU_ASSERT_TRUE(throws([] { table_of("A,b,37,65,20,2,10,,s\n"); }));    // time_saver is -1, 0 or 1
    CU_ASSERT_TRUE(throws([] { table_of("A,b,37,65,20,0,-1,,s\n"); }));    // stock cannot be negative
    CU_ASSERT_TRUE(throws([] { table_of(",b,37,65,20,0,10,,s\n"); }));     // no name
    CU_ASSERT_TRUE(throws([] { table_of("A,,37,65,20,0,10,,s\n"); }));     // no buffer
}

// A chart that does not carry stock or Time-Saver status is not a broken row:
// 0 stock and -1 time_saver are "not known", same as an enzyme with no row at all.
static void test_conditions_stock_and_time_saver_can_be_unknown(void) {
    Conditions c = table_of("A,b,37,65,20,-1,0,,s\n")[0];
    CU_ASSERT_DOUBLE_EQUAL(c.stock_u_per_ul, 0.0, EPS);
    CU_ASSERT_EQUAL(c.time_saver, -1);
}

// NEB's NEBuffer Performance Chart gives a heat-inactivation temperature with no
// time at all; that is a fact, not a broken row, and different from both "not
// known" (-1/-1) and "NEB lists none" (0/0).
static void test_conditions_inactivation_temperature_with_no_time(void) {
    Conditions c = table_of("A,b,37,65,-1,-1,0,,s\n")[0];
    CU_ASSERT_EQUAL(c.inactivate_c, 65);
    CU_ASSERT_EQUAL(c.inactivate_min, -1);
    Conditions unknown = table_of("A,b,37,-1,-1,-1,0,,s\n")[0];
    CU_ASSERT_EQUAL(unknown.inactivate_c, -1);
    CU_ASSERT_EQUAL(unknown.inactivate_min, -1);
    // still refused: a temperature with an explicit 0 minutes is a half fact,
    // not "not known" - that is what -1 is for.
    CU_ASSERT_TRUE(throws([] { table_of("A,b,37,65,0,-1,0,,s\n"); }));
    CU_ASSERT_TRUE(throws([] { table_of("A,b,37,-2,-1,-1,0,,s\n"); }));   // below -1
    CU_ASSERT_TRUE(throws([] { table_of("A,b,37,-1,0,-1,0,,s\n"); }));    // -1 paired with 0
}

static void test_conditions_duplicate_enzyme_rejected(void) {
    CU_ASSERT_TRUE(throws([] {
        table_of("EcoRI,b,37,65,20,0,10,,s\necori,b,37,65,20,0,10,,s\n");
    }));
}

// Whatever rows are in the shipped table, they have to have got through the loader.
static void test_builtin_table_loads(void) {
    CU_ASSERT_TRUE(!throws([] { builtin_conditions(); }));
    for (const Conditions &c : builtin_conditions()) {
        CU_ASSERT_TRUE(!c.source.empty());
        CU_ASSERT_TRUE(!c.buffer.empty());
    }
}

// ----------------------------------------------------------- enzyme specs ---

static void test_plain_name_resolves_case_insensitively(void) {
    Enzyme e = parse_enzyme_spec("ecori", {});
    CU_ASSERT_EQUAL(e.name, std::string("EcoRI"));      // reported as the table spells it
    CU_ASSERT_EQUAL(e.site, std::string("G^AATT_C"));   // as in Data/restriction_enzymes.csv
    CU_ASSERT_TRUE(e.buffer.empty());
    CU_ASSERT_EQUAL(e.incubation_c, 0);
    CU_ASSERT_EQUAL(e.inactivate_c, -1);                // unknown, not "none"
    CU_ASSERT_TRUE(!e.from_table && !e.given);
}

// The nicking enzymes are real NEB enzymes with reactions of their own.
static void test_nicking_enzyme_names_resolve(void) {
    CU_ASSERT_EQUAL(parse_enzyme_spec("nb.bsmi", {}).name, std::string("Nb.BsmI"));
}

static void test_unknown_enzyme_is_refused_with_the_way_out(void) {
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRZ", {}); }));
    CU_ASSERT_TRUE(contains(message_of([] { parse_enzyme_spec("EcoRZ", {}); }), ":BUFFER:TEMP_C"));
}

// An enzyme neither table lists can still be used, if the person supplies what
// the tables would have.
static void test_conditions_on_the_command_line(void) {
    Enzyme e = parse_enzyme_spec("MyEnz:rCutSmart:37", {});
    CU_ASSERT_EQUAL(e.name, std::string("MyEnz"));
    CU_ASSERT_EQUAL(e.buffer, std::string("rCutSmart"));
    CU_ASSERT_EQUAL(e.incubation_c, 37);
    CU_ASSERT_TRUE(e.given && !e.from_table);
    CU_ASSERT_TRUE(e.site.empty());
    CU_ASSERT_EQUAL(e.inactivate_c, -1);
}

static void test_inactivation_on_the_command_line(void) {
    Enzyme heat = parse_enzyme_spec("EcoRI:rCutSmart:37:65/20", {});
    CU_ASSERT_EQUAL(heat.inactivate_c, 65);
    CU_ASSERT_EQUAL(heat.inactivate_min, 20);
    Enzyme none = parse_enzyme_spec("EcoRI:rCutSmart:37:no", {});
    CU_ASSERT_EQUAL(none.inactivate_c, 0);
    CU_ASSERT_EQUAL(parse_enzyme_spec("EcoRI:rCutSmart:37:NO", {}).inactivate_c, 0);
}

// A source giving only a temperature (NEB's NEBuffer Performance Chart, for
// instance) can be typed as the temperature alone, with no fabricated time.
static void test_inactivation_temperature_only_on_the_command_line(void) {
    Enzyme e = parse_enzyme_spec("EcoRI:rCutSmart:37:65", {});
    CU_ASSERT_EQUAL(e.inactivate_c, 65);
    CU_ASSERT_EQUAL(e.inactivate_min, -1);
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRI:rCutSmart:37:0", {}); }));    // not positive
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRI:rCutSmart:37:warm", {}); }));
}

static void test_malformed_specs_are_refused(void) {
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRI:rCutSmart", {}); }));          // no temperature
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRI:a:37:65/20:extra", {}); }));   // too many
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRI:a:warm", {}); }));
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRI:a:0", {}); }));
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRI::37", {}); }));                // no buffer
    // "EcoRI:a:37:65" (temperature, no minutes) is not malformed: see
    // test_inactivation_temperature_only_on_the_command_line.
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRI:a:37:65/x", {}); }));
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("EcoRI:a:37:0/20", {}); }));
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec("", {}); }));
    CU_ASSERT_TRUE(throws([] { parse_enzyme_spec(":a:37", {}); }));
}

// A conditions row fills everything in; typing its own values back changes
// nothing, while a real difference is flagged so its source is not misattributed.
static void test_command_line_overrides_the_table_row(void) {
    std::vector<Conditions> t = table_of("EcoRI,rCutSmart,37,65,20,1,20,n,NEB chart\n");
    Enzyme row = parse_enzyme_spec("ecori", t);
    CU_ASSERT_TRUE(row.from_table && !row.given && !row.overrides_table);
    CU_ASSERT_EQUAL(row.buffer, std::string("rCutSmart"));
    CU_ASSERT_DOUBLE_EQUAL(row.stock_u_per_ul, 20.0, EPS);
    CU_ASSERT_EQUAL(row.source, std::string("NEB chart"));

    Enzyme same = parse_enzyme_spec("EcoRI:rCutSmart:37", t);
    CU_ASSERT_TRUE(same.from_table && same.given && !same.overrides_table);

    Enzyme other = parse_enzyme_spec("EcoRI:NEBuffer 2:37", t);
    CU_ASSERT_TRUE(other.overrides_table);
    CU_ASSERT_EQUAL(other.buffer, std::string("NEBuffer 2"));
    CU_ASSERT_EQUAL(other.inactivate_c, 65);              // not given, so still the row's

    CU_ASSERT_TRUE(parse_enzyme_spec("EcoRI:rCutSmart:37:80/20", t).overrides_table);
}

// The -HF versions are not in the enzyme table, only in a conditions table.
static void test_enzyme_only_in_the_conditions_table_resolves(void) {
    std::vector<Conditions> t = table_of("EcoRI-HF,rCutSmart,37,65,20,1,20,n,s\n");
    Enzyme e = parse_enzyme_spec("ecori-hf", t);
    CU_ASSERT_EQUAL(e.name, std::string("EcoRI-HF"));
    CU_ASSERT_TRUE(e.from_table);
    CU_ASSERT_TRUE(e.site.empty());
}

// --------------------------------------------- what the enzymes agree about ---

static void test_buffer_agreement(void) {
    Common same = common_buffer({known("A", "rCutSmart", 37), known("B", "rcutsmart", 37)});
    CU_ASSERT_TRUE(same.how == Agreement::agree);
    CU_ASSERT_EQUAL(same.text, std::string("rCutSmart"));
    CU_ASSERT_TRUE(common_buffer({known("A", "rCutSmart", 37), known("B", "NEBuffer 2", 37)}).how ==
                   Agreement::disagree);
    CU_ASSERT_TRUE(common_buffer({known("A", "rCutSmart", 37), plain("B")}).how == Agreement::unknown);
    CU_ASSERT_TRUE(common_buffer({plain("A")}).how == Agreement::unknown);
}

static void test_temperature_agreement(void) {
    Common same = common_temperature({known("A", "b", 37), known("B", "b", 37)});
    CU_ASSERT_TRUE(same.how == Agreement::agree);
    CU_ASSERT_EQUAL(same.text, std::string("37"));
    CU_ASSERT_TRUE(common_temperature({known("A", "b", 37), known("B", "b", 65)}).how ==
                   Agreement::disagree);
    CU_ASSERT_TRUE(common_temperature({known("A", "b", 37), plain("B")}).how == Agreement::unknown);
}

// One known disagreement wins over an unknown: it cannot become fine later.
static void test_disagreement_beats_unknown(void) {
    CU_ASSERT_TRUE(common_temperature({known("A", "b", 37), known("B", "b", 65), plain("C")}).how ==
                   Agreement::disagree);
}

static Enzyme with_inactivation(const char *name, int c, int minutes) {
    Enzyme e = known(name, "b", 37);
    e.inactivate_c = c;
    e.inactivate_min = minutes;
    return e;
}

static void test_inactivation_when_everyone_agrees(void) {
    std::vector<std::string> l = inactivation_lines({with_inactivation("A", 65, 20),
                                                     with_inactivation("B", 65, 20)});
    CU_ASSERT_EQUAL(l.size(), 1u);
    CU_ASSERT_EQUAL(l[0], std::string("Heat inactivate at 65 C for 20 min"));
}

static void test_inactivation_when_they_differ_is_per_enzyme(void) {
    std::vector<std::string> l = inactivation_lines({with_inactivation("A", 65, 20),
                                                     with_inactivation("B", 80, 20)});
    CU_ASSERT_EQUAL(l.size(), 2u);
    CU_ASSERT_TRUE(contains(l[0], "A: heat inactivate at 65 C for 20 min"));
    CU_ASSERT_TRUE(contains(l[1], "B: heat inactivate at 80 C for 20 min"));
}

// One that cannot be heat inactivated means heat is not the answer for the tube.
static void test_inactivation_when_one_cannot_be_heat_inactivated(void) {
    std::vector<std::string> l = inactivation_lines({with_inactivation("A", 65, 20),
                                                     with_inactivation("B", 0, 0)});
    bool says_cannot = false, says_heat_generally = false;
    for (const std::string &line : l) {
        if (contains(line, "B cannot be heat inactivated")) { says_cannot = true; }
        if (line == "Heat inactivate at 65 C for 20 min") { says_heat_generally = true; }
    }
    CU_ASSERT_TRUE(says_cannot);
    CU_ASSERT_TRUE(!says_heat_generally);
}

static void test_inactivation_unknown_is_said_not_guessed(void) {
    std::vector<std::string> l = inactivation_lines({plain("A")});
    CU_ASSERT_EQUAL(l.size(), 1u);
    CU_ASSERT_TRUE(contains(l[0], "A is not on file"));
    CU_ASSERT_TRUE(!contains(l[0], " C "));
}

// A temperature with no time (NEB's NEBuffer Performance Chart) prints the
// temperature and says the time is not on file - it does not invent one, and
// it is not folded into "unknown" since the temperature is a real fact.
static void test_inactivation_temperature_known_time_not(void) {
    Enzyme e = with_inactivation("A", 65, -1);
    std::vector<std::string> l = inactivation_lines({e});
    CU_ASSERT_EQUAL(l.size(), 1u);
    CU_ASSERT_TRUE(contains(l[0], "A at 65 C"));
    CU_ASSERT_TRUE(!contains(l[0], "not on file"));   // that phrasing is for fully unknown
    CU_ASSERT_TRUE(!contains(l[0], "for -1 min"));

    // Mixed with a fully known enzyme, each gets its own line rather than one
    // pretending to speak for both.
    std::vector<std::string> mixed =
        inactivation_lines({with_inactivation("A", 65, -1), with_inactivation("B", 65, 20)});
    CU_ASSERT_EQUAL(mixed.size(), 2u);
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
        {"NEB's table", {
            {"the 10, 25 and 50 uL rows", test_neb_small_volume_table},
            {"DNA limit between and beyond rows", test_dna_limit_between_and_beyond_the_rows},
        }},
        {"setup", {
            {"typical digest by hand", test_typical_setup_by_hand},
            {"genomic takes 20 units per ug", test_genomic_takes_twenty_units_per_ug},
            {"each enzyme gets the full dose", test_each_enzyme_gets_the_full_dose},
            {"stock sets the enzyme volume", test_stock_sets_the_enzyme_volume},
            {"default stock is the protocol's ratio", test_default_stock_is_the_protocols_ratio},
        }},
        {"refusals", {
            {"enzyme capped at 10%", test_enzyme_volume_is_capped_at_ten_percent},
            {"DNA that does not fit", test_dna_that_does_not_fit_is_refused},
            {"nonsense inputs", test_nonsense_inputs_are_refused},
        }},
        {"master mix and stop", {
            {"master mix by hand", test_master_mix_by_hand},
            {"needs a reaction", test_master_mix_needs_a_reaction},
            {"stop solution is a fifth", test_stop_solution_is_a_fifth_of_the_volume},
        }},
        {"conditions CSV", {
            {"a row is read", test_conditions_row_is_read},
            {"header only is empty", test_header_only_table_is_empty_not_an_error},
            {"columns in any order", test_conditions_columns_can_be_in_any_order},
            {"missing column", test_conditions_missing_column_rejected},
            {"bad number", test_conditions_bad_number_rejected},
            {"row without a source", test_conditions_row_without_a_source_rejected},
            {"inactivation both or neither", test_conditions_inactivation_is_both_or_neither},
            {"other fields checked", test_conditions_other_fields_are_checked},
            {"stock and time_saver can be unknown", test_conditions_stock_and_time_saver_can_be_unknown},
            {"inactivation temperature with no time", test_conditions_inactivation_temperature_with_no_time},
            {"duplicate enzyme", test_conditions_duplicate_enzyme_rejected},
            {"built-in table loads", test_builtin_table_loads},
        }},
        {"enzyme specs", {
            {"plain name, any case", test_plain_name_resolves_case_insensitively},
            {"nicking enzymes resolve", test_nicking_enzyme_names_resolve},
            {"unknown enzyme", test_unknown_enzyme_is_refused_with_the_way_out},
            {"conditions on the command line", test_conditions_on_the_command_line},
            {"inactivation on the command line", test_inactivation_on_the_command_line},
            {"inactivation temperature only on the command line",
             test_inactivation_temperature_only_on_the_command_line},
            {"malformed specs", test_malformed_specs_are_refused},
            {"command line vs the table row", test_command_line_overrides_the_table_row},
            {"enzyme only in the conditions table", test_enzyme_only_in_the_conditions_table_resolves},
        }},
        {"what the enzymes agree about", {
            {"buffer", test_buffer_agreement},
            {"temperature", test_temperature_agreement},
            {"disagreement beats unknown", test_disagreement_beats_unknown},
            {"inactivation, all agree", test_inactivation_when_everyone_agrees},
            {"inactivation, they differ", test_inactivation_when_they_differ_is_per_enzyme},
            {"inactivation, one cannot", test_inactivation_when_one_cannot_be_heat_inactivated},
            {"inactivation, unknown", test_inactivation_unknown_is_said_not_guessed},
            {"inactivation, temperature known time not", test_inactivation_temperature_known_time_not},
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
