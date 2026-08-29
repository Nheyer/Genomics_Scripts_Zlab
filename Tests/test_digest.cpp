// Unit tests for Restriction_Enzyme_Digest.
//
// Same arrangement as test_consensus.cpp: the program is one translation unit
// that owns its own main(), so we rename that main out of the way and pull the
// whole file in. The tests then call the real functions directly and the
// program source needs no changes to be testable.
//
// Ordering follows the upstream ROADMAP's instruction for a port: the
// invariants come first, because they are what catch a sign error or an
// off-by-one in the cut coordinate frame, and they are language independent.
// Several of them are written so they cannot be satisfied by re-deriving the
// implementation - see the comment on DEFINITE_MATCHES in particular.
#define main disabled_main
#include "Restriction_Enzyme_Digest.cpp"
#undef main

#include <CUnit/Basic.h>

#include <initializer_list>
#include <string>
#include <vector>

// ---------------------------------------------------------------- helpers --

static std::vector<int> cuts_of(const std::vector<SiteHit> &hits) {
    std::vector<int> out;
    for (const SiteHit &h : hits) { out.push_back(h.cut); }
    return out;
}

// One concrete ACGT spelling of a site, taking the first base each code allows.
// expand_ambiguity() would give the same string as its first element, but it
// refuses sites like XcmI's CCANNNNNNNNNTGG, and the bulk invariants below run
// over every enzyme in the table including that one.
static std::string first_concrete_site(const std::string &site) {
    static const char BASES[] = "ACGT";
    std::string out;
    for (char code : site) {
        uint8_t mask = mask_of(code);
        for (int bit = 0; bit < 4; bit++) {
            if (mask & (1u << bit)) { out += BASES[bit]; break; }
        }
    }
    return out;
}

// A deterministic pseudo-random sequence, so the bulk invariants run over
// something with real site density but stay reproducible across machines.
static std::string pseudo_random_sequence(int length, unsigned seed, const char *alphabet) {
    std::string out;
    out.reserve(length);
    unsigned state = seed;
    size_t span = std::string(alphabet).size();
    for (int i = 0; i < length; i++) {
        state = state * 1103515245u + 12345u;
        out += alphabet[(state >> 16) % span];
    }
    return out;
}

// ------------------------------------------------- invariant: IUPAC table --

// `--ambiguity definite` is the matching rule this tool has always used.
//
// This table is written out by hand rather than derived from IUPAC_BASES, so
// that changing the rule fails here instead of being mirrored by a
// re-derivation of itself. Read an entry as: a recognition site carrying this
// code matches exactly these bases in a target sequence. A code matches only
// where the target cannot resolve to anything the enzyme rejects - so the site
// is cut however the ambiguity resolves.
struct definite_row { char code; const char *matches; };
static const definite_row DEFINITE_MATCHES[] = {
    {'A', "A"},
    {'C', "C"},
    {'G', "G"},
    {'T', "T"},
    {'R', "AGR"},                // A or G
    {'Y', "CTY"},                // C or T
    {'S', "CGS"},                // C or G
    {'W', "ATW"},                // A or T
    {'K', "GTK"},                // G or T
    {'M', "ACM"},                // A or C
    {'B', "CGTYSKB"},            // not A
    {'D', "AGTRWKD"},            // not C
    {'H', "ACTYWMH"},            // not G
    {'V', "ACGRSMV"},            // not T
    {'N', "ACGTRYSWKMBDHVN"},    // anything
};
static const char PINNED_ALPHABET[] = "ACGTRYSWKMBDHVN";

// Which targets does this site code match, under one mode?
static std::string matches_for(char site_code, bool definite) {
    std::string out;
    for (int i = 0; PINNED_ALPHABET[i]; i++) {
        char target = PINNED_ALPHABET[i];
        if (match_base(mask_of(site_code), mask_of(target), definite)) { out += target; }
    }
    return out;
}

static std::string sorted_chars(const std::string &s) {
    std::string out = s;
    std::sort(out.begin(), out.end());
    return out;
}

static void test_table_covers_every_iupac_code(void) {
    std::string listed;
    for (const definite_row &row : DEFINITE_MATCHES) { listed += row.code; }
    CU_ASSERT_EQUAL(sorted_chars(listed), sorted_chars(PINNED_ALPHABET));
}

static void test_definite_matches_exactly_the_pinned_table(void) {
    for (const definite_row &row : DEFINITE_MATCHES) {
        std::string got = sorted_chars(matches_for(row.code, true));
        std::string want = sorted_chars(row.matches);
        if (got != want) {
            std::cerr << "\n  site code " << row.code << ": expected {" << want
                      << "} got {" << got << "}" << std::endl;
        }
        CU_ASSERT_EQUAL(got, want);
    }
}

// The anchor cases the docs promise: an enzyme N accepts a target N, an enzyme
// A does not; and enzyme T does not definitely match target Y, since Y may be C.
static void test_the_anchor_cases_the_docs_promise(void) {
    CU_ASSERT_TRUE(match_base(mask_of('N'), mask_of('N'), true));
    CU_ASSERT_FALSE(match_base(mask_of('A'), mask_of('N'), true));
    CU_ASSERT_FALSE(match_base(mask_of('T'), mask_of('Y'), true));
}

// Guards against the two modes being transposed: every definite match is also
// a possible match, and somewhere the two genuinely differ.
static void test_possible_is_strictly_looser_and_not_swapped(void) {
    int strictly_looser = 0;
    for (int i = 0; PINNED_ALPHABET[i]; i++) {
        char code = PINNED_ALPHABET[i];
        std::string definite = matches_for(code, true);
        std::string possible = matches_for(code, false);
        for (char c : definite) {
            CU_ASSERT_TRUE(possible.find(c) != std::string::npos);
        }
        if (definite.size() < possible.size()) { strictly_looser++; }
    }
    CU_ASSERT_TRUE(strictly_looser > 0);
}

// Enzyme T matches target Y under possible (Y may resolve to T); enzyme A does
// not, since Y is C or T and shares no base with A.
static void test_possible_anchor_cases(void) {
    CU_ASSERT_TRUE(match_base(mask_of('T'), mask_of('Y'), false));
    CU_ASSERT_FALSE(match_base(mask_of('A'), mask_of('Y'), false));
}

// A character that is not an IUPAC code matches under neither mode. Without
// this, mask 0 would be a subset of everything and match the whole alphabet.
static void test_non_iupac_target_never_matches(void) {
    CU_ASSERT_FALSE(match_base(mask_of('N'), mask_of('-'), true));
    CU_ASSERT_FALSE(match_base(mask_of('N'), mask_of('-'), false));
    CU_ASSERT_FALSE(match_base(mask_of('N'), mask_of('x'), true));
    CU_ASSERT_FALSE(match_base(mask_of('N'), mask_of('\0'), false));
}

// ------------------------------------------- invariant: fragments sum to N --

// The strongest check available in bulk: whatever the cuts are, the pieces have
// to add back up to the molecule. This is what catches sign errors and
// off-by-ones in the coordinate frame.
static void assert_fragments_sum(const std::string &seq, bool circular) {
    const EnzymeTables &db = builtin_tables();
    for (const Enzyme &enzyme : db.enzymes) {
        std::vector<SiteHit> hits = find_sites(seq, enzyme, AMBIGUITY_DEFINITE, circular);
        int dropped = 0;
        std::vector<int> cuts = normalise_cuts(cuts_of(hits), (int) seq.size(), circular, &dropped);
        std::vector<int> fragments = circular
            ? digest_circular((int) seq.size(), cuts, 1)
            : digest_linear((int) seq.size(), cuts, 1);
        long long total = 0;
        for (int f : fragments) { total += f; }
        if (total != (long long) seq.size()) {
            std::cerr << "\n  " << enzyme.name << " (" << to_neb_notation(enzyme) << ") "
                      << (circular ? "circular" : "linear") << ": fragments total " << total
                      << ", sequence is " << seq.size() << std::endl;
        }
        CU_ASSERT_EQUAL(total, (long long) seq.size());
    }
}

static void test_fragments_sum_to_length_linear(void) {
    assert_fragments_sum(pseudo_random_sequence(3000, 20260827u, "ACGT"), false);
}

static void test_fragments_sum_to_length_circular(void) {
    assert_fragments_sum(pseudo_random_sequence(3000, 20260827u, "ACGT"), true);
}

// Ambiguity codes in the target must not break the accounting either.
static void test_fragments_sum_with_ambiguous_sequence(void) {
    assert_fragments_sum(pseudo_random_sequence(1500, 99u, "ACGTRYSWKMBDHVN"), false);
    assert_fragments_sum(pseudo_random_sequence(1500, 99u, "ACGTRYSWKMBDHVN"), true);
}

// ------------------------------- invariant: overhangs are self-complementary --

// For a palindromic enzyme that cuts inside its own site, the two cut offsets
// have to mirror about the middle. That is forced by the symmetry of the
// molecule, so it is independent of the implementation and of REBASE.
static void test_palindromic_in_site_cutters_have_mirrored_cuts(void) {
    const EnzymeTables &db = builtin_tables();
    int checked = 0;
    for (const Enzyme &e : db.enzymes) {
        if (!is_palindromic(e)) { continue; }
        int length = (int) e.site.size();
        if (e.cut_top < 0 || e.cut_top > length) { continue; }   // cuts outside its site
        if (e.cut_bottom < 0 || e.cut_bottom > length) { continue; }
        if (e.cut_bottom != length - e.cut_top) {
            std::cerr << "\n  " << e.name << " (" << to_neb_notation(e) << "): cut_bottom "
                      << e.cut_bottom << " is not " << length << " - " << e.cut_top << std::endl;
        }
        CU_ASSERT_EQUAL(e.cut_bottom, length - e.cut_top);
        checked++;
    }
    CU_ASSERT_TRUE(checked > 50);   // the check is worthless if it matched nothing
}

// Every enzyme in the table must cut a sequence built from its own site.
static void test_every_enzyme_cuts_its_own_site(void) {
    const EnzymeTables &db = builtin_tables();
    for (const Enzyme &e : db.enzymes) {
        // Pad generously so a type IIS cut downstream still lands in-sequence.
        std::string filler(40, 'A');
        std::string concrete = first_concrete_site(e.site);
        std::string seq = filler + concrete + filler;
        std::vector<SiteHit> hits = find_sites(seq, e, AMBIGUITY_DEFINITE, false);
        if (hits.empty()) {
            std::cerr << "\n  " << e.name << " (" << to_neb_notation(e)
                      << ") found no site in " << seq << std::endl;
        }
        CU_ASSERT_TRUE(!hits.empty());
    }
}

// ...and must find it on the reverse complement strand too.
static void test_every_enzyme_finds_its_site_reverse_complemented(void) {
    const EnzymeTables &db = builtin_tables();
    for (const Enzyme &e : db.enzymes) {
        std::string filler(40, 'A');
        std::string concrete = first_concrete_site(e.site);
        std::string seq = filler + reverse_complement(concrete) + filler;
        std::vector<SiteHit> hits = find_sites(seq, e, AMBIGUITY_DEFINITE, false);
        if (hits.empty()) {
            std::cerr << "\n  " << e.name << " (" << to_neb_notation(e)
                      << ") missed its reverse complement" << std::endl;
        }
        CU_ASSERT_TRUE(!hits.empty());
    }
}

// -------------------------------------------------- the enzyme table itself --

static void test_table_sizes(void) {
    const EnzymeTables &db = builtin_tables();
    CU_ASSERT_EQUAL(db.enzymes.size(), 262u);
    CU_ASSERT_EQUAL(db.nicking.size(), 13u);
    CU_ASSERT_EQUAL(db.dual_cut.size(), 12u);
}

static void test_every_site_is_valid_iupac(void) {
    const EnzymeTables &db = builtin_tables();
    for (const Enzyme &e : db.enzymes) {
        CU_ASSERT_TRUE(!e.site.empty());
        for (char c : e.site) { CU_ASSERT_TRUE(mask_of(c) != 0); }
    }
}

// Known overhangs, from the enzymes' published end structures.
static void test_known_overhangs(void) {
    const EnzymeTables &db = builtin_tables();
    struct expectation { const char *name; int overhang; };
    static const expectation WANT[] = {
        {"EcoRI", 4},      // 5' AATT
        {"PstI", -4},      // 3' TGCA
        {"EcoRV", 0},      // blunt
        {"BsaI", 4},       // 5' four base, downstream
        {"KpnI", -4},
        {"SmaI", 0},
    };
    for (const expectation &want : WANT) {
        bool found = false;
        for (const Enzyme &e : db.enzymes) {
            if (e.name == want.name) {
                CU_ASSERT_EQUAL(overhang(e), want.overhang);
                found = true;
                break;
            }
        }
        CU_ASSERT_TRUE(found);
    }
}

// The NEB notation must survive a round trip through the parser for every
// enzyme in the table.
static void test_notation_round_trips_over_whole_db(void) {
    const EnzymeTables &db = builtin_tables();
    for (const Enzyme &e : db.enzymes) {
        std::string site;
        int top = 0, bottom = 0;
        parse_neb_notation(to_neb_notation(e), &site, &top, &bottom);
        CU_ASSERT_EQUAL(site, e.site);
        CU_ASSERT_EQUAL(top, e.cut_top);
        CU_ASSERT_EQUAL(bottom, e.cut_bottom);
    }
}

// ------------------------------------------------------- notation parsing --

static void test_parse_caret_notation(void) {
    std::string site;
    int top = 0, bottom = 0;
    parse_neb_notation("G^AATTC", &site, &top, &bottom);
    CU_ASSERT_EQUAL(site, std::string("GAATTC"));
    CU_ASSERT_EQUAL(top, 1);
    CU_ASSERT_EQUAL(bottom, 5);
}

// MboI cuts immediately before its site, so cut_top is 0.
static void test_parse_leading_caret(void) {
    std::string site;
    int top = 0, bottom = 0;
    parse_neb_notation("^GATC", &site, &top, &bottom);
    CU_ASSERT_EQUAL(site, std::string("GATC"));
    CU_ASSERT_EQUAL(top, 0);
    CU_ASSERT_EQUAL(bottom, 4);
}

static void test_parse_offset_notation(void) {
    std::string site;
    int top = 0, bottom = 0;
    parse_neb_notation("GGTCTC(1/5)", &site, &top, &bottom);
    CU_ASSERT_EQUAL(site, std::string("GGTCTC"));
    CU_ASSERT_EQUAL(top, 7);
    CU_ASSERT_EQUAL(bottom, 11);
}

static void test_parse_negative_offsets(void) {
    std::string site;
    int top = 0, bottom = 0;
    parse_neb_notation("GAATTC(-5/-1)", &site, &top, &bottom);
    CU_ASSERT_EQUAL(site, std::string("GAATTC"));
    CU_ASSERT_EQUAL(top, 1);
    CU_ASSERT_EQUAL(bottom, 5);
}

template <typename Fn>
static bool throws_invalid(Fn fn) {
    try { fn(); } catch (const std::exception &) { return true; }
    return false;
}

static void test_dual_cut_notation_rejected(void) {
    CU_ASSERT_TRUE(throws_invalid([] {
        std::string s; int t, b;
        parse_neb_notation("(10/15)ACNNNNGTAYC(12/7)", &s, &t, &b);
    }));
}

static void test_two_carets_rejected(void) {
    CU_ASSERT_TRUE(throws_invalid([] {
        std::string s; int t, b;
        parse_neb_notation("G^AA^TTC", &s, &t, &b);
    }));
}

static void test_site_without_cut_rejected(void) {
    CU_ASSERT_TRUE(throws_invalid([] {
        std::string s; int t, b;
        parse_neb_notation("GAATTC", &s, &t, &b);
    }));
}

// ------------------------------------------------------- ^/_ CSV notation --

static void test_split_markers_both_inside(void) {
    std::string site;
    std::vector<int> top, bottom;
    split_cut_markers("G^AATT_C", &site, &top, &bottom);
    CU_ASSERT_EQUAL(site, std::string("GAATTC"));
    CU_ASSERT_EQUAL(top.size(), 1u);
    CU_ASSERT_EQUAL(bottom.size(), 1u);
    CU_ASSERT_EQUAL(top[0], 1);
    CU_ASSERT_EQUAL(bottom[0], 5);
}

// BsaI: the padding Ns are spacing, not recognition, so they are trimmed and
// the offsets shift to match. This is the case most likely to be mis-ported.
static void test_split_markers_trims_padding_ns(void) {
    std::string site;
    std::vector<int> top, bottom;
    split_cut_markers("GGTCTCN^NNNN_", &site, &top, &bottom);
    CU_ASSERT_EQUAL(site, std::string("GGTCTC"));
    CU_ASSERT_EQUAL(top[0], 7);
    CU_ASSERT_EQUAL(bottom[0], 11);
}

// AatII, a 3' overhang: the bottom cut comes before the top one.
static void test_split_markers_three_prime_overhang(void) {
    std::string site;
    std::vector<int> top, bottom;
    split_cut_markers("G_ACGT^C", &site, &top, &bottom);
    CU_ASSERT_EQUAL(site, std::string("GACGTC"));
    CU_ASSERT_EQUAL(top[0], 5);
    CU_ASSERT_EQUAL(bottom[0], 1);
}

// --------------------------------------------------------- spec parsing ----

static void test_builtin_lookup_is_case_insensitive(void) {
    Enzyme lower = parse_enzyme_spec("ecori");
    Enzyme upper = parse_enzyme_spec("ECORI");
    CU_ASSERT_EQUAL(lower.name, std::string("EcoRI"));   // reported canonically
    CU_ASSERT_EQUAL(upper.name, std::string("EcoRI"));
    CU_ASSERT_EQUAL(lower.site, std::string("GAATTC"));
}

static void test_custom_spec_name_seq_offset(void) {
    Enzyme e = parse_enzyme_spec("PmeI:GTTTAAAC:4");
    CU_ASSERT_EQUAL(e.name, std::string("PmeI"));
    CU_ASSERT_EQUAL(e.cut_top, 4);
    CU_ASSERT_EQUAL(e.cut_bottom, 4);
}

static void test_custom_spec_two_offsets(void) {
    Enzyme e = parse_enzyme_spec("MyI:GGTCTC:7:11");
    CU_ASSERT_EQUAL(e.cut_top, 7);
    CU_ASSERT_EQUAL(e.cut_bottom, 11);
}

static void test_offset_out_of_range_rejected(void) {
    CU_ASSERT_TRUE(throws_invalid([] { parse_enzyme_spec("MyI:ACGT:99"); }));
}

static void test_non_integer_offset_rejected(void) {
    CU_ASSERT_TRUE(throws_invalid([] { parse_enzyme_spec("MyI:ACGT:x"); }));
}

static void test_unknown_enzyme_rejected(void) {
    CU_ASSERT_TRUE(throws_invalid([] { parse_enzyme_spec("NoSuchEnzyme"); }));
}

// A name we know but cannot simulate has to say why, not "unknown enzyme".
static void test_nicking_enzyme_named_not_unknown(void) {
    try {
        parse_enzyme_spec("Nb.BsmI");
        CU_FAIL("expected Nb.BsmI to be refused");
    } catch (const std::exception &e) {
        CU_ASSERT_TRUE(std::string(e.what()).find("nicks one strand") != std::string::npos);
    }
}

static void test_dual_cut_enzyme_named_not_unknown(void) {
    try {
        parse_enzyme_spec("BcgI");
        CU_FAIL("expected BcgI to be refused");
    } catch (const std::exception &e) {
        CU_ASSERT_TRUE(std::string(e.what()).find("both sides") != std::string::npos);
    }
}

// The legacy NAME:SEQ:OFFSET form cannot describe an asymmetric or outside-site
// cutter, and refuses rather than silently dropping the bottom-strand cut.
static void test_legacy_spec_refuses_to_downgrade(void) {
    CU_ASSERT_TRUE(throws_invalid([] { to_legacy_spec(parse_enzyme_spec("BsaI")); }));
    CU_ASSERT_EQUAL(to_legacy_spec(parse_enzyme_spec("EcoRI")),
                    std::string("EcoRI:GAATTC:1"));
}

// ------------------------------------------------------------ site finding --

static void test_two_sites(void) {
    Enzyme e = parse_enzyme_spec("EcoRI");
    std::vector<SiteHit> hits = find_sites("AAGAATTCAAAAGAATTCAA", e);
    CU_ASSERT_EQUAL(hits.size(), 2u);
    CU_ASSERT_EQUAL(hits[0].start, 2);
    CU_ASSERT_EQUAL(hits[1].start, 12);
}

static void test_no_site(void) {
    Enzyme e = parse_enzyme_spec("EcoRI");
    CU_ASSERT_EQUAL(find_sites("AAAACCCCGGGGTTTT", e).size(), 0u);
}

// Overlapping recognition sites must all be reported, not just the first.
static void test_overlapping_sites(void) {
    Enzyme e = parse_enzyme_spec("MyI:AAAA:1");
    std::vector<SiteHit> hits = find_sites("AAAAAA", e);
    CU_ASSERT_EQUAL(hits.size(), 3u);
}

// A fully symmetric enzyme is reported once per site, not twice.
static void test_palindromic_site_reported_once(void) {
    Enzyme e = parse_enzyme_spec("EcoRI");
    CU_ASSERT_EQUAL(find_sites("AAGAATTCAA", e).size(), 1u);
}

// A non-palindromic enzyme is found on the bottom strand, and cuts upstream of
// where it was found - the mirror rule.
static void test_non_palindromic_site_on_bottom_strand(void) {
    Enzyme e = parse_enzyme_spec("BsaI");           // GGTCTC(1/5), site len 6
    std::string rc = reverse_complement("GGTCTC");  // GAGACC
    std::string seq = std::string(20, 'A') + rc + std::string(20, 'A');
    std::vector<SiteHit> hits = find_sites(seq, e);
    CU_ASSERT_EQUAL(hits.size(), 1u);
    CU_ASSERT_EQUAL(hits[0].strand, '-');
    // start + len(site) - cut_bottom = 20 + 6 - 11 = 15, i.e. upstream.
    CU_ASSERT_EQUAL(hits[0].cut, 15);
}

static void test_non_palindromic_site_on_both_strands(void) {
    Enzyme e = parse_enzyme_spec("BsaI");
    std::string seq = std::string(20, 'A') + "GGTCTC" + std::string(20, 'A') +
                      reverse_complement("GGTCTC") + std::string(20, 'A');
    std::vector<SiteHit> hits = find_sites(seq, e);
    CU_ASSERT_EQUAL(hits.size(), 2u);
    int bottom = 0;
    for (const SiteHit &h : hits) {
        if (h.strand == '-') { bottom++; }
    }
    CU_ASSERT_EQUAL(bottom, 1);
}

// is_symmetric is strictly stronger than is_palindromic: it also needs the cuts
// to mirror, and it is what decides whether the reverse orientation is searched.
// Every built-in palindromic enzyme happens to be symmetric too, so this needs
// a custom spec to distinguish them.
static void test_palindromic_but_asymmetric_searches_both_strands(void) {
    Enzyme e = parse_enzyme_spec("MyI:GAATTC:1:4");
    CU_ASSERT_TRUE(is_palindromic(e));
    CU_ASSERT_FALSE(is_symmetric(e));
    std::vector<SiteHit> hits = find_sites("AAGAATTCAA", e);
    // Both orientations match the same palindromic site, and they cut in
    // different places, so both survive the (start, cut) dedup.
    CU_ASSERT_EQUAL(hits.size(), 2u);
}

// --------------------------------------------------- circular origin sites --

// A site straddling the origin is invisible on a linear scan and found on a
// circular one. This behaviour is promised by the docs, not incidental.
static void test_origin_spanning_site_found_only_when_circular(void) {
    // GAATTC split across the join: sequence ends "GAA", begins "TTC".
    std::string seq = "TTCAAAAAAAAAAAAAAAAGAA";
    Enzyme e = parse_enzyme_spec("EcoRI");
    CU_ASSERT_EQUAL(find_sites(seq, e, AMBIGUITY_DEFINITE, false).size(), 0u);
    std::vector<SiteHit> hits = find_sites(seq, e, AMBIGUITY_DEFINITE, true);
    CU_ASSERT_EQUAL(hits.size(), 1u);
    CU_ASSERT_EQUAL(hits[0].start, 19);
}

// The wrap must not double count a site that is already wholly inside.
static void test_circular_scan_does_not_double_count(void) {
    Enzyme e = parse_enzyme_spec("EcoRI");
    std::string seq = "AAGAATTCAA" + std::string(30, 'C');
    CU_ASSERT_EQUAL(find_sites(seq, e, AMBIGUITY_DEFINITE, true).size(), 1u);
}

// ------------------------------------------------------- cut normalisation --

static void test_linear_drops_cuts_outside_sequence(void) {
    int dropped = 0;
    std::vector<int> kept = normalise_cuts({-3, 0, 5, 10, 25}, 10, false, &dropped);
    // 0 and 10 are the sequence ends: cutting there yields no second piece.
    CU_ASSERT_EQUAL(kept.size(), 1u);
    CU_ASSERT_EQUAL(kept[0], 5);
    CU_ASSERT_EQUAL(dropped, 4);
}

// Python's % is a floor mod and never negative; C++'s is. A bottom-strand type
// IIS hit near the origin is exactly where that difference would show up.
static void test_circular_wraps_negative_cuts_the_python_way(void) {
    int dropped = 0;
    std::vector<int> kept = normalise_cuts({-1, -10, 3, 25}, 10, true, &dropped);
    // Nothing is ever dropped on a circle, and the results come back sorted:
    // -10 -> 0, 3 -> 3, 25 -> 5, -1 -> 9.
    CU_ASSERT_EQUAL(dropped, 0);
    CU_ASSERT_EQUAL(kept.size(), 4u);
    CU_ASSERT_EQUAL(kept[0], 0);
    CU_ASSERT_EQUAL(kept[1], 3);
    CU_ASSERT_EQUAL(kept[2], 5);
    CU_ASSERT_EQUAL(kept[3], 9);
}

// Spelled out separately so the arithmetic above cannot be read as a typo.
static void test_circular_modulo_values(void) {
    int dropped = 0;
    std::vector<int> kept = normalise_cuts({-1}, 10, true, &dropped);
    CU_ASSERT_EQUAL(kept.size(), 1u);
    CU_ASSERT_EQUAL(kept[0], 9);
    kept = normalise_cuts({25}, 10, true, &dropped);
    CU_ASSERT_EQUAL(kept[0], 5);
}

// A type IIS site near the end of a linear molecule is recognised but not cut.
static void test_type_iis_site_near_linear_end_is_not_cut(void) {
    Enzyme e = parse_enzyme_spec("BsaI");
    std::string seq = std::string(10, 'A') + "GGTCTC";   // cut would be at 17, len 16
    std::vector<SiteHit> hits = find_sites(seq, e);
    CU_ASSERT_TRUE(!hits.empty());
    int dropped = 0;
    std::vector<int> cuts = normalise_cuts(cuts_of(hits), (int) seq.size(), false, &dropped);
    CU_ASSERT_EQUAL(cuts.size(), 0u);
    CU_ASSERT_TRUE(dropped > 0);
}

// ------------------------------------------------------------- digestion ---

static void test_digest_linear_basic(void) {
    std::vector<int> fragments = digest_linear(100, {30, 70}, 1);
    CU_ASSERT_EQUAL(fragments.size(), 3u);
    CU_ASSERT_EQUAL(fragments[0], 30);
    CU_ASSERT_EQUAL(fragments[1], 40);
    CU_ASSERT_EQUAL(fragments[2], 30);
}

static void test_digest_linear_no_cuts(void) {
    std::vector<int> fragments = digest_linear(100, {}, 1);
    CU_ASSERT_EQUAL(fragments.size(), 1u);
    CU_ASSERT_EQUAL(fragments[0], 100);
}

static void test_digest_linear_min_fragment_filter(void) {
    std::vector<int> fragments = digest_linear(100, {10, 15, 70}, 20);
    for (int f : fragments) { CU_ASSERT_TRUE(f >= 20); }
}

static void test_digest_circular_basic(void) {
    std::vector<int> fragments = digest_circular(100, {30, 70}, 1);
    CU_ASSERT_EQUAL(fragments.size(), 2u);
    CU_ASSERT_EQUAL(fragments[0], 40);
    CU_ASSERT_EQUAL(fragments[1], 60);   // wraps the origin
}

// One cut on a circle linearises it: one fragment, the whole molecule.
static void test_digest_circular_single_cut(void) {
    std::vector<int> fragments = digest_circular(100, {30}, 1);
    CU_ASSERT_EQUAL(fragments.size(), 1u);
    CU_ASSERT_EQUAL(fragments[0], 100);
}

static void test_digest_circular_no_cuts(void) {
    std::vector<int> fragments = digest_circular(100, {}, 1);
    CU_ASSERT_EQUAL(fragments.size(), 1u);
    CU_ASSERT_EQUAL(fragments[0], 100);
}

// --------------------------------------------------------------- gaps ------

// A gap is not a base: the molecule reads straight through it, so AAT-ATT is
// really an SspI site once the gap is gone.
static void test_strip_gaps_reveals_a_site(void) {
    int removed = 0;
    std::string seq = strip_gaps("AAT-ATT", &removed);
    CU_ASSERT_EQUAL(removed, 1);
    CU_ASSERT_EQUAL(seq, std::string("AATATT"));
    Enzyme e = parse_enzyme_spec("SspI");
    CU_ASSERT_EQUAL(find_sites(seq, e).size(), 1u);
}

static void test_strip_gaps_counts_every_gap(void) {
    int removed = 0;
    std::string seq = strip_gaps("--ACGT--ACGT--", &removed);
    CU_ASSERT_EQUAL(removed, 6);
    CU_ASSERT_EQUAL(seq, std::string("ACGTACGT"));
}

// ------------------------------------------------------ ambiguity in target --

// An ambiguous target base only cuts definitely when every resolution cuts.
static void test_ambiguous_target_definite_versus_possible(void) {
    Enzyme e = parse_enzyme_spec("EcoRI");             // GAATTC
    std::string seq = "AAGAATTNAA";                    // last site base is N
    CU_ASSERT_EQUAL(find_sites(seq, e, AMBIGUITY_DEFINITE, false).size(), 0u);
    CU_ASSERT_EQUAL(find_sites(seq, e, AMBIGUITY_POSSIBLE, false).size(), 1u);
}

// An enzyme whose own site carries a code matches every base that code covers.
static void test_ambiguous_enzyme_finds_all_variants(void) {
    Enzyme e = parse_enzyme_spec("HincII");            // GTYRAC
    for (const std::string &variant : expand_ambiguity("GTYRAC")) {
        std::string seq = "AA" + variant + "AA";
        CU_ASSERT_EQUAL(find_sites(seq, e).size(), 1u);
    }
}

static void test_expand_ambiguity_enumerates_in_acgt_order(void) {
    std::vector<std::string> variants = expand_ambiguity("GTYRAC");
    CU_ASSERT_EQUAL(variants.size(), 4u);
    CU_ASSERT_EQUAL(variants[0], std::string("GTCAAC"));
    CU_ASSERT_EQUAL(variants[1], std::string("GTCGAC"));
    CU_ASSERT_EQUAL(variants[2], std::string("GTTAAC"));
    CU_ASSERT_EQUAL(variants[3], std::string("GTTGAC"));
}

static void test_expand_ambiguity_caps_runaway_sites(void) {
    // XcmI's nine Ns alone would give 262144 strings.
    CU_ASSERT_TRUE(throws_invalid([] { expand_ambiguity("CCANNNNNNNNNTGG"); }));
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
        {"invariant: definite matching is pinned", {
            {"table covers every IUPAC code", test_table_covers_every_iupac_code},
            {"definite matches the pinned table", test_definite_matches_exactly_the_pinned_table},
            {"the anchor cases the docs promise", test_the_anchor_cases_the_docs_promise},
            {"possible is looser, not swapped", test_possible_is_strictly_looser_and_not_swapped},
            {"possible anchor cases", test_possible_anchor_cases},
            {"non IUPAC target never matches", test_non_iupac_target_never_matches},
        }},
        {"invariant: fragments sum to sequence length", {
            {"linear, whole database", test_fragments_sum_to_length_linear},
            {"circular, whole database", test_fragments_sum_to_length_circular},
            {"ambiguous sequence", test_fragments_sum_with_ambiguous_sequence},
        }},
        {"invariant: enzyme geometry", {
            {"palindromic in site cutters mirror", test_palindromic_in_site_cutters_have_mirrored_cuts},
            {"every enzyme cuts its own site", test_every_enzyme_cuts_its_own_site},
            {"every enzyme finds its reverse complement", test_every_enzyme_finds_its_site_reverse_complemented},
        }},
        {"enzyme table", {
            {"table sizes", test_table_sizes},
            {"every site is valid IUPAC", test_every_site_is_valid_iupac},
            {"known overhangs", test_known_overhangs},
            {"notation round trips over whole db", test_notation_round_trips_over_whole_db},
        }},
        {"NEB notation", {
            {"caret notation", test_parse_caret_notation},
            {"leading caret (MboI)", test_parse_leading_caret},
            {"offset notation", test_parse_offset_notation},
            {"negative offsets", test_parse_negative_offsets},
            {"dual cut rejected", test_dual_cut_notation_rejected},
            {"two carets rejected", test_two_carets_rejected},
            {"site without a cut rejected", test_site_without_cut_rejected},
        }},
        {"CSV ^/_ notation", {
            {"both markers inside the site", test_split_markers_both_inside},
            {"trims padding Ns (BsaI)", test_split_markers_trims_padding_ns},
            {"3' overhang (AatII)", test_split_markers_three_prime_overhang},
        }},
        {"enzyme specs", {
            {"builtin lookup is case insensitive", test_builtin_lookup_is_case_insensitive},
            {"NAME:SEQ:OFFSET", test_custom_spec_name_seq_offset},
            {"NAME:SEQ:TOP:BOTTOM", test_custom_spec_two_offsets},
            {"offset out of range rejected", test_offset_out_of_range_rejected},
            {"non integer offset rejected", test_non_integer_offset_rejected},
            {"unknown enzyme rejected", test_unknown_enzyme_rejected},
            {"nicking enzyme named, not unknown", test_nicking_enzyme_named_not_unknown},
            {"dual cut enzyme named, not unknown", test_dual_cut_enzyme_named_not_unknown},
            {"legacy spec refuses to downgrade", test_legacy_spec_refuses_to_downgrade},
        }},
        {"site finding", {
            {"two sites", test_two_sites},
            {"no site", test_no_site},
            {"overlapping sites all reported", test_overlapping_sites},
            {"palindromic site reported once", test_palindromic_site_reported_once},
            {"bottom strand site cuts upstream", test_non_palindromic_site_on_bottom_strand},
            {"sites on both strands", test_non_palindromic_site_on_both_strands},
            {"palindromic but asymmetric searches both", test_palindromic_but_asymmetric_searches_both_strands},
        }},
        {"circular molecules", {
            {"origin spanning site needs --circular", test_origin_spanning_site_found_only_when_circular},
            {"wrap does not double count", test_circular_scan_does_not_double_count},
        }},
        {"cut normalisation", {
            {"linear drops out of range cuts", test_linear_drops_cuts_outside_sequence},
            {"circular wraps negative cuts", test_circular_wraps_negative_cuts_the_python_way},
            {"circular modulo values", test_circular_modulo_values},
            {"type IIS near a linear end is not cut", test_type_iis_site_near_linear_end_is_not_cut},
        }},
        {"digestion", {
            {"linear basic", test_digest_linear_basic},
            {"linear with no cuts", test_digest_linear_no_cuts},
            {"linear min fragment filter", test_digest_linear_min_fragment_filter},
            {"circular basic", test_digest_circular_basic},
            {"circular single cut linearises", test_digest_circular_single_cut},
            {"circular with no cuts", test_digest_circular_no_cuts},
        }},
        {"gaps", {
            {"stripping reveals a site", test_strip_gaps_reveals_a_site},
            {"counts every gap", test_strip_gaps_counts_every_gap},
        }},
        {"ambiguity codes", {
            {"definite versus possible in the target", test_ambiguous_target_definite_versus_possible},
            {"enzyme site codes match all variants", test_ambiguous_enzyme_finds_all_variants},
            {"expansion is in ACGT order", test_expand_ambiguity_enumerates_in_acgt_order},
            {"expansion caps runaway sites", test_expand_ambiguity_caps_runaway_sites},
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
