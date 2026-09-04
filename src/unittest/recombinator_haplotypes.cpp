/** \file
 *
 * Unit tests for the kmer presence iterators in `Haplotypes::Subchain`.
 */

#include "catch.hpp"

#include "../recombinator_haplotypes.hpp"

#include <vector>

namespace vg {
namespace unittest {

//------------------------------------------------------------------------------

namespace {

/*
  Builds a subchain with the given kmer presence matrix, which is given as one
  string of '0' and '1' characters per sequence. The kmers and the sequences are
  dummy values, as the iterators only care about their number.
*/
Haplotypes::Subchain build_subchain(const std::vector<std::string>& matrix, size_t kmer_count) {
    Haplotypes::Subchain result;
    result.type = Haplotypes::Subchain::normal;
    result.start = gbwt::Node::encode(1, false);
    result.end = gbwt::Node::encode(2, false);

    for (size_t kmer_id = 0; kmer_id < kmer_count; kmer_id++) {
        result.kmers.push_back(kmer_id + 1);
    }
    result.kmer_counts = sdsl::int_vector<0>(kmer_count, 0, 64);
    for (size_t i = 0; i < matrix.size(); i++) {
        result.sequences.push_back({ static_cast<std::uint32_t>(i), 0 });
    }

    result.kmers_present = sdsl::bit_vector(matrix.size() * kmer_count, 0);
    for (size_t i = 0; i < matrix.size(); i++) {
        REQUIRE(matrix[i].length() == kmer_count);
        for (size_t j = 0; j < kmer_count; j++) {
            if (matrix[i][j] == '1') {
                result.kmers_present[i * kmer_count + j] = 1;
                result.kmer_counts[j]++;
            }
        }
    }

    return result;
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Iterating over the kmers in a sequence", "[haplotypes]") {
    std::vector<std::string> matrix {
        "00000",
        "11111",
        "10110",
        "01001"
    };
    Haplotypes::Subchain subchain = build_subchain(matrix, 5);

    SECTION("every kmer is visited once in order") {
        for (size_t i = 0; i < matrix.size(); i++) {
            std::vector<bool> found;
            subchain.for_each_kmer(i, [&](size_t kmer_id, bool is_present) {
                REQUIRE(kmer_id == found.size());
                found.push_back(is_present);
            });
            REQUIRE(found.size() == subchain.kmers.size());
            for (size_t j = 0; j < found.size(); j++) {
                bool should_be_present = (matrix[i][j] == '1');
                REQUIRE(found[j] == should_be_present);
            }
        }
    }

    SECTION("no kmers") {
        Haplotypes::Subchain empty = build_subchain({ "", "", "" }, 0);
        for (size_t i = 0; i < 3; i++) {
            size_t calls = 0;
            empty.for_each_kmer(i, [&](size_t, bool) { calls++; });
            REQUIRE(calls == 0);
        }
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Scoring the haplotypes in a subchain", "[haplotypes]") {
    std::vector<std::string> matrix {
        "00000",
        "11111",
        "10110",
        "01001"
    };
    Haplotypes::Subchain subchain = build_subchain(matrix, 5);

    // Arbitrary weights for the kmers.
    std::vector<double> weights { 1.0, -2.0, 0.5, 3.0, -1.5 };
    auto kmer_score = [&](size_t kmer_id, bool is_present) -> double {
        return (is_present ? 1.0 : -1.0) * weights[kmer_id];
    };

    SECTION("every haplotype is scored once in order") {
        std::vector<double> scores;
        subchain.score_sequences(
            kmer_score,
            [&](size_t haplotype_id, double score) {
                REQUIRE(haplotype_id == scores.size());
                scores.push_back(score);
            }
        );
        REQUIRE(scores.size() == subchain.sequences.size());
        for (size_t i = 0; i < scores.size(); i++) {
            double expected = 0.0;
            for (size_t j = 0; j < weights.size(); j++) {
                expected += (matrix[i][j] == '1' ? 1.0 : -1.0) * weights[j];
            }
            REQUIRE(scores[i] == expected);
        }
    }

    SECTION("every kmer of every haplotype is visited once") {
        std::vector<size_t> visits(subchain.kmers.size(), 0);
        subchain.score_sequences(
            [&](size_t kmer_id, bool) -> double { visits[kmer_id]++; return 0.0; },
            [&](size_t, double) {}
        );
        for (size_t kmer_id = 0; kmer_id < visits.size(); kmer_id++) {
            REQUIRE(visits[kmer_id] == subchain.sequences.size());
        }
    }

    SECTION("no kmers") {
        Haplotypes::Subchain empty = build_subchain({ "", "", "" }, 0);
        std::vector<double> scores;
        empty.score_sequences(
            [&](size_t, bool) -> double { REQUIRE(false); return 1.0; },
            [&](size_t haplotype_id, double score) {
                REQUIRE(haplotype_id == scores.size());
                scores.push_back(score);
            }
        );
        REQUIRE(scores.size() == 3);
        for (double score : scores) {
            REQUIRE(score == 0.0);
        }
    }

    SECTION("no haplotypes") {
        Haplotypes::Subchain empty = build_subchain({}, 5);
        size_t calls = 0;
        empty.score_sequences(
            [&](size_t, bool) -> double { REQUIRE(false); return 1.0; },
            [&](size_t, double) { calls++; }
        );
        REQUIRE(calls == 0);
    }
}

//------------------------------------------------------------------------------

} // namespace unittest
} // namespace vg
