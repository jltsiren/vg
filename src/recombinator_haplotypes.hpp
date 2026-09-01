#ifndef VG_RECOMBINATOR_HAPLOTYPES_HPP_INCLUDED
#define VG_RECOMBINATOR_HAPLOTYPES_HPP_INCLUDED

/** \file recombinator_haplotypes.hpp
 * Haplotype information files used in haplotype sampling.
 */

#include "gbwt_helper.hpp"
#include "gbwtgraph_helper.hpp"
#include "hash_map.hpp"

#include <functional>
#include <iostream>

namespace vg {

//------------------------------------------------------------------------------

// FIXME: version 7 with compressed k-mer matrices. and tests

// FIXME: the haplotype scoring interface must be parameterized over score types (int64_t, double)

// FIXME: version 7: each top-level chain should know the default scoring model for that chain

/**
 * A representation of the haplotypes in a graph.
 *
 * The graph is partitioned into top-level chains, which are further partitioned
 * into subchains. Each subchain contains a set of kmers and a collection of
 * sequences. Each sequence is defined by a bitvector marking the kmers that are
 * present.
 *
 * At the moment, the kmers are minimizers with a single occurrence in the graph.
 * The requirement is that each kmer is specific to a single subchain and does
 * not occur anywhere else in either orientation. (If no haplotype crosses a
 * snarl, that snarl is broken into a suffix and a prefix, and those subchains
 * may share kmers.)
 *
 * NOTE: This assumes that the top-level chains are linear, not cyclical.
 *
 * Versions:
 * * Version 6: Tags for storing metadata, such as graph names. Compatible with
 *   version 5.
 *
 * * Version 5: Every path in the graph is assigned to a construction job.
 *   This allows including reference paths that do not visit any snarls in the
 *   sampled graph. Not compatible with earlier versions.
 *
 * * Version 4: Subchains can have fragmented haplotypes instead of a single
 *   GBWT sequence always crossing from start to end. Compatible with version 3.
 *
 * * Version 3: Subchains use smaller integers when possible. Compatible with
 *   version 2.
 *
 * * Version 2: Top-level chains include a contig name. Compatible with version 1.
 *
 * * Version 1: Initial version.
 */
class Haplotypes {
public:
    /// Default target length of a subchain.
    constexpr static size_t SUBCHAIN_LENGTH = 10000;

    /// Block size (in kmers) for reading KFF files.
    constexpr static size_t KFF_BLOCK_SIZE = 1000000;

    /// The amount of progress information that should be printed to stderr.
    enum Verbosity : size_t {
        /// No progress information.
        verbosity_silent = 0,

        /// Basic information.
        verbosity_basic = 1,

        /// Basic information and detailed statistics.
        verbosity_detailed = 2,

        /// Basic information, detailed statistics, and debug information.
        verbosity_debug = 3,

        /// Hidden level; potentially tens of thousands of lines of debugging information.
        verbosity_extra_debug = 4
    };

    /// Header of the serialized file.
    struct Header {
        constexpr static std::uint32_t MAGIC_NUMBER = 0x4C504148; // "HAPL"
        constexpr static std::uint32_t VERSION = 6;
        constexpr static std::uint32_t VERSION_WITH_TAGS = 6;
        constexpr static std::uint32_t MIN_VERSION = 5;
        constexpr static std::uint64_t DEFAULT_K = 29;

        /// A magic number that identifies the file.
        std::uint32_t magic_number = MAGIC_NUMBER;

        /// Version of the file.
        std::uint32_t version = VERSION;

        /// Number of top-level chains in the graph.
        std::uint64_t top_level_chains = 0;

        /// Number of GBWT construction jobs for the chains.
        std::uint64_t construction_jobs = 0;

        /// Total number of subchains in all chains.
        std::uint64_t total_subchains = 0;

        /// Total number of kmers in all subchains.
        std::uint64_t total_kmers = 0;

        /// Length of the kmers.
        std::uint64_t k = DEFAULT_K;
    };

    /// A GBWT sequence as (sequence identifier, offset in a node).
    typedef std::pair<gbwt::size_type, gbwt::size_type> sequence_type;

    /// A more space-efficient representation of `sequence_type`.
    typedef std::pair<std::uint32_t, std::uint32_t> compact_sequence_type;

    /// Representation of a subchain.
    struct Subchain {
        /// Subchain types.
        enum subchain_t : std::uint64_t {
            /// Normal subchain with two boundary nodes.
            normal = 0,

            /// A prefix with only an end node.
            prefix = 1,

            /// A suffix with only a start node.
            suffix = 2,

            /// A full haplotype with no boundary nodes.
            full_haplotype = 3
        };

        /// An encoded kmer.
        typedef gbwtgraph::Key64::value_type kmer_type;

        /// The type of this subchain.
        subchain_t type;

        /// Boundary nodes, or `gbwt::ENDMARKER` if not present.
        gbwt::node_type start, end;

        /// A vector of distinct kmers. For each kmer, list the kmer itself and the number
        /// of haplotypes it appears in.
        std::vector<kmer_type> kmers;

        /// Number of haplotypes each kmer appears in.
        sdsl::int_vector<0> kmer_counts;

        /// Sequences as (GBWT sequence id, offset in the relevant node).
        std::vector<compact_sequence_type> sequences;

        // FIXME: compressed representation
        /// A bit vector marking the presence of kmers in the sequences.
        /// Sequence `i` contains kmer `j` if and only if `kmers_present[i * kmers.size() + j] == 1`.
        sdsl::bit_vector kmers_present;

        /// Calls `callback(kmer_id, is_present)` for each kmer in the given sequence,
        /// identified by its index in `sequences`, in the order the kmers are stored
        /// in `kmers`.
        void for_each_kmer(size_t sequence, const std::function<void(size_t, bool)>& callback) const;

        /// Calls `callback(kmer_id, num_present)` for each kmer in the given sequences,
        /// identified by their indices in `sequences`, in the order the kmers are stored
        /// in `kmers`.
        void for_each_kmer(size_t first, size_t second, const std::function<void(size_t, size_t)>& callback) const;

        /// Scores every haplotype in the subchain. The score of a haplotype is the sum of
        /// `kmer_score(kmer_id, is_present)` over all kmers, and it is reported by calling
        /// `haplotype_score(sequence, score)`, where `sequence` is the index in `sequences`.
        /// Iteration order and the number of `kmer_score` calls are unspecified.
        void score_haplotypes(
            const std::function<double(size_t, bool)>& kmer_score,
            const std::function<void(size_t, double)>& haplotype_score
        ) const;

        /// Returns the number of haplotypes in the subchain.
        size_t num_haplotypes() const { return this->sequences.size(); }

        /// Returns the number of kmers in the subchain.
        size_t num_kmers() const { return this->kmers.size(); }

        /// Returns the number of kmers present in the given haplotype.
        size_t num_present(size_t i) const;

        /// Returns the size of the kmer presence matrix.
        size_t total_kmers() const { return this->sequences.size() * this->kmers.size(); }

        /// Returns the total number of kmers present in the haplotypes.
        size_t total_present() const;

        /// Returns the start node as a GBWTGraph handle.
        handle_t start_handle() const { return gbwtgraph::GBWTGraph::node_to_handle(this->start); }

        /// Returns the end node as a GBWTGraph handle.
        handle_t end_handle() const { return gbwtgraph::GBWTGraph::node_to_handle(this->end); }

        /// Returns `true` if the subchain has a start node.
        bool has_start() const { return (this->type == normal || this->type == suffix); }

        /// Returns `true` if the subchain has an end node.
        bool has_end() const { return (this->type == normal || this->type == prefix); }

        /// Returns a string representation of the type and the boundary nodes.
        std::string to_string() const;

        /// Returns (sequence identifier, offset in a node) for the given sequence.
        sequence_type get_sequence(size_t i) const {
            return { this->sequences[i].first, this->sequences[i].second };
        }

        /// Returns the distance from the last base of `start` to the first base of
        /// `end` over the given sequence. Returns 0 if the subchain is not normal or
        /// if the sequence does not exist.
        size_t distance(const gbwtgraph::GBZ& gbz, size_t i) const;

        /// Returns an estimate of the badness of the subchain.
        /// The ideal value is 0.0, and higher values indicate worse subchains.
        /// The estimate is based on the following factors:
        /// * Length of the subchain.
        /// * Number of haplotypes relative to the expected number.
        /// * Information content of the kmers (disabled).
        double badness(const gbwtgraph::GBZ& gbz) const;

        /// Serializes the object to a stream in the Simple-SDS format.
        void simple_sds_serialize(std::ostream& out) const;

        /// Loads the object from a stream in the Simple-SDS format.
        void simple_sds_load(std::istream& in);

        /// Returns the size of the object in elements.
        size_t simple_sds_size() const;
    };

    /// Representation of a top-level chain.
    struct TopLevelChain {
        /// Offset in the child list of the root snarl.
        size_t offset;

        /// GBWT construction job for this chain.
        size_t job_id;

        /// Contig name corresponding to the chain.
        std::string contig_name;

        /// Subchains in the order they appear in.
        std::vector<Subchain> subchains;

        /// Serializes the object to a stream in the Simple-SDS format.
        void simple_sds_serialize(std::ostream& out) const;

        /// Loads the object from a stream in the Simple-SDS format.
        void simple_sds_load(std::istream& in);

        /// Returns the size of the object in elements.
        size_t simple_sds_size() const;
    };

    /// Returns the number of weakly connected components.
    size_t components() const { return this->header.top_level_chains; }

    /// Returns the number of GBWT construction jobs.
    size_t jobs() const { return this->header.construction_jobs; }

    /// Returns the length of the kmers.
    size_t k() const { return this->header.k; }

    /// Returns the number of kmers in the subchains.
    size_t kmers() const { return this->header.total_kmers; }

    Header header;

    // Tags for storing metadata, such as graph names.
    gbwt::Tags tags;

    // Job ids for each path in the GBWTGraph, or `jobs()` if the path is empty.
    std::vector<size_t> jobs_for_paths;

    std::vector<TopLevelChain> chains;

    /**
      * Returns a mapping from kmers to their counts in the given KFF file.
      * The counts include both the kmer and the reverse complement.
      *
      * Reads the KFF file using OpenMP threads. Exits with `std::exit()` if
      * the file cannot be opened and throws `std::runtime_error` if the kmer
      * counts cannot be used.
     */
    hash_map<Subchain::kmer_type, size_t> kmer_counts(const std::string& kff_file, Verbosity verbosity) const;

    /// Serializes the object to a stream in the Simple-SDS format.
    /// I/O errors can be detected by checking the stream state.
    void simple_sds_serialize(std::ostream& out) const;

    /// Serializes the object to a file in the Simple-SDS format.
    /// Prints an error message and exits the program on failure.
    void serialize_to(const std::string& filename) const;

    /// Loads the object from a stream in the Simple-SDS format.
    /// I/O errors can be detected by checking the stream state.
    /// Throws `sdsl::simple_sds::InvalidData` if data is unacceptable.
    void simple_sds_load(std::istream& in);

    /// Loads the object from a file in the Simple-SDS format.
    /// Prints an error message and exits the program on failure.
    void load_from(const std::string& filename);

    /// Returns the size of the object in elements.
    size_t simple_sds_size() const;

    /**
     * Assigns each reference and generic path in the graph to a GBWT construction job.
     *
     * For each path handle from 0 to gbz.named_paths() - 1, we assign the path to
     * the given construction job, or jobs() if the path is empty.
     */
    std::vector<size_t> assign_reference_paths(const gbwtgraph::GBZ& gbz, Verbosity verbosity) const;

    /**
     * Sets graph name (pggname) information based on the given GBZ graph.
     */
    void set_graph_name(const gbwtgraph::GBZ& gbz);

    /**
     * Returns the graph name (pggname) information stored in the tags.
     */
    gbwtgraph::GraphName graph_name() const { return gbwtgraph::GraphName(this->tags); }
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_RECOMBINATOR_HAPLOTYPES_HPP_INCLUDED
