#include "recombinator_haplotypes.hpp"

#include "kff.hpp"

#include <cmath>

//#define debug_selected_haplotypes

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t Haplotypes::SUBCHAIN_LENGTH;
constexpr size_t Haplotypes::KFF_BLOCK_SIZE;

constexpr std::uint32_t Haplotypes::Header::MAGIC_NUMBER;
constexpr std::uint32_t Haplotypes::Header::VERSION;
constexpr std::uint32_t Haplotypes::Header::VERSION_WITH_TAGS;
constexpr std::uint32_t Haplotypes::Header::MIN_VERSION;
constexpr std::uint64_t Haplotypes::Header::DEFAULT_K;

//------------------------------------------------------------------------------

hash_map<Haplotypes::Subchain::kmer_type, size_t>::iterator
find_kmer(hash_map<Haplotypes::Subchain::kmer_type, size_t>& counts, Haplotypes::Subchain::kmer_type kmer, size_t k) {
    Haplotypes::Subchain::kmer_type rc = minimizer_reverse_complement(kmer, k);
    auto forward = counts.find(kmer);
    auto reverse = counts.find(rc);
    return (forward != counts.end() ? forward : reverse);
}

hash_map<Haplotypes::Subchain::kmer_type, size_t> Haplotypes::kmer_counts(const std::string& kff_file, Verbosity verbosity) const {
    double start = gbwt::readTimer();
    if (verbosity >= verbosity_basic) {
        std::cerr << "Reading kmer counts" << std::endl;
    }

    // Open and validate the kmer count file.
    ParallelKFFReader reader(kff_file);

    // Populate the map with the kmers we are interested in.
    double checkpoint = gbwt::readTimer();
    hash_map<Subchain::kmer_type, size_t> result;
    result.reserve(this->header.total_kmers);
    for (size_t chain_id = 0; chain_id < this->chains.size(); chain_id++) {
        const TopLevelChain& chain = this->chains[chain_id];
        for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
            const Subchain& subchain = chain.subchains[subchain_id];
            for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                result[subchain.kmers[kmer_id]] = 0;
            }
        }
    }
    if (verbosity >= verbosity_detailed) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Initialized the hash map with " << result.size() << " kmers in " << seconds << " seconds" << std::endl;
    }

    // Read the KFF file and add the counts using multiple threads.
    checkpoint = gbwt::readTimer();
    size_t kmer_count = 0;
    #pragma omp parallel
    {
        #pragma omp task
        {
            while (true) {
                std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> block = reader.read(KFF_BLOCK_SIZE);
                if (block.empty()) {
                    break;
                }
                std::vector<std::pair<hash_map<Subchain::kmer_type, size_t>::iterator, size_t>> buffer;
                for (auto kmer : block) {
                    auto iter = find_kmer(result, kmer.first, this->k());
                    if (iter != result.end()) {
                        buffer.push_back({ iter, kmer.second });
                    }
                }
                #pragma omp critical
                {
                    for (auto to_update : buffer) {
                        to_update.first->second += to_update.second;
                    }
                    kmer_count += block.size();
                }
            }
        }
    }
    if (verbosity >= verbosity_detailed) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Read " << kmer_count << " kmers in " << seconds << " seconds" << std::endl;
    }

    if (verbosity >= verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Read the kmer counts in " << seconds << " seconds" << std::endl;
    }
    return result;
}

//------------------------------------------------------------------------------

std::string Haplotypes::Subchain::to_string() const {
    std::string result;
    switch (this->type) {
    case normal:
        result.append("normal");
        break;
    case prefix:
        result.append("prefix");
        break;
    case suffix:
        result.append("suffix");
        break;
    case full_haplotype:
        result.append("full");
        break;
    default:
        result.append("invalid");
        break;
    }

    result.append(" from ");
    result.append(to_string_gbwtgraph(this->start));
    result.append(" to ");
    result.append(to_string_gbwtgraph(this->end));

    return result;
}

size_t Haplotypes::Subchain::distance(const gbwtgraph::GBZ& gbz, size_t i) const {
    if (this->type != normal || i >= this->sequences.size()) {
        return 0;
    }

    size_t result = 1;
    gbwt::edge_type curr(this->start, this->sequences[i].second);
    while (true) {
        curr = gbz.index.LF(curr);
        if (curr.first == gbwt::ENDMARKER || curr.first == this->end) {
            break;
        }
        result += gbz.graph.get_length(gbwtgraph::GBWTGraph::node_to_handle(curr.first));
    }

    return result;
}

// TODO: What is the right formula? Should subchain length be a parameter?
double Haplotypes::Subchain::badness(const gbwtgraph::GBZ& gbz) const {
    double result = 0.0;

    // Factor 1: Subchain length, ideally over a reference path.
    if (this->type == normal) {
        size_t selected = 0;
        for (size_t i = 0; i < this->sequences.size(); i++) {
            gbwt::size_type path_id = gbwt::Path::id(this->sequences[i].first);
            auto found = gbz.graph.id_to_path.find(path_id);
            if (found != gbz.graph.id_to_path.end()) {
                selected = i; break;
            }
        }
        size_t length = this->distance(gbz, selected);
        if (length > SUBCHAIN_LENGTH) {
            result += std::log(static_cast<double>(length) / static_cast<double>(SUBCHAIN_LENGTH));
        }
    }

    // Factor 2: Number of haplotypes relative to the expected number.
    size_t expected_haplotypes = gbz.index.metadata.haplotypes();
    size_t haplotypes = this->sequences.size();
    if (haplotypes < expected_haplotypes) {
        result += std::log(static_cast<double>(expected_haplotypes) / static_cast<double>(haplotypes));
    }

    // Factor 3: Information content of the kmers.
    // Disabled for the moment.
    /* double expected_entropy = 4.0 * std::log(static_cast<double>(haplotypes));
    double entropy = 0.0;
    for (size_t i = 0; i < this->kmer_counts.size(); i++) {
        double p = static_cast<double>(this->kmer_counts[i]) / static_cast<double>(haplotypes);
        entropy -= p * std::log(p);
    }
    if (entropy < expected_entropy) {
        result += expected_entropy - entropy;
    } */

    return result;
}

void Haplotypes::Subchain::simple_sds_serialize(std::ostream& out) const {
    sdsl::simple_sds::serialize_value<std::uint64_t>(this->type, out);
    sdsl::simple_sds::serialize_value<gbwt::node_type>(this->start, out);
    sdsl::simple_sds::serialize_value<gbwt::node_type>(this->end, out);
    sdsl::simple_sds::serialize_vector(this->kmers, out);
    this->kmer_counts.simple_sds_serialize(out);
    sdsl::simple_sds::serialize_vector(this->sequences, out);
    this->kmers_present.simple_sds_serialize(out);
}

void load_subchain_header(Haplotypes::Subchain& subchain, std::istream& in) {
    std::uint64_t temp = sdsl::simple_sds::load_value<std::uint64_t>(in);
    switch (temp) {
    case Haplotypes::Subchain::normal: // Fall through.
    case Haplotypes::Subchain::prefix: // Fall through.
    case Haplotypes::Subchain::suffix: // Fall through.
    case Haplotypes::Subchain::full_haplotype:
        subchain.type = static_cast<Haplotypes::Subchain::subchain_t>(temp);
        break;
    default:
        throw sdsl::simple_sds::InvalidData("Invalid subchain type: " + std::to_string(temp));
    }

    subchain.start = sdsl::simple_sds::load_value<gbwt::node_type>(in);
    subchain.end = sdsl::simple_sds::load_value<gbwt::node_type>(in);
    bool should_have_start = (subchain.type == Haplotypes::Subchain::normal || subchain.type == Haplotypes::Subchain::suffix);
    bool should_have_end = (subchain.type == Haplotypes::Subchain::normal || subchain.type == Haplotypes::Subchain::prefix);
    if ((subchain.start != gbwt::ENDMARKER) != should_have_start) {
        throw sdsl::simple_sds::InvalidData("Subchain start node " + std::to_string(subchain.start) + " does not match type " + std::to_string(temp));
    }
    if ((subchain.end != gbwt::ENDMARKER) != should_have_end) {
        throw sdsl::simple_sds::InvalidData("Subchain end node " + std::to_string(subchain.end) + " does not match type" + std::to_string(temp));
    }
}

void load_subchain_kmers_present(Haplotypes::Subchain& subchain, std::istream& in) {
    subchain.kmers_present.simple_sds_load(in);
    if (subchain.kmers_present.size() != subchain.kmers.size() * subchain.sequences.size()) {
        throw sdsl::simple_sds::InvalidData("Invalid length for the kmer presence bitvector in subchain from " +
            std::to_string(subchain.start) + " to " + std::to_string(subchain.end));
    }
}

void Haplotypes::Subchain::simple_sds_load(std::istream& in) {
    load_subchain_header(*this, in);

    this->kmers = sdsl::simple_sds::load_vector<kmer_type>(in);
    this->kmer_counts.simple_sds_load(in);
    this->sequences = sdsl::simple_sds::load_vector<compact_sequence_type>(in);

    load_subchain_kmers_present(*this, in);
}

size_t Haplotypes::Subchain::simple_sds_size() const {
    size_t result = sdsl::simple_sds::value_size<std::uint64_t>() + 2 * sdsl::simple_sds::value_size<gbwt::node_type>();
    result += sdsl::simple_sds::vector_size(this->kmers);
    result += this->kmer_counts.simple_sds_size();
    result += sdsl::simple_sds::vector_size(this->sequences);
    result += this->kmers_present.simple_sds_size();
    return result;
}

void Haplotypes::TopLevelChain::simple_sds_serialize(std::ostream& out) const {
    sdsl::simple_sds::serialize_value<size_t>(this->offset, out);
    sdsl::simple_sds::serialize_value<size_t>(this->job_id, out);
    sdsl::simple_sds::serialize_string(this->contig_name, out);
    sdsl::simple_sds::serialize_value<size_t>(this->subchains.size(), out);
    for (auto& subchain : this->subchains) {
        subchain.simple_sds_serialize(out);
    }
}

void Haplotypes::TopLevelChain::simple_sds_load(std::istream& in) {
    this->offset = sdsl::simple_sds::load_value<size_t>(in);
    this->job_id = sdsl::simple_sds::load_value<size_t>(in);
    this->contig_name = sdsl::simple_sds::load_string(in);
    size_t subchain_count = sdsl::simple_sds::load_value<size_t>(in);
    this->subchains.resize(subchain_count);
    for (size_t i = 0; i < subchain_count; i++) {
        this->subchains[i].simple_sds_load(in);
    }
}

size_t Haplotypes::TopLevelChain::simple_sds_size() const {
    size_t result = 3 * sdsl::simple_sds::value_size<size_t>();
    result += sdsl::simple_sds::string_size(this->contig_name);
    for (auto& subchain : this->subchains) {
        result += subchain.simple_sds_size();
    }
    return result;
}

void Haplotypes::simple_sds_serialize(std::ostream& out) const {
    sdsl::simple_sds::serialize_value<Header>(this->header, out);
    this->tags.simple_sds_serialize(out);
    sdsl::simple_sds::serialize_vector(this->jobs_for_paths, out);
    for (auto& chain : this->chains) {
        chain.simple_sds_serialize(out);
    }
}

void Haplotypes::serialize_to(const std::string& filename) const {
    try {
        sdsl::simple_sds::serialize_to(*this, filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [Haplotypes] Serialization to " << filename << " failed: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void Haplotypes::simple_sds_load(std::istream& in) {
    this->header = sdsl::simple_sds::load_value<Header>(in);
    if (this->header.magic_number != Header::MAGIC_NUMBER) {
        throw sdsl::simple_sds::InvalidData("Haplotypes::simple_sds_load(): Expected magic number " + std::to_string(Header::MAGIC_NUMBER) +
            ", got " + std::to_string(this->header.magic_number));
    }
    if (this->header.version < Header::MIN_VERSION || this->header.version > Header::VERSION) {
        std::string msg = "Haplotypes::simple_sds_load(): Expected version " + std::to_string(Header::MIN_VERSION)
            + " to " + std::to_string(Header::VERSION) + ", got version " + std::to_string(this->header.version);
        throw sdsl::simple_sds::InvalidData(msg);
    }

    // Tags are present from version 6 onwards.
    if (this->header.version >= Header::VERSION_WITH_TAGS) {
        this->tags.simple_sds_load(in);
    }

    this->jobs_for_paths = sdsl::simple_sds::load_vector<size_t>(in);

    this->chains.resize(this->header.top_level_chains);
    for (auto& chain : this->chains) {
        chain.simple_sds_load(in);
    }

    // Update to the current version.
    this->header.version = Header::VERSION;
}

void Haplotypes::load_from(const std::string& filename) {
    try {
        sdsl::simple_sds::load_from(*this, filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [Haplotypes] Loading from " << filename << " failed: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

size_t Haplotypes::simple_sds_size() const {
    size_t result = sdsl::simple_sds::value_size<Header>();
    result += sdsl::simple_sds::vector_size(this->jobs_for_paths);
    for (auto& chain : this->chains) {
        result += chain.simple_sds_size();
    }
    return result;
}

std::vector<size_t> Haplotypes::assign_reference_paths(const gbwtgraph::GBZ& gbz, Verbosity verbosity) const {
    if (verbosity >= verbosity_basic) {
        std::cerr << "Assigning reference paths to GBWT construction jobs" << std::endl;
    }
    double start = gbwt::readTimer();

    if (this->jobs_for_paths.size() != gbz.index.metadata.paths()) {
        std::string msg = "Haplotypes::assign_reference_paths(): Haplotype information was built for "
            + std::to_string(this->jobs_for_paths.size()) + " paths, but the graph contains "
            + std::to_string(gbz.index.metadata.paths()) + " paths";
        throw std::runtime_error(msg);
    }

    // All paths are initially unassigned.
    std::vector<size_t> result (gbz.named_paths(), this->jobs());
    size_t unassigned = 0;
    for (size_t i = 0; i < gbz.named_paths(); i++) {
        gbwt::size_type path_id = gbz.graph.named_paths[i].id;
        result[i] = this->jobs_for_paths[path_id];
        if (result[i] >= this->jobs()) {
            unassigned++;
        }
    }

    if (verbosity >= verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Assigned " << (result.size() - unassigned) << " reference paths (" << unassigned << " unassigned) in " << seconds << " seconds" << std::endl;
    }
    return result;
}

void Haplotypes::set_graph_name(const gbwtgraph::GBZ& gbz) {
    gbwtgraph::GraphName name = gbz.graph_name();
    name.set_tags(this->tags);
}

//------------------------------------------------------------------------------

} // namespace vg
