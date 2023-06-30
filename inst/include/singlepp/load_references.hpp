#ifndef SINGLEPP_LOADERS_HPP
#define SINGLEPP_LOADERS_HPP

#include "macros.hpp"

#include "byteme/PerByte.hpp"
#include "byteme/RawFileReader.hpp"
#include "tatami/tatami.hpp"

#ifdef SINGLEPP_USE_ZLIB
#include "byteme/GzipFileReader.hpp"
#include "byteme/ZlibBufferReader.hpp"
#endif

#include "Markers.hpp"

#include <string>
#include <vector>
#include <cctype>
#include <type_traits>

/**
 * @file load_references.hpp
 *
 * @brief Load reference datasets from a few expected formats.
 */

namespace singlepp {

/** 
 * @cond
 */
template<bool parallel_>
using SuperPerByte = typename std::conditional<parallel_, byteme::PerByte<char>, byteme::PerByteParallel<char> >::type;

template<bool parallel_, class Reader_>
std::vector<int> load_labels_internal(Reader_& reader) {
    bool non_empty = false;
    int current = 0;
    std::vector<int> labels;
    bool remaining = true;

    SuperPerByte<parallel_> pb(&reader);
    bool okay = pb.valid();
    while (okay) {
        char x = pb.get();
        okay = pb.advance();

        if (x == '\n') {
            if (!non_empty) {
                throw std::runtime_error("label index must be an integer");
            }
            labels.push_back(current);
            current = 0;
            non_empty = false;
        } else if (std::isdigit(x)) {
            non_empty = true;
            current *= 10;
            current += (x - '0');
        } else {
            throw std::runtime_error("label index must be an integer");
        }
    } 

    if (non_empty) {
        labels.push_back(current);
    }

    return labels;
}
/** 
 * @endcond
 */

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param path Path to a text file containing the labels.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return Vector containing the label index for each reference profile.
 *
 * The file should contain one line per profile, containing an integer label index for that profile.
 * Label indices refer to another array containing the actual names of the labels.
 * The total number of lines should be equal to the number of profiles in the dataset.
 * The file should not contain any header.
 */
template<bool parallel_ = false>
std::vector<int> load_labels_from_text_file(const char* path, size_t buffer_size = 65536) {
    byteme::RawFileReader reader(path, buffer_size);
    return load_labels_internal<parallel_>(reader);
}

#ifdef SINGLEPP_USE_ZLIB

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param path Path to a Gzip-compressed file containing the labels.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return Vector containing the label index for each reference profile.
 *
 * See `load_labels_from_text_file()` for details about the format.
 */
template<bool parallel_ = false>
std::vector<int> load_labels_from_gzip_file(const char* path, size_t buffer_size = 65536) {
    byteme::GzipFileReader reader(path, buffer_size);
    return load_labels_internal<parallel_>(reader);
}

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param[in] buffer Pointer to an array containing a Zlib/Gzip-compressed string of labels.
 * @param len Length of the array for `buffer`.
 * @param buffer_size Size of the buffer to use when decompressing the buffer.
 *
 * @return Vector containing the label index for each reference profile.
 *
 * See `load_labels_from_text_file()` for details about the format.
 */
template<bool parallel_ = false>
inline std::vector<int> load_labels_from_zlib_buffer(const unsigned char* buffer, size_t len, size_t buffer_size = 65536) {
    byteme::ZlibBufferReader reader(buffer, len, 3, buffer_size);
    return load_labels_internal<parallel_>(reader);
}

#endif

/**
 * @cond
 */
template<bool parallel_, class Reader_>
std::vector<std::string> load_names_internal(Reader_& reader) {
    std::string current;
    std::vector<std::string> names;

    SuperPerByte<parallel_> pb(&reader);
    bool okay = pb.valid();
    while (okay) {
        char x = pb.get();
        okay = pb.advance();

        if (x == '\n') {
            names.emplace_back(std::move(current));
            current.clear();
        } else {
            current += x;
        }
    }

    if (!current.empty()) { // absence of trailing newline is okay.
        names.emplace_back(std::move(current));
    }

    return names;
};
/** 
 * @endcond
 */

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param path Path to a text file containing the label names.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return Vector of strings containing the name for each reference profile.
 *
 * The file should contain one line per label, containing a (non-quoted) string with the full name for that label.
 * The total number of lines should be equal to the number of unique labels in the dataset.
 * The file should not contain any header.
 */
template<bool parallel_ = false>
std::vector<std::string> load_label_names_from_text_file(const char* path, size_t buffer_size = 65536) {
    byteme::RawFileReader reader(path, buffer_size);
    return load_names_internal<parallel_>(reader);
}

#ifdef SINGLEPP_USE_ZLIB

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param path Path to a Gzip-compressed file containing the label names.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return Vector of strings containing the name for each label.
 *
 * See `load_label_names_from_text_file()` for details about the format.
 */
template<bool parallel_ = false>
std::vector<std::string> load_label_names_from_gzip_file(const char* path, size_t buffer_size = 65536) {
    byteme::GzipFileReader reader(path, buffer_size);
    return load_names_internal<parallel_>(reader);
}

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param[in] buffer Pointer to an array containing a Zlib/Gzip-compressed string of label names.
 * @param len Length of the array for `buffer`.
 * @param buffer_size Size of the buffer to use when decompressing the buffer.
 *
 * @return Vector of strings containing the name for each label.
 *
 * See `load_label_names_from_text_file()` for details about the format.
 */
template<bool parallel_ = false>
std::vector<std::string> load_label_names_from_zlib_buffer(const unsigned char* buffer, size_t len, size_t buffer_size = 65536) {
    byteme::ZlibBufferReader reader(buffer, len, 3, buffer_size);
    return load_names_internal<parallel_>(reader);
}

#endif

/** 
 * @cond
 */
template<bool parallel_, class Reader_>
std::pair<std::vector<std::string>, std::vector<std::string> > load_features_internal(Reader_& reader) {
    std::string current;
    std::vector<std::string> ensembl, symbols;

    SuperPerByte<parallel_> pb(&reader);
    bool okay = pb.valid();
    while (okay) {
        current.clear();

        // Pulling out the Ensembl ID.
        do {
            char x = pb.get();
            if (x == ',') { // don't advance yet, so that okay check below doesn't trigger if the symbol is empty and file is not newline terminated.
                break;
            } else if (x == '\n') {
                okay = false; // hit the error below.
                break;
            }
            current += x;
            okay = pb.advance();
        } while (okay);

        if (!okay) {
            throw std::runtime_error("two comma-separated fields (Ensembl ID and symbol) expected on each line");
        }
        okay = pb.advance();

        ensembl.push_back(std::move(current));
        current.clear();

        // Now pulling out the gene symbol.
        while (okay) {
            char x = pb.get();
            okay = pb.advance(); 
            if (x == '\n') {
                break;
            } 
            current += x;
        }
   
        symbols.push_back(std::move(current)); // still gets added if file is not newline terminated.
        current.clear();
    }

    return std::make_pair(std::move(ensembl), std::move(symbols));
}
/** 
 * @endcond
 */

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param path Path to a text file containing the feature annotation.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return Pair of vectors, each of length equal to the number of features.
 * The first contains Ensembl IDs while the second contains gene symbols.
 *
 * The file should contain one line per feature, with total number of lines equal to the number of features in the dataset.
 * Each line should contain two strings separated by a comma.
 * The first string should be the Ensembl ID while the second string should be the gene symbol; either string may be empty.
 * The file should not contain any header.
 */
template<bool parallel_ = false>
std::pair<std::vector<std::string>, std::vector<std::string> > load_features_from_text_file(const char* path, size_t buffer_size = 65536) {
    byteme::RawFileReader reader(path, buffer_size);
    return load_features_internal<parallel_>(reader);
}

#ifdef SINGLEPP_USE_ZLIB

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param path Path to a Gzip-compressed file containing the feature annotation.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return Pair of vectors, each of length equal to the number of features.
 * The first contains Ensembl IDs while the second contains gene symbols.
 *
 * See `load_features_from_text_file()` for details about the format.
 */
template<bool parallel_ = false>
std::pair<std::vector<std::string>, std::vector<std::string> > load_features_from_gzip_file(const char* path, size_t buffer_size = 65536) {
    byteme::GzipFileReader reader(path, buffer_size);
    return load_features_internal<parallel_>(reader);
}

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param[in] buffer Pointer to an array containing a Zlib/Gzip-compressed string containing the feature annotation.
 * @param len Length of the array for `buffer`.
 * @param buffer_size Size of the buffer to use when decompressing the buffer.
 *
 * @return Pair of vectors, each of length equal to the number of features.
 * The first contains Ensembl IDs while the second contains gene symbols.
 *
 * See `load_features_from_text_file()` for details about the format.
 */
template<bool parallel_ = false>
std::pair<std::vector<std::string>, std::vector<std::string> > load_features_from_zlib_buffer(const unsigned char* buffer, size_t len, size_t buffer_size = 65536) {
    byteme::ZlibBufferReader reader(buffer, len, 3, buffer_size);
    return load_features_internal<parallel_>(reader);
}

#endif

/**
 * Matrix of ranks as a dense column-major matrix.
 * Each column corresponds to a sample while each row corresponds to a feature.
 * Each column contains the ranked expression values for all features.
 *
 * @tparam Data Numeric type for data in the matrix interface.
 * @tparam Index Integer type for indices in the matrix interface.
 */
template<typename Data = double, typename Index = int>
using RankMatrix = tatami::DenseColumnMatrix<Data, Index, std::vector<int> >;

/** 
 * @cond
 */
template<typename Data, typename Index, bool parallel_, class Reader>
RankMatrix<Data, Index> load_rankings_internal(Reader& reader) {
    size_t nfeatures = 0;
    size_t line = 0;
    std::vector<int> values;

    int field = 0;
    bool non_empty = false;
    int current = 0;

    bool has_nfeatures = false;
    auto check_nfeatures = [&]() -> void {
        if (!has_nfeatures) {
            has_nfeatures = true;
            nfeatures = field + 1;
        } else if (field + 1 != nfeatures) {
            throw std::runtime_error("number of fields on each line should be equal to the number of features");
        }
    };

    SuperPerByte<parallel_> pb(&reader);
    bool okay = pb.valid();
    while (okay) {
        char x = pb.get();
        okay = pb.advance();

        if (x == '\n') {
            check_nfeatures();
            if (!non_empty) {
                throw std::runtime_error("fields should not be empty");
            }
            values.push_back(current);
            current = 0;
            field = 0;
            non_empty = false;
            ++line;

        } else if (x == ',') {
            if (!non_empty) {
                throw std::runtime_error("fields should not be empty");
            }
            values.push_back(current);
            current = 0;
            ++field;
            non_empty = false;

        } else if (std::isdigit(x)) {
            non_empty = true;
            current *= 10;
            current += (x - '0');

        } else {
            throw std::runtime_error("fields should only contain integer ranks");
        }
    }

    if (field || non_empty) { // aka no terminating newline.
        check_nfeatures();
        if (!non_empty) {
            throw std::runtime_error("fields should not be empty");
        }
        values.push_back(current);
        ++line;
    }

    return RankMatrix<Data, Index>(nfeatures, line, std::move(values));
};
/** 
 * @endcond
 */

/**
 * @tparam Data Numeric type for data in the matrix interface.
 * @tparam Index Integer type for indices in the matrix interface.
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 * 
 * @param path Path to a text file containing the ranking matrix.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return A `RankMatrix` containing the feature rankings for each reference profile.
 * Each column corresponds to a reference profile while each row corresponds to a feature.
 *
 * The file should contain one line per reference profile, with the total number of lines equal to the number of profiles in the dataset.
 * Each line should contain the rank of each feature's expression within that profile, separated by commas.
 * The number of comma-separated fields on each line should be equal to the number of features.
 * Ranks should be strictly integer - tied ranks should default to the minimum rank among the index set of ties.
 */
template<typename Data = double, typename Index = int, bool parallel_ = false>
RankMatrix<Data, Index> load_rankings_from_text_file(const char* path, size_t buffer_size = 65536) {
    byteme::RawFileReader reader(path, buffer_size);
    return load_rankings_internal<Data, Index, parallel_>(reader);
}

#ifdef SINGLEPP_USE_ZLIB

/**
 * @tparam Data Numeric type for data in the matrix interface.
 * @tparam Index Integer type for indices in the matrix interface.
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param path Path to a Gzip-compressed file containing the ranking matrix.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return A `RankMatrix` containing the feature rankings for each reference profile.
 * Each column corresponds to a reference profile while each row corresponds to a feature.
 *
 * See `load_rankings_from_text_file()` for details about the format.
 */
template<typename Data = double, typename Index = int, bool parallel_ = false>
RankMatrix<Data, Index> load_rankings_from_gzip_file(const char* path, size_t buffer_size = 65536) {
    byteme::GzipFileReader reader(path, buffer_size);
    return load_rankings_internal<Data, Index, parallel_>(reader);
}

/**
 * @tparam Data Numeric type for data in the matrix interface.
 * @tparam Index Integer type for indices in the matrix interface.
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param[in] buffer Pointer to an array containing a Zlib/Gzip-compressed string containing the ranking matrix.
 * @param len Length of the array for `buffer`.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return A `RankMatrix` containing the feature rankings for each reference profile.
 * Each column corresponds to a reference profile while each row corresponds to a feature.
 *
 * See `load_rankings_from_text_file()` for details about the format.
 */
template<typename Data = double, typename Index = int, bool parallel_ = false>
RankMatrix<Data, Index> load_rankings_from_zlib_buffer(const unsigned char* buffer, size_t len, size_t buffer_size = 65536) {
    byteme::ZlibBufferReader reader(buffer, len, 3, buffer_size);
    return load_rankings_internal<Data, Index, parallel_>(reader);
}

#endif

/** 
 * @cond
 */
template<bool parallel_, class Reader>
Markers load_markers_internal(size_t nfeatures, size_t nlabels, Reader& reader) {
    Markers markers(nlabels);
    for (auto& m : markers) {
        m.resize(nlabels);
    }

    std::vector<int> values;
    SuperPerByte<parallel_> pb(&reader);
    bool okay = pb.valid();
    while (okay) {

        // Processing the label IDs.
        size_t first = 0, second = 0;
        for (int l = 0; l < 2; ++l) {
            auto& current = (l == 0 ? first : second);
            bool non_empty = false;

            do {
                char x = pb.get();
                okay = pb.advance();

                if (x == '\t') {
                    if (!non_empty) {
                        throw std::runtime_error("empty field detected in the label indices");
                    }
                    break;
                } else if (x == '\n') {
                    okay = false; // hit the error below.
                    break;
                } else if (!std::isdigit(x)) {
                    throw std::runtime_error("label indices should be integers");
                }

                non_empty = true;
                current *= 10;
                current += (x - '0');
            } while (okay);

            if (!okay) {
                throw std::runtime_error("expected at least three tab-separated fields on each line");
            }
            if (current >= markers.size()) {
                throw std::runtime_error("label index out of range");
            }
        }

        // Processing the actual gene indices.
        bool non_empty = false;
        int current = 0;
        while (okay) {
            char x = pb.get();
            okay = pb.advance();

            if (std::isdigit(x)) {
                non_empty = true;
                current *= 10;
                current += (x - '0');

            } else if (x == '\t') {
                if (!non_empty) {
                    throw std::runtime_error("gene index fields should not be empty");
                }
                values.push_back(current);
                current = 0;
                non_empty = false;

            } else if (x == '\n') {
                break;

            } else {
                throw std::runtime_error("gene index fields should be integers");
            }
        }

        // Adding the last element. We don't do this inside the newline check,
        // as we need to account for cases where the file is not newline-terminated.
        if (!non_empty) {
            throw std::runtime_error("gene index fields should not be empty");
        }
        values.push_back(current);

        for (auto v : values) {
            if (static_cast<size_t>(v) >= nfeatures) {
                throw std::runtime_error("gene index out of range");
            }
        }

        auto& store = markers[first][second];
        if (!store.empty()) {
            throw std::runtime_error("multiple marker sets listed for a single pairwise comparison");
        }
        store.swap(values); // implicit clear of 'values', as store is empty.
    }

    return markers;
}
/** 
 * @endcond
 */

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param path Path to a text file containing the marker lists.
 * @param nfeatures Total number of features in the dataset.
 * @param nlabels Number of labels in the dataset.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return A `Markers` object containing the markers from each pairwise comparison between labels.
 *
 * The file should contain one line per pairwise comparison between labels.
 * Each line should at least 3 tab-delimited fields - the index of the first label, the index of the second label, 
 * and then the indices of the features selected as marker genes for the first label relative to the second.
 * Any (non-zero) number of marker indices may be reported provided they are ordered by marker strength.
 * The total number of lines in this file should be equal to the total number of pairwise comparisons between different labels, including permutations.
 */
template<bool parallel_ = false>
Markers load_markers_from_text_file(const char* path, size_t nfeatures, size_t nlabels, size_t buffer_size = 65536) {
    byteme::RawFileReader reader(path, buffer_size);
    return load_markers_internal<parallel_>(nfeatures, nlabels, reader);
}

#ifdef SINGLEPP_USE_ZLIB

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param path Path to a Gzip-compressed file containing the marker lists.
 * @param nfeatures Total number of features in the dataset.
 * @param nlabels Number of labels in the dataset.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return A `Markers` object containing the markers from each pairwise comparison between labels.
 *
 * See `load_markers_from_text_file()` for details about the format.
 */
template<bool parallel_ = false>
Markers load_markers_from_gzip_file(const char* path, size_t nfeatures, size_t nlabels, size_t buffer_size = 65536) {
    byteme::GzipFileReader reader(path, buffer_size);
    return load_markers_internal<parallel_>(nfeatures, nlabels, reader);
}

/**
 * @tparam parallel_ Whether file loading and parsing should be parallelized.
 *
 * @param[in] buffer Pointer to an array containing a Zlib/Gzip-compressed string containing the marker lists.
 * @param len Length of the array for `buffer`.
 * @param nfeatures Total number of features in the dataset.
 * @param nlabels Number of labels in the dataset.
 * @param buffer_size Size of the buffer to use when reading the file.
 *
 * @return A `Markers` object containing the markers from each pairwise comparison between labels.
 *
 * See `load_markers_from_text_file()` for details about the format.
 */
template<bool parallel_ = false>
Markers load_markers_from_zlib_buffer(const unsigned char* buffer, size_t len, size_t nfeatures, size_t nlabels, size_t buffer_size = 65536) {
    byteme::ZlibBufferReader reader(buffer, len, 3, buffer_size);
    return load_markers_internal<parallel_>(nfeatures, nlabels, reader);
}

#endif

}

#endif
