#include "gf256.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Quarantined prototype.
//
// This file intentionally is not part of the mr-lrc-generator CMake target and
// has no public declaration. The construction below targets all-symbol MR-LRCs:
// global/heavy parity symbols are members of local groups. That is a different
// topology from the current data-local generator, whose local parity rows must
// only depend on data symbols in the same group.
//
// Keep this code unreachable until the project has an explicit all-symbol LRC
// model, checker, and output format.

namespace mrlrc::all_symbol_skew {
namespace {

struct LocalGroupLayout {
    std::vector<int> data;
    std::vector<int> global;
    int local_parity = 0;
    int local_row_start = 0;
};

struct Candidate {
    int data = 0;
    int local_rows = 0;
    int global_parity = 0;
    int symbols = 0;
    std::vector<LocalGroupLayout> groups;
    std::vector<uint8_t> matrix;
};

struct Result {
    bool supported = false;
    std::string message;
    Candidate candidate;
};

std::size_t idx(int row, int col, int cols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(cols) +
           static_cast<std::size_t>(col);
}

void reset_identity(Candidate *candidate)
{
    std::fill(candidate->matrix.begin(), candidate->matrix.end(), 0);
    for (int row = 0; row < candidate->data; row++) {
        candidate->matrix[idx(row, row, candidate->data)] = 1;
    }
}

std::vector<int> distinct_prime_factors(int value)
{
    std::vector<int> factors;
    for (int divisor = 2; divisor * divisor <= value; divisor++) {
        if (value % divisor != 0) {
            continue;
        }
        factors.push_back(divisor);
        while (value % divisor == 0) {
            value /= divisor;
        }
    }
    if (value > 1) {
        factors.push_back(value);
    }
    return factors;
}

std::vector<uint8_t> subfield_elements(int degree)
{
    int field_size = 1 << degree;
    std::vector<uint8_t> elements;
    elements.reserve(static_cast<std::size_t>(field_size));
    for (int value = 0; value < 256; value++) {
        uint8_t element = static_cast<uint8_t>(value);
        if (gf256_pow(element, static_cast<unsigned int>(field_size)) == element) {
            elements.push_back(element);
        }
    }
    return elements;
}

uint8_t multiplicative_generator(const std::vector<uint8_t> &field)
{
    int order = static_cast<int>(field.size()) - 1;
    if (order <= 1) {
        return 1;
    }

    std::vector<int> factors = distinct_prime_factors(order);
    for (uint8_t candidate : field) {
        if (candidate == 0) {
            continue;
        }

        bool is_generator = true;
        for (int factor : factors) {
            if (gf256_pow(candidate, static_cast<unsigned int>(order / factor)) == 1) {
                is_generator = false;
                break;
            }
        }
        if (is_generator) {
            return candidate;
        }
    }

    return 0;
}

bool span_contains(const std::vector<uint8_t> &span, uint8_t value)
{
    return std::find(span.begin(), span.end(), value) != span.end();
}

std::vector<uint8_t> add_to_span(const std::vector<uint8_t> &span,
                                 const std::vector<uint8_t> &base_field,
                                 uint8_t basis_element)
{
    std::vector<uint8_t> expanded = span;
    for (uint8_t existing : span) {
        for (uint8_t coefficient : base_field) {
            expanded.push_back(static_cast<uint8_t>(existing ^ gf256_mul(coefficient, basis_element)));
        }
    }

    std::sort(expanded.begin(), expanded.end());
    expanded.erase(std::unique(expanded.begin(), expanded.end()), expanded.end());
    return expanded;
}

std::vector<uint8_t> extension_basis(const std::vector<uint8_t> &base_field,
                                     const std::vector<uint8_t> &extension_field,
                                     int degree)
{
    std::vector<uint8_t> basis;
    std::vector<uint8_t> span(1, 0);

    for (uint8_t candidate : extension_field) {
        if (candidate == 0 || span_contains(span, candidate)) {
            continue;
        }
        basis.push_back(candidate);
        span = add_to_span(span, base_field, candidate);
        if (static_cast<int>(basis.size()) == degree) {
            return basis;
        }
    }

    return {};
}

int choose_base_field_degree(int minimum_q0, int extension_degree)
{
    const int candidate_degrees[] = {1, 2, 4, 8};
    for (int base_degree : candidate_degrees) {
        int q0 = 1 << base_degree;
        int total_degree = base_degree * extension_degree;
        if (q0 >= minimum_q0 && total_degree <= 8 && 8 % total_degree == 0) {
            return base_degree;
        }
    }
    return 0;
}

std::vector<uint8_t> multiply_matrices(const std::vector<uint8_t> &left,
                                       int left_rows,
                                       int shared,
                                       const std::vector<uint8_t> &right,
                                       int right_cols)
{
    std::vector<uint8_t> product(static_cast<std::size_t>(left_rows) *
                                     static_cast<std::size_t>(right_cols),
                                 0);
    for (int row = 0; row < left_rows; row++) {
        for (int col = 0; col < right_cols; col++) {
            uint8_t value = 0;
            for (int inner = 0; inner < shared; inner++) {
                value ^= gf256_mul(left[idx(row, inner, shared)], right[idx(inner, col, right_cols)]);
            }
            product[idx(row, col, right_cols)] = value;
        }
    }
    return product;
}

Result build_candidate(const std::vector<std::vector<int>> &data_groups,
                       int data,
                       int local_parity,
                       int global_parity)
{
    Result result;
    int groups = static_cast<int>(data_groups.size());
    if (groups <= 0) {
        result.message = "all-symbol skew construction requires at least one group";
        return result;
    }
    if (global_parity <= 0) {
        result.message = "all-symbol skew construction requires global parity";
        return result;
    }
    if ((data + global_parity) % groups != 0) {
        result.message = "all-symbol skew construction requires equal paper groups";
        return result;
    }

    int payload_per_group = (data + global_parity) / groups;
    int paper_group_size = payload_per_group + local_parity;
    int skew_degree = std::min(payload_per_group, global_parity);
    int minimum_q0 = std::max(groups + 1, paper_group_size);
    int base_degree = choose_base_field_degree(minimum_q0, skew_degree);
    if (base_degree == 0) {
        result.message = "all-symbol skew construction does not fit GF(256)";
        return result;
    }

    int q0 = 1 << base_degree;
    int extension_field_degree = base_degree * skew_degree;
    std::vector<uint8_t> base_field = subfield_elements(base_degree);
    std::vector<uint8_t> extension_field = subfield_elements(extension_field_degree);
    std::vector<uint8_t> basis = extension_basis(base_field, extension_field, skew_degree);
    uint8_t gamma = multiplicative_generator(extension_field);
    if (static_cast<int>(basis.size()) != skew_degree || gamma == 0 ||
        static_cast<int>(base_field.size()) < paper_group_size) {
        result.message = "all-symbol skew construction failed to initialize field data";
        return result;
    }

    int local_rows = groups * local_parity;
    int total_parity = local_rows + global_parity;
    int symbols = data + total_parity;

    Candidate candidate;
    candidate.data = data;
    candidate.local_rows = local_rows;
    candidate.global_parity = global_parity;
    candidate.symbols = symbols;
    candidate.groups.resize(static_cast<std::size_t>(groups));
    candidate.matrix.assign(static_cast<std::size_t>(symbols) * static_cast<std::size_t>(data), 0);

    std::vector<uint8_t> alpha(static_cast<std::size_t>(paper_group_size));
    for (int col = 0; col < paper_group_size; col++) {
        alpha[static_cast<std::size_t>(col)] = base_field[static_cast<std::size_t>(col)];
    }

    std::vector<uint8_t> beta(static_cast<std::size_t>(paper_group_size), 0);
    for (int col = 0; col < paper_group_size; col++) {
        uint8_t value = 0;
        uint8_t alpha_value = alpha[static_cast<std::size_t>(col)];
        for (int basis_index = 0; basis_index < skew_degree; basis_index++) {
            uint8_t coefficient =
                gf256_pow(alpha_value, static_cast<unsigned int>(local_parity + basis_index));
            value ^= gf256_mul(coefficient, basis[static_cast<std::size_t>(basis_index)]);
        }
        beta[static_cast<std::size_t>(col)] = value;
    }

    std::vector<std::vector<int>> group_symbols(static_cast<std::size_t>(groups));
    int next_global = 0;
    for (int group_index = 0; group_index < groups; group_index++) {
        LocalGroupLayout group;
        group.data = data_groups[static_cast<std::size_t>(group_index)];
        group.local_parity = local_parity;
        group.local_row_start = data + group_index * local_parity;

        std::vector<int> symbols = group.data;
        int heavy_in_group = payload_per_group - static_cast<int>(group.data.size());
        for (int heavy = 0; heavy < heavy_in_group; heavy++) {
            int symbol = data + local_rows + next_global;
            symbols.push_back(symbol);
            group.global.push_back(symbol);
            next_global++;
        }
        for (int local = 0; local < local_parity; local++) {
            symbols.push_back(group.local_row_start + local);
        }

        candidate.groups[static_cast<std::size_t>(group_index)] = std::move(group);
        group_symbols[static_cast<std::size_t>(group_index)] = std::move(symbols);
    }

    std::vector<uint8_t> parity_check(static_cast<std::size_t>(total_parity) *
                                          static_cast<std::size_t>(symbols),
                                      0);
    for (int group_index = 0; group_index < groups; group_index++) {
        const auto &symbols_in_group = group_symbols[static_cast<std::size_t>(group_index)];
        for (int local = 0; local < local_parity; local++) {
            int row = group_index * local_parity + local;
            for (int col = 0; col < paper_group_size; col++) {
                parity_check[idx(row, symbols_in_group[static_cast<std::size_t>(col)], symbols)] =
                    gf256_pow(alpha[static_cast<std::size_t>(col)], static_cast<unsigned int>(local));
            }
        }
    }

    int first_global_check = groups * local_parity;
    for (int group_index = 0; group_index < groups; group_index++) {
        const auto &symbols_in_group = group_symbols[static_cast<std::size_t>(group_index)];
        uint8_t conjugacy_term = gf256_pow(gamma, static_cast<unsigned int>(group_index));
        uint8_t conjugacy_factor = 1;
        std::vector<uint8_t> beta_power = beta;

        for (int global = 0; global < global_parity; global++) {
            int row = first_global_check + global;
            for (int col = 0; col < paper_group_size; col++) {
                parity_check[idx(row, symbols_in_group[static_cast<std::size_t>(col)], symbols)] =
                    gf256_mul(conjugacy_factor, beta_power[static_cast<std::size_t>(col)]);
            }

            conjugacy_factor = gf256_mul(conjugacy_factor, conjugacy_term);
            conjugacy_term = gf256_pow(conjugacy_term, static_cast<unsigned int>(q0));
            for (auto &value : beta_power) {
                value = gf256_pow(value, static_cast<unsigned int>(q0));
            }
        }
    }

    std::vector<uint8_t> parity_block(static_cast<std::size_t>(total_parity) *
                                          static_cast<std::size_t>(total_parity),
                                      0);
    for (int row = 0; row < total_parity; row++) {
        for (int col = 0; col < total_parity; col++) {
            parity_block[idx(row, col, total_parity)] =
                parity_check[idx(row, data + col, symbols)];
        }
    }

    std::vector<uint8_t> parity_inverse;
    if (!gf256_invert_matrix(parity_block, &parity_inverse, total_parity)) {
        result.message = "all-symbol skew construction has singular parity block";
        return result;
    }

    std::vector<uint8_t> data_block(static_cast<std::size_t>(total_parity) *
                                        static_cast<std::size_t>(data),
                                    0);
    for (int row = 0; row < total_parity; row++) {
        for (int col = 0; col < data; col++) {
            data_block[idx(row, col, data)] = parity_check[idx(row, col, symbols)];
        }
    }

    std::vector<uint8_t> parity_rows =
        multiply_matrices(parity_inverse, total_parity, total_parity, data_block, data);

    reset_identity(&candidate);
    for (int parity = 0; parity < total_parity; parity++) {
        int row = data + parity;
        for (int col = 0; col < data; col++) {
            candidate.matrix[idx(row, col, data)] = parity_rows[idx(parity, col, data)];
        }
    }

    result.supported = true;
    result.message = "all-symbol skew prototype built over GF(" +
                     std::to_string(1 << extension_field_degree) + ")";
    result.candidate = std::move(candidate);
    return result;
}

} // namespace
} // namespace mrlrc::all_symbol_skew
