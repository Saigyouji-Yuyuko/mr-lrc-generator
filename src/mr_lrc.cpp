#include "mr_lrc.h"

#include "gf256.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mrlrc {
namespace {

constexpr uint64_t kMaxThreadCount = 256;

class ValidationError : public std::runtime_error {
public:
    explicit ValidationError(const std::string &message) : std::runtime_error(message) {}
};

uint64_t splitmix64_value(uint64_t value)
{
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

uint64_t seed_for_attempt(uint64_t seed, uint64_t attempt)
{
    if (attempt == 1) {
        return seed;
    }
    return splitmix64_value(seed ^ (attempt * UINT64_C(0x9e3779b97f4a7c15)));
}

struct Rng {
    uint64_t state;

    explicit Rng(uint64_t seed) : state(seed) {}

    uint64_t next()
    {
        uint64_t z = (state += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27U)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31U);
    }

    uint64_t range(uint64_t upper)
    {
        if (upper <= 1) {
            return 0;
        }

        uint64_t threshold = static_cast<uint64_t>(-upper) % upper;
        for (;;) {
            uint64_t x = next();
            if (x >= threshold) {
                return x % upper;
            }
        }
    }

    uint8_t nonzero_byte()
    {
        return static_cast<uint8_t>(1U + range(255));
    }

    uint8_t byte()
    {
        return static_cast<uint8_t>(range(256));
    }
};

std::size_t idx(int row, int col, int cols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(cols) +
           static_cast<std::size_t>(col);
}

template <typename T>
void shuffle(Rng &rng, std::vector<T> *items)
{
    for (std::size_t i = items->size(); i > 1; i--) {
        std::size_t j = static_cast<std::size_t>(rng.range(i));
        std::swap((*items)[i - 1], (*items)[j]);
    }
}

std::vector<uint8_t> cauchy_matrix(int rows, int cols, Rng &rng)
{
    if (rows < 0 || cols < 0 || rows + cols > 256) {
        throw ValidationError("cauchy matrix requires rows + cols <= 256 over GF(256)");
    }

    std::vector<uint8_t> matrix(static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols), 0);
    if (rows == 0 || cols == 0) {
        return matrix;
    }

    std::vector<uint8_t> pool(256);
    for (int i = 0; i < 256; i++) {
        pool[static_cast<std::size_t>(i)] = static_cast<uint8_t>(i);
    }
    shuffle(rng, &pool);

    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            matrix[idx(row, col, cols)] =
                gf256_inv(static_cast<uint8_t>(pool[static_cast<std::size_t>(row)] ^
                                               pool[static_cast<std::size_t>(rows + col)]));
        }
    }

    return matrix;
}

std::vector<uint8_t> vandermonde_matrix(int rows, int cols, Rng &rng)
{
    if (rows < 0 || cols < 0 || rows > 255 || cols > 255) {
        throw ValidationError("vandermonde matrix requires rows <= 255 and cols <= 255 over GF(256)");
    }

    std::vector<uint8_t> matrix(static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols), 0);
    if (rows == 0 || cols == 0) {
        return matrix;
    }

    std::vector<uint8_t> points(255);
    for (int i = 0; i < 255; i++) {
        points[static_cast<std::size_t>(i)] = static_cast<uint8_t>(i + 1);
    }
    shuffle(rng, &points);

    std::vector<uint8_t> column_scale(static_cast<std::size_t>(cols));
    for (int col = 0; col < cols; col++) {
        column_scale[static_cast<std::size_t>(col)] = rng.nonzero_byte();
    }

    for (int row = 0; row < rows; row++) {
        uint8_t x = points[static_cast<std::size_t>(row)];
        for (int col = 0; col < cols; col++) {
            uint8_t v = gf256_pow(x, static_cast<unsigned int>(col));
            matrix[idx(row, col, cols)] = gf256_mul(column_scale[static_cast<std::size_t>(col)], v);
        }
    }

    return matrix;
}

std::vector<uint8_t> random_matrix(int rows, int cols, Rng &rng)
{
    if (rows < 0 || cols < 0) {
        throw ValidationError("random matrix dimensions must be non-negative");
    }

    std::vector<uint8_t> matrix(static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols), 0);
    for (auto &value : matrix) {
        value = rng.byte();
    }
    return matrix;
}

std::vector<uint8_t> build_dense_matrix(MatrixFamily family, int rows, int cols, Rng &rng)
{
    switch (family) {
    case MatrixFamily::Cauchy:
        return cauchy_matrix(rows, cols, rng);
    case MatrixFamily::Vandermonde:
        return vandermonde_matrix(rows, cols, rng);
    case MatrixFamily::Random:
        return random_matrix(rows, cols, rng);
    }

    throw ValidationError("unknown matrix construction method");
}

std::vector<LocalGroup> build_groups(const Params &params)
{
    if (params.groups <= 0) {
        throw ValidationError("--groups must be positive");
    }
    if (params.groups > params.data) {
        throw ValidationError("--groups cannot exceed --data");
    }
    if (params.local_parity <= 0) {
        throw ValidationError("--local-parity must be positive");
    }

    std::vector<LocalGroup> groups;
    groups.reserve(static_cast<std::size_t>(params.groups));

    int base = params.data / params.groups;
    int extra = params.data % params.groups;
    int start = 0;
    int local_row_start = params.data;
    for (int group_index = 0; group_index < params.groups; group_index++) {
        int size = base + (group_index < extra ? 1 : 0);
        if (params.local_family == MatrixFamily::Cauchy && size + params.local_parity > 256) {
            throw ValidationError("each local group needs data_count + local_parity <= 256");
        }
        if (params.local_family == MatrixFamily::Vandermonde &&
            (size > 255 || params.local_parity > 255)) {
            throw ValidationError("vandermonde local rows need group_data <= 255 and local_parity <= 255");
        }

        LocalGroup group;
        group.local_parity = params.local_parity;
        group.local_row_start = local_row_start;
        group.data.reserve(static_cast<std::size_t>(size));
        for (int i = 0; i < size; i++) {
            group.data.push_back(start + i);
        }

        groups.push_back(std::move(group));
        start += size;
        local_row_start += params.local_parity;
    }

    return groups;
}

Code make_empty_code(const Params &params)
{
    if (params.data <= 0) {
        throw ValidationError("--data must be positive");
    }
    if (params.global_parity < 0) {
        throw ValidationError("--global-parity must be non-negative");
    }
    if (params.random_limit == 0) {
        throw ValidationError("--random-limit must be positive");
    }
    if (params.thread_count == 0) {
        throw ValidationError("--thread-count must be positive");
    }
    if (params.thread_count > kMaxThreadCount) {
        throw ValidationError("--thread-count cannot exceed 256");
    }

    Code code;
    code.data = params.data;
    code.groups = build_groups(params);
    code.local_rows = params.groups * params.local_parity;
    code.global_parity = params.global_parity;
    code.total_parity = code.local_rows + code.global_parity;
    code.seed = params.seed;
    code.local_family = params.local_family;
    code.global_family = params.global_family;
    code.construction = params.construction;

    if (code.global_family == MatrixFamily::Cauchy && code.global_parity > 0 &&
        code.data + code.global_parity > 256) {
        throw ValidationError("cauchy global rows need data + global_parity <= 256 over GF(256)");
    }
    if (code.global_family == MatrixFamily::Vandermonde) {
        if (code.data > 255) {
            throw ValidationError("vandermonde global rows need data <= 255 over GF(256)");
        }
        if (code.global_parity > 255) {
            throw ValidationError("vandermonde global rows need global_parity <= 255 over GF(256)");
        }
    }

    code.symbols = code.data + code.total_parity;
    code.matrix.assign(static_cast<std::size_t>(code.symbols) * static_cast<std::size_t>(code.data), 0);
    return code;
}

void reset_identity(Code *code)
{
    std::fill(code->matrix.begin(), code->matrix.end(), 0);
    for (int row = 0; row < code->data; row++) {
        code->matrix[idx(row, row, code->data)] = 1;
    }
}

void build_local_rows(Code *code, Rng &rng)
{
    for (const auto &group : code->groups) {
        int group_data = static_cast<int>(group.data.size());
        std::vector<uint8_t> local =
            build_dense_matrix(code->local_family, group.local_parity, group_data, rng);

        for (int local_row = 0; local_row < group.local_parity; local_row++) {
            int matrix_row = group.local_row_start + local_row;
            for (int col = 0; col < group_data; col++) {
                int data_col = group.data[static_cast<std::size_t>(col)];
                code->matrix[idx(matrix_row, data_col, code->data)] =
                    local[idx(local_row, col, group_data)];
            }
        }
    }
}

void build_global_rows(Code *code, Rng &rng)
{
    int first_global = code->data + code->local_rows;
    std::vector<uint8_t> global = build_dense_matrix(code->global_family, code->global_parity, code->data, rng);
    for (int row = 0; row < code->global_parity; row++) {
        for (int col = 0; col < code->data; col++) {
            code->matrix[idx(first_global + row, col, code->data)] = global[idx(row, col, code->data)];
        }
    }
}

void build_candidate(Code *code, Rng &rng)
{
    reset_identity(code);
    build_local_rows(code, rng);
    build_global_rows(code, rng);
}

std::vector<int> group_symbol_ids(const Code &code, int group_index)
{
    const auto &group = code.groups[static_cast<std::size_t>(group_index)];
    std::vector<int> symbols = group.data;
    for (int local = 0; local < group.local_parity; local++) {
        symbols.push_back(group.local_row_start + local);
    }
    return symbols;
}

struct RrefResult {
    int rank = 0;
    std::vector<int> pivot_columns;
    std::vector<uint8_t> matrix;
};

RrefResult gf256_rref(std::vector<uint8_t> matrix, int rows, int cols)
{
    RrefResult result;
    result.matrix = std::move(matrix);
    if (rows <= 0 || cols <= 0) {
        return result;
    }

    int rank = 0;
    for (int col = 0; col < cols && rank < rows; col++) {
        int pivot = -1;
        for (int row = rank; row < rows; row++) {
            if (result.matrix[idx(row, col, cols)] != 0) {
                pivot = row;
                break;
            }
        }
        if (pivot < 0) {
            continue;
        }

        if (pivot != rank) {
            for (int c = 0; c < cols; c++) {
                std::swap(result.matrix[idx(rank, c, cols)], result.matrix[idx(pivot, c, cols)]);
            }
        }

        uint8_t inv_pivot = gf256_inv(result.matrix[idx(rank, col, cols)]);
        for (int c = 0; c < cols; c++) {
            result.matrix[idx(rank, c, cols)] = gf256_mul(result.matrix[idx(rank, c, cols)], inv_pivot);
        }

        for (int row = 0; row < rows; row++) {
            if (row == rank) {
                continue;
            }
            uint8_t factor = result.matrix[idx(row, col, cols)];
            if (factor == 0) {
                continue;
            }
            for (int c = 0; c < cols; c++) {
                result.matrix[idx(row, c, cols)] ^=
                    gf256_mul(factor, result.matrix[idx(rank, c, cols)]);
            }
        }

        result.pivot_columns.push_back(col);
        rank++;
    }

    result.rank = rank;
    return result;
}

int gf256_rank(const std::vector<uint8_t> &matrix, int rows, int cols)
{
    return gf256_rref(matrix, rows, cols).rank;
}

std::vector<uint8_t> gf256_nullspace_basis(const std::vector<uint8_t> &matrix, int rows, int cols)
{
    RrefResult rref = gf256_rref(matrix, rows, cols);
    std::vector<uint8_t> is_pivot(static_cast<std::size_t>(cols), 0);
    for (int pivot_col : rref.pivot_columns) {
        is_pivot[static_cast<std::size_t>(pivot_col)] = 1;
    }

    std::vector<int> free_columns;
    for (int col = 0; col < cols; col++) {
        if (is_pivot[static_cast<std::size_t>(col)] == 0) {
            free_columns.push_back(col);
        }
    }

    int nullity = static_cast<int>(free_columns.size());
    std::vector<uint8_t> basis(static_cast<std::size_t>(cols) * static_cast<std::size_t>(nullity), 0);
    for (int basis_col = 0; basis_col < nullity; basis_col++) {
        int free_col = free_columns[static_cast<std::size_t>(basis_col)];
        basis[idx(free_col, basis_col, nullity)] = 1;
        for (int pivot_row = 0; pivot_row < rref.rank; pivot_row++) {
            int pivot_col = rref.pivot_columns[static_cast<std::size_t>(pivot_row)];
            basis[idx(pivot_col, basis_col, nullity)] = rref.matrix[idx(pivot_row, free_col, cols)];
        }
    }
    return basis;
}

std::vector<uint8_t> transpose_matrix(const std::vector<uint8_t> &matrix, int rows, int cols)
{
    std::vector<uint8_t> transposed(static_cast<std::size_t>(cols) * static_cast<std::size_t>(rows), 0);
    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            transposed[idx(col, row, rows)] = matrix[idx(row, col, cols)];
        }
    }
    return transposed;
}

std::vector<uint8_t> column_space_key(const std::vector<uint8_t> &matrix, int rows, int cols)
{
    if (cols == 0) {
        return {};
    }
    std::vector<uint8_t> transposed = transpose_matrix(matrix, rows, cols);
    RrefResult rref = gf256_rref(std::move(transposed), cols, rows);
    std::vector<uint8_t> key;
    key.reserve(static_cast<std::size_t>(rref.rank) * static_cast<std::size_t>(rows));
    for (int row = 0; row < rref.rank; row++) {
        for (int col = 0; col < rows; col++) {
            key.push_back(rref.matrix[idx(row, col, rows)]);
        }
    }
    return key;
}

std::vector<uint8_t> append_columns(const std::vector<uint8_t> &left, int rows, int left_cols,
                                    const std::vector<uint8_t> &right, int right_cols)
{
    std::vector<uint8_t> joined(static_cast<std::size_t>(rows) *
                                    static_cast<std::size_t>(left_cols + right_cols),
                                0);
    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < left_cols; col++) {
            joined[idx(row, col, left_cols + right_cols)] = left[idx(row, col, left_cols)];
        }
        for (int col = 0; col < right_cols; col++) {
            joined[idx(row, left_cols + col, left_cols + right_cols)] = right[idx(row, col, right_cols)];
        }
    }
    return joined;
}

std::vector<uint8_t> project_rows(const std::vector<uint8_t> &matrix, int rows, int cols,
                                  const std::vector<int> &kept_rows)
{
    (void)rows;
    std::vector<uint8_t> projected(static_cast<std::size_t>(kept_rows.size()) *
                                       static_cast<std::size_t>(cols),
                                   0);
    for (int row = 0; row < static_cast<int>(kept_rows.size()); row++) {
        int source_row = kept_rows[static_cast<std::size_t>(row)];
        for (int col = 0; col < cols; col++) {
            projected[idx(row, col, cols)] = matrix[idx(source_row, col, cols)];
        }
    }
    return projected;
}

template <typename Fn>
bool enumerate_combinations(const std::vector<int> &items, int need, int start, std::vector<int> *chosen,
                            Fn &&fn)
{
    if (static_cast<int>(chosen->size()) == need) {
        return fn(*chosen);
    }

    int remaining = need - static_cast<int>(chosen->size());
    for (int i = start; i <= static_cast<int>(items.size()) - remaining; i++) {
        chosen->push_back(items[static_cast<std::size_t>(i)]);
        if (!enumerate_combinations(items, need, i + 1, chosen, fn)) {
            chosen->pop_back();
            return false;
        }
        chosen->pop_back();
    }
    return true;
}

struct ResidualBlock {
    int extra = 0;
    std::vector<int> erased_symbols;
    std::vector<uint8_t> matrix;
};

struct GroupResidualBlocks {
    std::vector<std::vector<ResidualBlock>> by_extra;
};

struct DualParent {
    bool reachable = false;
    int previous_sum = 0;
    int extra = 0;
    const ResidualBlock *block = nullptr;
};

struct DualGroupOptions {
    std::vector<const ResidualBlock *> by_extra;
};

class ResidualVerifier {
public:
    explicit ResidualVerifier(const Code &code, CheckResult *result)
        : code_(code), result_(result), group_erased_(code.groups.size()), h_(code.global_parity)
    {
    }

    void run()
    {
        precompute_residual_blocks();
        if (!should_stop()) {
            if (h_ == 2) {
                verify_h2_fast_path();
            } else if (h_ == 3 || h_ == 4) {
                verify_by_dual_index();
            } else {
                std::vector<uint8_t> empty;
                verify_from_group(0, 0, empty);
            }
        }
        result_->is_mr = result_->failures == 0;
        result_->message =
            result_->is_mr
                ? "candidate is recoverable for every maximal residual erasure pattern"
                : "candidate is singular for at least one maximal residual erasure pattern";
    }

private:
    bool should_stop() const
    {
        return result_->failures > 0;
    }

    std::vector<uint8_t> build_local_submatrix(int group_index, const std::vector<int> &symbols) const
    {
        const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
        int local_rows = group.local_parity;
        std::vector<uint8_t> matrix(static_cast<std::size_t>(local_rows) *
                                        static_cast<std::size_t>(symbols.size()),
                                    0);
        for (int col = 0; col < static_cast<int>(symbols.size()); col++) {
            int symbol = symbols[static_cast<std::size_t>(col)];
            if (symbol < code_.data) {
                for (int local = 0; local < local_rows; local++) {
                    matrix[idx(local, col, static_cast<int>(symbols.size()))] =
                        code_.matrix[idx(group.local_row_start + local, symbol, code_.data)];
                }
                continue;
            }

            int local_symbol = symbol - group.local_row_start;
            if (local_symbol >= 0 && local_symbol < local_rows) {
                matrix[idx(local_symbol, col, static_cast<int>(symbols.size()))] = 1;
            }
        }
        return matrix;
    }

    std::vector<uint8_t> build_global_submatrix(const std::vector<int> &symbols) const
    {
        std::vector<uint8_t> matrix(static_cast<std::size_t>(h_) * static_cast<std::size_t>(symbols.size()),
                                    0);
        int first_global = code_.data + code_.local_rows;
        for (int col = 0; col < static_cast<int>(symbols.size()); col++) {
            int symbol = symbols[static_cast<std::size_t>(col)];
            if (symbol >= code_.data) {
                continue;
            }
            for (int global = 0; global < h_; global++) {
                matrix[idx(global, col, static_cast<int>(symbols.size()))] =
                    code_.matrix[idx(first_global + global, symbol, code_.data)];
            }
        }
        return matrix;
    }

    std::vector<uint8_t> build_residual_block(int group_index, const std::vector<int> &erased_symbols,
                                              int extra) const
    {
        const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
        std::vector<uint8_t> local = build_local_submatrix(group_index, erased_symbols);
        std::vector<uint8_t> nullspace =
            gf256_nullspace_basis(local, group.local_parity, static_cast<int>(erased_symbols.size()));
        std::vector<uint8_t> global = build_global_submatrix(erased_symbols);
        std::vector<uint8_t> residual(static_cast<std::size_t>(h_) * static_cast<std::size_t>(extra), 0);

        for (int row = 0; row < h_; row++) {
            for (int col = 0; col < extra; col++) {
                uint8_t value = 0;
                for (int inner = 0; inner < static_cast<int>(erased_symbols.size()); inner++) {
                    value ^= gf256_mul(global[idx(row, inner, static_cast<int>(erased_symbols.size()))],
                                       nullspace[idx(inner, col, extra)]);
                }
                residual[idx(row, col, extra)] = value;
            }
        }
        return residual;
    }

    void mark_failure_with_erased(const std::vector<std::vector<int>> &group_erased,
                                  const std::vector<int> &kept_global_rows)
    {
        std::vector<uint8_t> erased(static_cast<std::size_t>(code_.symbols), 0);
        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            if (!group_erased[static_cast<std::size_t>(group_index)].empty()) {
                for (int symbol : group_erased[static_cast<std::size_t>(group_index)]) {
                    erased[static_cast<std::size_t>(symbol)] = 1;
                }
                continue;
            }

            std::vector<int> symbols = group_symbol_ids(code_, group_index);
            const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
            for (int i = 0; i < group.local_parity; i++) {
                erased[static_cast<std::size_t>(symbols[static_cast<std::size_t>(i)])] = 1;
            }
        }

        std::vector<uint8_t> keep_global(static_cast<std::size_t>(h_), 0);
        for (int row : kept_global_rows) {
            keep_global[static_cast<std::size_t>(row)] = 1;
        }
        int first_global = code_.data + code_.local_rows;
        for (int global = 0; global < h_; global++) {
            if (keep_global[static_cast<std::size_t>(global)] == 0) {
                erased[static_cast<std::size_t>(first_global + global)] = 1;
            }
        }

        result_->failures = 1;
        result_->first_failed_erased = erased_symbols_from_mask(code_, erased);
    }

    std::vector<int> first_kept_global_rows(int count) const
    {
        std::vector<int> kept;
        kept.reserve(static_cast<std::size_t>(count));
        for (int row = 0; row < count && row < h_; row++) {
            kept.push_back(row);
        }
        return kept;
    }

    void mark_precompute_failure(int group_index, const std::vector<int> &erased_symbols, int extra)
    {
        std::vector<std::vector<int>> group_erased(code_.groups.size());
        group_erased[static_cast<std::size_t>(group_index)] = erased_symbols;
        mark_failure_with_erased(group_erased, first_kept_global_rows(extra));
    }

    void precompute_group(int group_index)
    {
        const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
        std::vector<int> symbols = group_symbol_ids(code_, group_index);
        int max_extra = std::min(h_, static_cast<int>(symbols.size()) - group.local_parity);
        GroupResidualBlocks group_blocks;
        group_blocks.by_extra.resize(static_cast<std::size_t>(h_ + 1));

        for (int extra = 0; extra <= max_extra; extra++) {
            int need = group.local_parity + extra;
            std::vector<int> chosen;
            enumerate_combinations(symbols, need, 0, &chosen, [&](const std::vector<int> &erased_symbols) {
                if (should_stop()) {
                    return false;
                }

                std::vector<uint8_t> local = build_local_submatrix(group_index, erased_symbols);
                result_->patterns_checked++;
                if (gf256_rank(local, group.local_parity, need) < group.local_parity) {
                    mark_precompute_failure(group_index, erased_symbols, extra);
                    return false;
                }

                if (extra == 0) {
                    return true;
                }

                ResidualBlock block;
                block.extra = extra;
                block.erased_symbols = erased_symbols;
                block.matrix = build_residual_block(group_index, erased_symbols, extra);
                result_->patterns_checked++;
                if (gf256_rank(block.matrix, h_, extra) < extra) {
                    mark_precompute_failure(group_index, erased_symbols, extra);
                    return false;
                }
                group_blocks.by_extra[static_cast<std::size_t>(extra)].push_back(std::move(block));
                return true;
            });
            if (should_stop()) {
                return;
            }
        }

        for (int extra = 1; extra <= max_extra; extra++) {
            auto &blocks = group_blocks.by_extra[static_cast<std::size_t>(extra)];
            std::set<std::vector<uint8_t>> seen_spaces;
            std::vector<ResidualBlock> unique_blocks;
            unique_blocks.reserve(blocks.size());
            for (auto &block : blocks) {
                std::vector<uint8_t> key = column_space_key(block.matrix, h_, extra);
                if (seen_spaces.insert(std::move(key)).second) {
                    unique_blocks.push_back(std::move(block));
                }
            }
            blocks = std::move(unique_blocks);
        }

        residual_blocks_.push_back(std::move(group_blocks));
    }

    void precompute_residual_blocks()
    {
        residual_blocks_.clear();
        residual_blocks_.reserve(code_.groups.size());
        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            precompute_group(group_index);
            if (should_stop()) {
                return;
            }
        }
    }

    std::string state_key(int group_index, int used, const std::vector<uint8_t> &current) const
    {
        std::string key = std::to_string(group_index) + ":" + std::to_string(used) + ":";
        std::vector<uint8_t> space = column_space_key(current, h_, used);
        if (!space.empty()) {
            key.append(reinterpret_cast<const char *>(space.data()), space.size());
        }
        return key;
    }

    void check_global_completions(const std::vector<uint8_t> &current, int used, int start,
                                  std::vector<int> *kept_rows)
    {
        if (should_stop()) {
            return;
        }
        if (static_cast<int>(kept_rows->size()) == used) {
            result_->patterns_checked++;
            if (used == 0) {
                return;
            }

            std::vector<uint8_t> projected = project_rows(current, h_, used, *kept_rows);
            if (gf256_rank(projected, used, used) < used) {
                mark_failure_with_erased(group_erased_, *kept_rows);
            }
            return;
        }

        int remaining = used - static_cast<int>(kept_rows->size());
        for (int row = start; row <= h_ - remaining; row++) {
            kept_rows->push_back(row);
            check_global_completions(current, used, row + 1, kept_rows);
            kept_rows->pop_back();
            if (should_stop()) {
                return;
            }
        }
    }

    bool verify_global_completions(const std::vector<uint8_t> &current, int used)
    {
        std::vector<int> kept_rows;
        check_global_completions(current, used, 0, &kept_rows);
        return !should_stop();
    }

    std::vector<uint8_t> block_vector(const ResidualBlock &block) const
    {
        std::vector<uint8_t> vector(static_cast<std::size_t>(h_), 0);
        for (int row = 0; row < h_; row++) {
            vector[static_cast<std::size_t>(row)] = block.matrix[idx(row, 0, block.extra)];
        }
        return vector;
    }

    bool verify_h2_fast_path()
    {
        std::unordered_map<std::string, std::pair<int, const ResidualBlock *>> seen_lines;
        for (int group_index = 0; group_index < static_cast<int>(residual_blocks_.size()); group_index++) {
            const auto &group_blocks = residual_blocks_[static_cast<std::size_t>(group_index)];
            if (group_blocks.by_extra.size() <= 1) {
                continue;
            }

            for (const auto &block : group_blocks.by_extra[1]) {
                result_->patterns_checked++;
                std::vector<uint8_t> vector = block_vector(block);
                for (int row = 0; row < h_; row++) {
                    if (vector[static_cast<std::size_t>(row)] == 0) {
                        std::vector<std::vector<int>> group_erased(code_.groups.size());
                        group_erased[static_cast<std::size_t>(group_index)] = block.erased_symbols;
                        mark_failure_with_erased(group_erased, std::vector<int>{row});
                        return false;
                    }
                }

                std::string line = canonical_dual_key(vector);
                auto found = seen_lines.find(line);
                if (found != seen_lines.end()) {
                    int previous_group = found->second.first;
                    if (previous_group != group_index) {
                        std::vector<std::vector<int>> group_erased(code_.groups.size());
                        group_erased[static_cast<std::size_t>(previous_group)] =
                            found->second.second->erased_symbols;
                        group_erased[static_cast<std::size_t>(group_index)] = block.erased_symbols;
                        mark_failure_with_erased(group_erased, std::vector<int>{0, 1});
                        return false;
                    }
                } else {
                    seen_lines.emplace(std::move(line), std::make_pair(group_index, &block));
                }
            }
        }

        return true;
    }

    int support_size(const std::vector<uint8_t> &y) const
    {
        int support = 0;
        for (uint8_t value : y) {
            if (value != 0) {
                support++;
            }
        }
        return support;
    }

    std::vector<int> kept_rows_for_dual_failure(const std::vector<uint8_t> &y, int used) const
    {
        std::vector<int> kept;
        kept.reserve(static_cast<std::size_t>(used));
        for (int row = 0; row < h_; row++) {
            if (y[static_cast<std::size_t>(row)] != 0) {
                kept.push_back(row);
            }
        }
        for (int row = 0; row < h_ && static_cast<int>(kept.size()) < used; row++) {
            if (y[static_cast<std::size_t>(row)] == 0) {
                kept.push_back(row);
            }
        }
        return kept;
    }

    void mark_dual_failure(const std::vector<std::vector<DualParent>> &parents, int stage, int used,
                           const std::vector<uint8_t> &y)
    {
        std::vector<std::vector<int>> group_erased(code_.groups.size());
        int sum = used;
        for (int current_stage = stage; current_stage > 0; current_stage--) {
            const DualParent &parent = parents[static_cast<std::size_t>(current_stage)]
                                              [static_cast<std::size_t>(sum)];
            if (parent.extra > 0 && parent.block != nullptr) {
                group_erased[static_cast<std::size_t>(current_stage - 1)] =
                    parent.block->erased_symbols;
            }
            sum = parent.previous_sum;
        }

        mark_failure_with_erased(group_erased, kept_rows_for_dual_failure(y, used));
    }

    bool check_dual_options(const std::vector<uint8_t> &y, const std::vector<DualGroupOptions> &options)
    {
        result_->patterns_checked++;
        int min_used = support_size(y);
        std::vector<uint8_t> reachable(static_cast<std::size_t>(h_ + 1), 0);
        reachable[0] = 1;

        std::vector<std::vector<DualParent>> parents(
            code_.groups.size() + 1,
            std::vector<DualParent>(static_cast<std::size_t>(h_ + 1)));
        parents[0][0].reachable = true;

        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            const auto &option = options[static_cast<std::size_t>(group_index)].by_extra;

            std::vector<uint8_t> next(static_cast<std::size_t>(h_ + 1), 0);
            for (int sum = 0; sum <= h_; sum++) {
                if (reachable[static_cast<std::size_t>(sum)] == 0) {
                    continue;
                }

                if (next[static_cast<std::size_t>(sum)] == 0) {
                    next[static_cast<std::size_t>(sum)] = 1;
                    DualParent &parent = parents[static_cast<std::size_t>(group_index + 1)]
                                                [static_cast<std::size_t>(sum)];
                    parent.reachable = true;
                    parent.previous_sum = sum;
                    parent.extra = 0;
                    parent.block = nullptr;
                }

                for (int extra = 1; extra <= h_ - sum; extra++) {
                    const ResidualBlock *block = option[static_cast<std::size_t>(extra)];
                    if (block == nullptr) {
                        continue;
                    }
                    int next_sum = sum + extra;
                    if (next[static_cast<std::size_t>(next_sum)] != 0) {
                        continue;
                    }
                    next[static_cast<std::size_t>(next_sum)] = 1;
                    DualParent &parent = parents[static_cast<std::size_t>(group_index + 1)]
                                                [static_cast<std::size_t>(next_sum)];
                    parent.reachable = true;
                    parent.previous_sum = sum;
                    parent.extra = extra;
                    parent.block = block;
                }
            }

            for (int used = min_used; used <= h_; used++) {
                if (next[static_cast<std::size_t>(used)] != 0) {
                    mark_dual_failure(parents, group_index + 1, used, y);
                    return false;
                }
            }

            reachable = std::move(next);
        }

        return true;
    }

    std::string canonical_dual_key(const std::vector<uint8_t> &y) const
    {
        std::string key;
        key.resize(static_cast<std::size_t>(h_));
        int pivot = -1;
        for (int row = 0; row < h_; row++) {
            if (y[static_cast<std::size_t>(row)] != 0) {
                pivot = row;
                break;
            }
        }
        if (pivot < 0) {
            return {};
        }

        uint8_t inv_pivot = gf256_inv(y[static_cast<std::size_t>(pivot)]);
        for (int row = 0; row < h_; row++) {
            key[static_cast<std::size_t>(row)] =
                static_cast<char>(gf256_mul(y[static_cast<std::size_t>(row)], inv_pivot));
        }
        return key;
    }

    std::vector<uint8_t> left_nullspace_basis(const ResidualBlock &block) const
    {
        std::vector<uint8_t> transposed(static_cast<std::size_t>(block.extra) *
                                            static_cast<std::size_t>(h_),
                                        0);
        for (int row = 0; row < h_; row++) {
            for (int col = 0; col < block.extra; col++) {
                transposed[idx(col, row, h_)] = block.matrix[idx(row, col, block.extra)];
            }
        }
        return gf256_nullspace_basis(transposed, block.extra, h_);
    }

    std::vector<uint8_t> dual_from_coefficients(const std::vector<uint8_t> &basis, int dim,
                                                const std::vector<uint8_t> &coefficients) const
    {
        std::vector<uint8_t> y(static_cast<std::size_t>(h_), 0);
        for (int row = 0; row < h_; row++) {
            uint8_t value = 0;
            for (int col = 0; col < dim; col++) {
                value ^= gf256_mul(basis[idx(row, col, dim)],
                                   coefficients[static_cast<std::size_t>(col)]);
            }
            y[static_cast<std::size_t>(row)] = value;
        }
        return y;
    }

    std::vector<uint8_t> dual_vector_from_key(const std::string &key) const
    {
        std::vector<uint8_t> y(static_cast<std::size_t>(h_), 0);
        for (int row = 0; row < h_; row++) {
            y[static_cast<std::size_t>(row)] = static_cast<uint8_t>(key[static_cast<std::size_t>(row)]);
        }
        return y;
    }

    void add_dual_index_option(int group_index, const ResidualBlock &block, const std::vector<uint8_t> &basis,
                               int dim, const std::vector<uint8_t> &coefficients)
    {
        std::vector<uint8_t> y = dual_from_coefficients(basis, dim, coefficients);
        std::string key = canonical_dual_key(y);
        if (key.empty()) {
            return;
        }

        auto inserted = dual_index_.emplace(std::move(key), std::vector<DualGroupOptions>{});
        std::vector<DualGroupOptions> &options = inserted.first->second;
        if (inserted.second) {
            options.resize(code_.groups.size());
            for (auto &group_options : options) {
                group_options.by_extra.assign(static_cast<std::size_t>(h_ + 1), nullptr);
            }
        }

        auto &by_extra = options[static_cast<std::size_t>(group_index)].by_extra;
        if (by_extra[static_cast<std::size_t>(block.extra)] == nullptr) {
            by_extra[static_cast<std::size_t>(block.extra)] = &block;
        }
    }

    bool enumerate_annihilator_suffix(const std::vector<uint8_t> &basis, int dim, int pos,
                                      std::vector<uint8_t> *coefficients, int group_index,
                                      const ResidualBlock &block)
    {
        if (should_stop()) {
            return false;
        }
        if (pos == dim) {
            add_dual_index_option(group_index, block, basis, dim, *coefficients);
            return true;
        }

        for (int value = 0; value < 256; value++) {
            (*coefficients)[static_cast<std::size_t>(pos)] = static_cast<uint8_t>(value);
            if (!enumerate_annihilator_suffix(basis, dim, pos + 1, coefficients, group_index, block)) {
                return false;
            }
        }
        return true;
    }

    bool index_annihilator_projective_vectors(int group_index, const ResidualBlock &block)
    {
        std::vector<uint8_t> basis = left_nullspace_basis(block);
        int dim = h_ == 0 ? 0 : static_cast<int>(basis.size()) / h_;
        if (dim == 0) {
            return true;
        }

        std::vector<uint8_t> coefficients(static_cast<std::size_t>(dim), 0);
        for (int pivot = 0; pivot < dim; pivot++) {
            std::fill(coefficients.begin(), coefficients.end(), 0);
            coefficients[static_cast<std::size_t>(pivot)] = 1;
            if (!enumerate_annihilator_suffix(basis, dim, pivot + 1, &coefficients, group_index, block)) {
                return false;
            }
        }
        return true;
    }

    bool build_dual_index()
    {
        dual_index_.clear();
        for (int group_index = 0; group_index < static_cast<int>(residual_blocks_.size()); group_index++) {
            const auto &group_blocks = residual_blocks_[static_cast<std::size_t>(group_index)];
            for (int extra = 1; extra < static_cast<int>(group_blocks.by_extra.size()) && extra <= h_;
                 extra++) {
                for (const auto &block : group_blocks.by_extra[static_cast<std::size_t>(extra)]) {
                    if (!index_annihilator_projective_vectors(group_index, block)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool verify_by_dual_index()
    {
        if (!build_dual_index()) {
            return false;
        }

        for (const auto &entry : dual_index_) {
            std::vector<uint8_t> y = dual_vector_from_key(entry.first);
            if (!check_dual_options(y, entry.second)) {
                return false;
            }
        }
        return true;
    }

    bool verify_from_group(int group_index, int used, const std::vector<uint8_t> &current)
    {
        if (should_stop()) {
            return false;
        }
        if (group_index == static_cast<int>(code_.groups.size())) {
            return verify_global_completions(current, used);
        }

        std::string key = state_key(group_index, used, current);
        if (verified_states_.find(key) != verified_states_.end()) {
            return true;
        }

        const auto &group_blocks = residual_blocks_[static_cast<std::size_t>(group_index)];
        int max_extra = std::min(h_ - used, static_cast<int>(group_blocks.by_extra.size()) - 1);
        for (int extra = max_extra; extra >= 1; extra--) {
            const auto &blocks = group_blocks.by_extra[static_cast<std::size_t>(extra)];
            for (const auto &block : blocks) {
                group_erased_[static_cast<std::size_t>(group_index)] = block.erased_symbols;
                std::vector<uint8_t> next = append_columns(current, h_, used, block.matrix, extra);
                result_->patterns_checked++;
                if (gf256_rank(next, h_, used + extra) < used + extra) {
                    mark_failure_with_erased(group_erased_, first_kept_global_rows(used + extra));
                    group_erased_[static_cast<std::size_t>(group_index)].clear();
                    return false;
                }

                if (!verify_from_group(group_index + 1, used + extra, next)) {
                    group_erased_[static_cast<std::size_t>(group_index)].clear();
                    return false;
                }
                group_erased_[static_cast<std::size_t>(group_index)].clear();
            }
            if (should_stop()) {
                return false;
            }
        }

        if (!verify_from_group(group_index + 1, used, current)) {
            return false;
        }

        verified_states_.insert(std::move(key));
        return true;
    }

    const Code &code_;
    CheckResult *result_;
    std::vector<GroupResidualBlocks> residual_blocks_;
    std::vector<std::vector<int>> group_erased_;
    std::set<std::string> verified_states_;
    std::unordered_map<std::string, std::vector<DualGroupOptions>> dual_index_;
    int h_ = 0;
};

class FullErasureGate {
public:
    explicit FullErasureGate(const Code &code, CheckResult *result)
        : code_(code), result_(result), erased_(static_cast<std::size_t>(code.symbols), 0),
          choose_counts_(code.groups.size(), 0)
    {
    }

    void run()
    {
        enumerate_extra_counts(0, 0);
        result_->is_mr = result_->failures == 0;
        result_->message = result_->is_mr
                               ? "candidate passed the full maximal erasure gate"
                               : "candidate failed the full maximal erasure gate";
    }

private:
    bool should_stop() const
    {
        return result_->failures > 0;
    }

    void check_current_pattern()
    {
        if (should_stop()) {
            return;
        }

        std::vector<uint8_t> submatrix;
        submatrix.reserve(static_cast<std::size_t>(code_.data) * static_cast<std::size_t>(code_.data));
        int survivors = 0;
        for (int row = 0; row < code_.symbols; row++) {
            if (erased_[static_cast<std::size_t>(row)] != 0) {
                continue;
            }
            survivors++;
            for (int col = 0; col < code_.data; col++) {
                submatrix.push_back(code_.matrix[idx(row, col, code_.data)]);
            }
        }

        result_->patterns_checked++;
        if (survivors != code_.data || !gf256_matrix_is_invertible(submatrix, code_.data)) {
            result_->failures = 1;
            result_->first_failed_erased = erased_symbols_from_mask(code_, erased_);
        }
    }

    void enumerate_global_combinations(int start, int chosen)
    {
        if (should_stop()) {
            return;
        }
        if (chosen == global_choose_) {
            check_current_pattern();
            return;
        }

        int remaining = global_choose_ - chosen;
        for (int idx_global = start; idx_global <= code_.global_parity - remaining; idx_global++) {
            int symbol = code_.data + code_.local_rows + idx_global;
            erased_[static_cast<std::size_t>(symbol)] = 1;
            enumerate_global_combinations(idx_global + 1, chosen + 1);
            erased_[static_cast<std::size_t>(symbol)] = 0;
            if (should_stop()) {
                return;
            }
        }
    }

    void enumerate_group_combinations(int group_index, const std::vector<int> &symbols, int start, int chosen)
    {
        if (should_stop()) {
            return;
        }

        int need = choose_counts_[static_cast<std::size_t>(group_index)];
        if (chosen == need) {
            enumerate_group_choices(group_index + 1);
            return;
        }

        int remaining = need - chosen;
        for (int i = start; i <= static_cast<int>(symbols.size()) - remaining; i++) {
            erased_[static_cast<std::size_t>(symbols[static_cast<std::size_t>(i)])] = 1;
            enumerate_group_combinations(group_index, symbols, i + 1, chosen + 1);
            erased_[static_cast<std::size_t>(symbols[static_cast<std::size_t>(i)])] = 0;
            if (should_stop()) {
                return;
            }
        }
    }

    void enumerate_group_choices(int group_index)
    {
        if (should_stop()) {
            return;
        }
        if (group_index == static_cast<int>(code_.groups.size())) {
            enumerate_global_combinations(0, 0);
            return;
        }

        std::vector<int> symbols = group_symbol_ids(code_, group_index);
        enumerate_group_combinations(group_index, symbols, 0, 0);
    }

    void enumerate_extra_counts(int group_index, int extra_used)
    {
        if (should_stop()) {
            return;
        }
        if (group_index == static_cast<int>(code_.groups.size())) {
            global_choose_ = code_.global_parity - extra_used;
            if (global_choose_ >= 0 && global_choose_ <= code_.global_parity) {
                enumerate_group_choices(0);
            }
            return;
        }

        std::vector<int> symbols = group_symbol_ids(code_, group_index);
        const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
        int max_extra = std::min(static_cast<int>(symbols.size()) - group.local_parity,
                                 code_.global_parity - extra_used);
        for (int extra = 0; extra <= max_extra; extra++) {
            choose_counts_[static_cast<std::size_t>(group_index)] = group.local_parity + extra;
            enumerate_extra_counts(group_index + 1, extra_used + extra);
            if (should_stop()) {
                return;
            }
        }
    }

    const Code &code_;
    CheckResult *result_;
    std::vector<uint8_t> erased_;
    std::vector<int> choose_counts_;
    int global_choose_ = 0;
};

} // namespace

const char *family_name(MatrixFamily family)
{
    switch (family) {
    case MatrixFamily::Cauchy:
        return "cauchy";
    case MatrixFamily::Vandermonde:
        return "vandermonde";
    case MatrixFamily::Random:
        return "random";
    }
    return "unknown";
}

bool parse_family(const std::string &name, MatrixFamily *family)
{
    if (name == "cauchy") {
        *family = MatrixFamily::Cauchy;
        return true;
    }
    if (name == "vandermonde" || name == "rs") {
        *family = MatrixFamily::Vandermonde;
        return true;
    }
    if (name == "random" || name == "rand") {
        *family = MatrixFamily::Random;
        return true;
    }
    return false;
}

std::vector<int> erased_symbols_from_mask(const Code &code, const std::vector<uint8_t> &erased)
{
    std::vector<int> symbols;
    for (int symbol = 0; symbol < code.symbols; symbol++) {
        if (erased[static_cast<std::size_t>(symbol)] != 0) {
            symbols.push_back(symbol);
        }
    }
    return symbols;
}

std::string symbol_label(const Code &code, int symbol)
{
    if (symbol < code.data) {
        return "D" + std::to_string(symbol);
    }

    for (std::size_t group_index = 0; group_index < code.groups.size(); group_index++) {
        const auto &group = code.groups[group_index];
        if (symbol >= group.local_row_start && symbol < group.local_row_start + group.local_parity) {
            return "L" + std::to_string(group_index) + "_" +
                   std::to_string(symbol - group.local_row_start);
        }
    }

    return "G" + std::to_string(symbol - code.data - code.local_rows);
}

CheckResult check_mr(const Code &code)
{
    CheckResult result;
    ResidualVerifier verifier(code, &result);
    verifier.run();
    if (!result.is_mr) {
        return result;
    }

    CheckResult gate;
    FullErasureGate full_gate(code, &gate);
    full_gate.run();
    result.patterns_checked += gate.patterns_checked;
    if (!gate.is_mr) {
        result.is_mr = false;
        result.failures = gate.failures;
        result.first_failed_erased = std::move(gate.first_failed_erased);
        result.message = "residual checker passed, but the full maximal erasure gate failed";
        return result;
    }

    result.message =
        "candidate passed residual verification and the full maximal erasure gate";
    return result;
}

GenerateResult generate(const Params &params)
{
    GenerateResult result;
    try {
        Code base_code = make_empty_code(params);
        result.status = 2;

        std::atomic<uint64_t> next_attempt{1};
        std::atomic<uint64_t> completed_attempts{0};
        std::atomic<uint64_t> total_patterns_checked{0};
        std::atomic<uint64_t> failed_candidates{0};
        std::atomic<bool> stop{false};
        std::atomic<bool> progress_done{false};
        std::mutex result_mutex;
        std::mutex progress_mutex;
        std::condition_variable progress_cv;
        std::string worker_error;
        std::string progress_error;

        std::thread progress_worker;
        if (params.progress_callback && params.step_time != 0) {
            progress_worker = std::thread([&]() {
                std::unique_lock<std::mutex> lock(progress_mutex);
                for (;;) {
                    if (progress_cv.wait_for(lock, std::chrono::seconds(params.step_time),
                                             [&]() { return progress_done.load(std::memory_order_relaxed); })) {
                        break;
                    }

                    uint64_t completed = completed_attempts.load(std::memory_order_relaxed);
                    lock.unlock();
                    try {
                        params.progress_callback(completed);
                    } catch (const std::exception &ex) {
                        {
                            std::lock_guard<std::mutex> result_lock(result_mutex);
                            if (progress_error.empty()) {
                                progress_error = ex.what();
                            }
                        }
                        stop.store(true, std::memory_order_relaxed);
                        progress_done.store(true, std::memory_order_relaxed);
                        progress_cv.notify_all();
                        break;
                    } catch (...) {
                        {
                            std::lock_guard<std::mutex> result_lock(result_mutex);
                            if (progress_error.empty()) {
                                progress_error = "unknown progress callback error";
                            }
                        }
                        stop.store(true, std::memory_order_relaxed);
                        progress_done.store(true, std::memory_order_relaxed);
                        progress_cv.notify_all();
                        break;
                    }
                    lock.lock();
                }
            });
        }

        auto run_attempt = [&](uint64_t attempt) {
            Code code = base_code;
            code.attempt = attempt;
            Rng rng(seed_for_attempt(params.seed, attempt));
            build_candidate(&code, rng);
            CheckResult check = check_mr(code);
            total_patterns_checked.fetch_add(check.patterns_checked, std::memory_order_relaxed);
            completed_attempts.fetch_add(1, std::memory_order_relaxed);

            if (check.is_mr) {
                bool expected = false;
                if (stop.compare_exchange_strong(expected, true)) {
                    code.patterns_checked = check.patterns_checked;
                    std::lock_guard<std::mutex> lock(result_mutex);
                    result.status = 0;
                    result.message = "MR-LRC matrix found on attempt " + std::to_string(attempt);
                    result.code = std::move(code);
                    result.check = std::move(check);
                }
                return;
            }

            failed_candidates.fetch_add(1, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(result_mutex);
            if (result.code.matrix.empty()) {
                code.patterns_checked = check.patterns_checked;
                result.code = std::move(code);
                result.check = std::move(check);
            }
        };

        uint64_t worker_count = std::min(params.thread_count, params.random_limit);
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(worker_count));

        for (uint64_t worker = 0; worker < worker_count; worker++) {
            workers.emplace_back([&]() {
                try {
                    for (;;) {
                        if (stop.load(std::memory_order_relaxed)) {
                            break;
                        }
                        uint64_t attempt = next_attempt.fetch_add(1, std::memory_order_relaxed);
                        if (attempt > params.random_limit) {
                            break;
                        }
                        run_attempt(attempt);
                    }
                } catch (const std::exception &ex) {
                    {
                        std::lock_guard<std::mutex> lock(result_mutex);
                        if (worker_error.empty()) {
                            worker_error = ex.what();
                        }
                    }
                    stop.store(true, std::memory_order_relaxed);
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lock(result_mutex);
                        if (worker_error.empty()) {
                            worker_error = "unknown worker error";
                        }
                    }
                    stop.store(true, std::memory_order_relaxed);
                }
            });
        }

        for (auto &worker : workers) {
            worker.join();
        }

        progress_done.store(true, std::memory_order_relaxed);
        progress_cv.notify_all();
        if (progress_worker.joinable()) {
            progress_worker.join();
        }

        if (!progress_error.empty()) {
            result.status = -1;
            result.message = progress_error;
            return result;
        }

        if (!worker_error.empty()) {
            result.status = -1;
            result.message = worker_error;
            return result;
        }

        result.attempts_done = completed_attempts.load(std::memory_order_relaxed);
        if (result.status == 0) {
            return result;
        }

        result.status = 2;
        result.message = "no MR-LRC matrix found after " +
                         std::to_string(completed_attempts.load(std::memory_order_relaxed)) + " attempts";
        result.check.patterns_checked = total_patterns_checked.load(std::memory_order_relaxed);
        result.check.failures = failed_candidates.load(std::memory_order_relaxed);
        if (!result.code.matrix.empty()) {
            result.code.patterns_checked = result.check.patterns_checked;
        }
        return result;
    } catch (const std::exception &ex) {
        result.status = -1;
        result.message = ex.what();
        return result;
    }
}

} // namespace mrlrc
