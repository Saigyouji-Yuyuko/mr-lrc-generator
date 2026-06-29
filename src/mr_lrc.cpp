#include "mr_lrc.h"

#include "gf256.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <cstring>
#include <initializer_list>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <ankerl/unordered_dense.h>

namespace mrlrc {
namespace {

constexpr uint64_t kMaxThreadCount = 256;
constexpr uint64_t kH4StreamingDualLimit = 8192;

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

template <typename T, std::size_t N>
void shuffle_array_prefix(Rng &rng, std::array<T, N> *items, std::size_t count)
{
    for (std::size_t i = count; i > 1; i--) {
        std::size_t j = static_cast<std::size_t>(rng.range(i));
        std::swap((*items)[i - 1], (*items)[j]);
    }
}

struct CauchyKeyBytes {
    std::array<uint8_t, 256> bytes{};
    std::size_t size = 0;

    void clear() { size = 0; }

    void append(const uint8_t *data, std::size_t count)
    {
        if (count > bytes.size() - size) {
            throw ValidationError("cauchy dedup key exceeds 256 bytes; use --no-cauchy-dedup");
        }
        std::memcpy(bytes.data() + size, data, count);
        size += count;
    }
};

template <std::size_t Bytes>
struct Key {
    std::array<uint8_t, Bytes> bytes{};
    uint16_t size = 0;

    explicit Key(const CauchyKeyBytes &key)
    {
        if (key.size > Bytes) {
            throw ValidationError("cauchy dedup key does not fit selected fixed key size");
        }
        size = static_cast<uint16_t>(key.size);
        std::memcpy(bytes.data(), key.bytes.data(), key.size);
    }

    bool operator==(const Key &other) const
    {
        return size == other.size &&
               std::memcmp(bytes.data(), other.bytes.data(), size) == 0;
    }
};

template <std::size_t Bytes>
struct KeyHash {
    using is_avalanching = void;

    uint64_t operator()(const Key<Bytes> &key) const
    {
        uint64_t h = UINT64_C(0x8ddca3a1f0f23a6d) ^ key.size;
        for (std::size_t i = 0; i < key.size; i++) {
            h = splitmix64_value(h ^ key.bytes[i]);
        }
        return h;
    }
};

class CauchyDeduper {
public:
    virtual ~CauchyDeduper() = default;
    virtual bool insert(const CauchyKeyBytes &key) = 0;
};

template <std::size_t Bytes>
class DenseCauchyDeduper final : public CauchyDeduper {
public:
    bool insert(const CauchyKeyBytes &key) override
    {
        return keys_.insert(Key<Bytes>(key)).second;
    }

private:
    ankerl::unordered_dense::set<Key<Bytes>, KeyHash<Bytes>> keys_;
};

std::unique_ptr<CauchyDeduper> make_cauchy_deduper(uint64_t key_bytes)
{
    if (key_bytes <= 32) {
        return std::make_unique<DenseCauchyDeduper<32>>();
    }
    if (key_bytes <= 64) {
        return std::make_unique<DenseCauchyDeduper<64>>();
    }
    if (key_bytes <= 128) {
        return std::make_unique<DenseCauchyDeduper<128>>();
    }
    if (key_bytes <= 256) {
        return std::make_unique<DenseCauchyDeduper<256>>();
    }
    throw ValidationError("cauchy dedup key exceeds 256 bytes; use --no-cauchy-dedup");
}

void append_cauchy_canonical_key(const std::array<uint8_t, 256> &pool, int rows, int cols,
                                 CauchyKeyBytes *key)
{
    if (key == nullptr || rows == 0 || cols == 0) {
        return;
    }

    auto block_size = static_cast<std::size_t>(rows + cols);
    std::array<uint8_t, 256> best{};
    std::array<uint8_t, 256> candidate{};
    std::array<uint8_t, 256> normalized_rows{};
    bool have_best = false;

    for (int shift = 0; shift < 256; shift++) {
        uint8_t t = static_cast<uint8_t>(shift);
        for (int row = 0; row < rows; row++) {
            normalized_rows[static_cast<std::size_t>(row)] =
                static_cast<uint8_t>(pool[static_cast<std::size_t>(row)] ^ t);
        }
        std::sort(normalized_rows.begin(), normalized_rows.begin() + rows);

        std::copy(normalized_rows.begin(), normalized_rows.begin() + rows, candidate.begin());
        for (int col = 0; col < cols; col++) {
            candidate[static_cast<std::size_t>(rows + col)] =
                static_cast<uint8_t>(pool[static_cast<std::size_t>(rows + col)] ^ t);
        }

        if (!have_best ||
            std::lexicographical_compare(candidate.begin(), candidate.begin() + block_size,
                                         best.begin(), best.begin() + block_size)) {
            std::copy(candidate.begin(), candidate.begin() + block_size, best.begin());
            have_best = true;
        }
    }

    key->append(best.data() + 1, block_size - 1);
}

template <typename SetValue>
void fill_cauchy_entries(int rows, int cols, Rng &rng, CauchyKeyBytes *canonical_key,
                         SetValue &&set_value)
{
    if (rows < 0 || cols < 0 || rows + cols > 256) {
        throw ValidationError("cauchy matrix needs rows + cols <= 256 over GF(256)");
    }

    std::array<uint8_t, 256> pool{};
    for (int i = 0; i < 256; i++) {
        pool[static_cast<std::size_t>(i)] = static_cast<uint8_t>(i);
    }
    shuffle_array_prefix(rng, &pool, pool.size());
    append_cauchy_canonical_key(pool, rows, cols, canonical_key);

    for (int row = 0; row < rows; row++) {
        uint8_t x = pool[static_cast<std::size_t>(row)];
        for (int col = 0; col < cols; col++) {
            uint8_t y = pool[static_cast<std::size_t>(rows + col)];
            set_value(row, col, gf256_inv(static_cast<uint8_t>(x ^ y)));
        }
    }
}

template <typename SetValue>
void fill_column_multiplier_cauchy_entries(int rows, int cols, Rng &rng, SetValue &&set_value)
{
    if (rows < 0 || cols < 0 || rows + cols > 256) {
        throw ValidationError("cauchy matrix needs rows + cols <= 256 over GF(256)");
    }

    std::array<uint8_t, 256> pool{};
    for (int i = 0; i < 256; i++) {
        pool[static_cast<std::size_t>(i)] = static_cast<uint8_t>(i);
    }
    shuffle_array_prefix(rng, &pool, pool.size());

    std::array<uint8_t, 256> column_scale{};
    for (int col = 0; col < cols; col++) {
        column_scale[static_cast<std::size_t>(col)] = rng.nonzero_byte();
    }

    for (int row = 0; row < rows; row++) {
        uint8_t x = pool[static_cast<std::size_t>(row)];
        for (int col = 0; col < cols; col++) {
            uint8_t y = pool[static_cast<std::size_t>(rows + col)];
            uint8_t value = gf256_inv(static_cast<uint8_t>(x ^ y));
            set_value(row, col, gf256_mul(value, column_scale[static_cast<std::size_t>(col)]));
        }
    }
}

template <typename SetValue>
void fill_vandermonde_entries(int rows, int cols, Rng &rng, SetValue &&set_value)
{
    if (rows < 0 || cols < 0 || rows > 255 || cols > 255) {
        throw ValidationError("vandermonde matrix needs rows, cols <= 255 over GF(256)");
    }

    std::array<uint8_t, 256> points{};
    for (int i = 0; i < 255; i++) {
        points[static_cast<std::size_t>(i)] = static_cast<uint8_t>(i + 1);
    }
    shuffle_array_prefix(rng, &points, 255);

    std::array<uint8_t, 256> column_scale{};
    for (int col = 0; col < cols; col++) {
        column_scale[static_cast<std::size_t>(col)] = rng.nonzero_byte();
    }

    for (int row = 0; row < rows; row++) {
        uint8_t x = points[static_cast<std::size_t>(row)];
        for (int col = 0; col < cols; col++) {
            uint8_t v = gf256_pow(x, static_cast<unsigned int>(col));
            set_value(row, col, gf256_mul(column_scale[static_cast<std::size_t>(col)], v));
        }
    }
}

template <typename SetValue>
void fill_random_entries(int rows, int cols, Rng &rng, bool nonzero, SetValue &&set_value)
{
    if (rows < 0 || cols < 0) {
        throw ValidationError("random matrix dimensions must be non-negative");
    }
    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            set_value(row, col, nonzero ? rng.nonzero_byte() : rng.byte());
        }
    }
}

template <typename SetValue>
void fill_dense_matrix_entries(MatrixFamily family, int rows, int cols, Rng &rng,
                               CauchyKeyBytes *cauchy_key, SetValue &&set_value)
{
    switch (family) {
    case MatrixFamily::Cauchy:
        fill_cauchy_entries(rows, cols, rng, cauchy_key, std::forward<SetValue>(set_value));
        return;
    case MatrixFamily::ColumnMultiplierCauchy:
        fill_column_multiplier_cauchy_entries(rows, cols, rng,
                                              std::forward<SetValue>(set_value));
        return;
    case MatrixFamily::Vandermonde:
        fill_vandermonde_entries(rows, cols, rng, std::forward<SetValue>(set_value));
        return;
    case MatrixFamily::Random:
        fill_random_entries(rows, cols, rng, false, std::forward<SetValue>(set_value));
        return;
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
        if (params.local_family == MatrixFamily::ColumnMultiplierCauchy) {
            throw ValidationError("column_multiplier_cauchy is only supported for --global-method");
        }
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

    if ((code.global_family == MatrixFamily::Cauchy ||
         code.global_family == MatrixFamily::ColumnMultiplierCauchy) &&
        code.global_parity > 0 && code.data + code.global_parity > 256) {
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

void build_local_rows(Code *code, Rng &rng, CauchyKeyBytes *cauchy_key)
{
    for (const auto &group : code->groups) {
        int group_data = static_cast<int>(group.data.size());
        auto set_local = [&](int local_row, int col, uint8_t value) {
            int matrix_row = group.local_row_start + local_row;
            int data_col = group.data[static_cast<std::size_t>(col)];
            code->matrix[idx(matrix_row, data_col, code->data)] = value;
        };

        if (code->local_family == MatrixFamily::Random) {
            fill_random_entries(group.local_parity, group_data, rng, true, set_local);
        } else {
            fill_dense_matrix_entries(code->local_family, group.local_parity, group_data, rng,
                                      code->local_family == MatrixFamily::Cauchy ? cauchy_key
                                                                                 : nullptr,
                                      set_local);
        }
    }
}

void build_global_rows(Code *code, Rng &rng, CauchyKeyBytes *cauchy_key)
{
    int first_global = code->data + code->local_rows;
    fill_dense_matrix_entries(
        code->global_family, code->global_parity, code->data, rng,
        code->global_family == MatrixFamily::Cauchy ? cauchy_key : nullptr,
        [&](int row, int col, uint8_t value) {
            code->matrix[idx(first_global + row, col, code->data)] = value;
        });
}

void build_candidate(Code *code, Rng &rng, CauchyKeyBytes *cauchy_key)
{
    reset_identity(code);
    if (cauchy_key != nullptr) {
        cauchy_key->clear();
    }
    build_local_rows(code, rng, cauchy_key);
    build_global_rows(code, rng, cauchy_key);
}

bool use_cauchy_dedup(const Params &params)
{
    return params.cauchy_dedup && params.local_family == MatrixFamily::Cauchy &&
           params.global_family == MatrixFamily::Cauchy;
}

uint64_t cauchy_dedup_key_bytes(const Code &code)
{
    uint64_t bytes = 0;
    for (const auto &group : code.groups) {
        if (group.local_parity > 0 && !group.data.empty()) {
            bytes += static_cast<uint64_t>(group.local_parity) +
                     static_cast<uint64_t>(group.data.size()) - 1;
        }
    }
    if (code.global_parity > 0 && code.data > 0) {
        bytes += static_cast<uint64_t>(code.global_parity) +
                 static_cast<uint64_t>(code.data) - 1;
    }
    return bytes;
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

bool erased_pattern_is_recoverable(const Code &code, const std::vector<uint8_t> &erased)
{
    std::vector<uint8_t> submatrix;
    submatrix.reserve(static_cast<std::size_t>(code.data) * static_cast<std::size_t>(code.data));
    int survivors = 0;
    for (int row = 0; row < code.symbols; row++) {
        if (erased[static_cast<std::size_t>(row)] != 0) {
            continue;
        }
        survivors++;
        for (int col = 0; col < code.data; col++) {
            submatrix.push_back(code.matrix[idx(row, col, code.data)]);
        }
    }

    return survivors == code.data && gf256_matrix_is_invertible(submatrix, code.data);
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

int gf256_rank_small(std::array<uint8_t, 16> matrix, int rows, int cols)
{
    if (rows <= 0 || cols <= 0) {
        return 0;
    }

    int rank = 0;
    for (int col = 0; col < cols && rank < rows; col++) {
        int pivot = -1;
        for (int row = rank; row < rows; row++) {
            if (matrix[idx(row, col, cols)] != 0) {
                pivot = row;
                break;
            }
        }
        if (pivot < 0) {
            continue;
        }

        if (pivot != rank) {
            for (int c = 0; c < cols; c++) {
                std::swap(matrix[idx(rank, c, cols)], matrix[idx(pivot, c, cols)]);
            }
        }

        uint8_t inv_pivot = gf256_inv(matrix[idx(rank, col, cols)]);
        for (int c = 0; c < cols; c++) {
            matrix[idx(rank, c, cols)] = gf256_mul(matrix[idx(rank, c, cols)], inv_pivot);
        }

        for (int row = 0; row < rows; row++) {
            if (row == rank) {
                continue;
            }
            uint8_t factor = matrix[idx(row, col, cols)];
            if (factor == 0) {
                continue;
            }
            for (int c = 0; c < cols; c++) {
                matrix[idx(row, c, cols)] ^= gf256_mul(factor, matrix[idx(rank, c, cols)]);
            }
        }
        rank++;
    }
    return rank;
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

bool projective_column_key_packed(const std::array<uint8_t, 16> &matrix, int rows, int cols,
                                  int column, uint32_t *key)
{
    if (key == nullptr || rows <= 0 || rows > 4) {
        return false;
    }

    int pivot = -1;
    for (int row = 0; row < rows; row++) {
        if (matrix[static_cast<std::size_t>(idx(row, column, cols))] != 0) {
            pivot = row;
            break;
        }
    }
    if (pivot < 0) {
        return false;
    }

    uint8_t inv_pivot =
        gf256_inv(matrix[static_cast<std::size_t>(idx(pivot, column, cols))]);
    uint32_t packed = 0;
    for (int row = 0; row < rows; row++) {
        uint8_t value =
            gf256_mul(matrix[static_cast<std::size_t>(idx(row, column, cols))], inv_pivot);
        packed |= static_cast<uint32_t>(value) << (8U * static_cast<unsigned int>(row));
    }
    *key = packed;
    return true;
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

std::vector<uint8_t> build_local_erasure_submatrix(const Code &code, int group_index,
                                                   const std::vector<int> &symbols)
{
    const auto &group = code.groups[static_cast<std::size_t>(group_index)];
    int local_rows = group.local_parity;
    std::vector<uint8_t> matrix(static_cast<std::size_t>(local_rows) *
                                    static_cast<std::size_t>(symbols.size()),
                                0);
    for (int col = 0; col < static_cast<int>(symbols.size()); col++) {
        int symbol = symbols[static_cast<std::size_t>(col)];
        if (symbol < code.data) {
            for (int local = 0; local < local_rows; local++) {
                matrix[idx(local, col, static_cast<int>(symbols.size()))] =
                    code.matrix[idx(group.local_row_start + local, symbol, code.data)];
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

std::vector<uint8_t> build_global_erasure_submatrix(const Code &code, const std::vector<int> &symbols)
{
    int global_rows = code.global_parity;
    std::vector<uint8_t> matrix(static_cast<std::size_t>(global_rows) *
                                    static_cast<std::size_t>(symbols.size()),
                                0);
    int first_global = code.data + code.local_rows;
    for (int col = 0; col < static_cast<int>(symbols.size()); col++) {
        int symbol = symbols[static_cast<std::size_t>(col)];
        if (symbol >= code.data) {
            continue;
        }
        for (int global = 0; global < global_rows; global++) {
            matrix[idx(global, col, static_cast<int>(symbols.size()))] =
                code.matrix[idx(first_global + global, symbol, code.data)];
        }
    }
    return matrix;
}

uint8_t local_erasure_coefficient(const Code &code, int group_index, int symbol)
{
    const auto &group = code.groups[static_cast<std::size_t>(group_index)];
    if (symbol < code.data) {
        return code.matrix[idx(group.local_row_start, symbol, code.data)];
    }
    if (symbol == group.local_row_start) {
        return 1;
    }
    return 0;
}

uint8_t global_erasure_coefficient(const Code &code, int global_row, int symbol)
{
    if (symbol >= code.data) {
        return 0;
    }
    int first_global = code.data + code.local_rows;
    return code.matrix[idx(first_global + global_row, symbol, code.data)];
}

std::vector<uint8_t> build_residual_erasure_block_a1(const Code &code, int group_index,
                                                     const std::vector<int> &erased_symbols,
                                                     int extra)
{
    int global_rows = code.global_parity;
    std::vector<uint8_t> residual(static_cast<std::size_t>(global_rows) *
                                      static_cast<std::size_t>(extra),
                                  0);

    int pivot = erased_symbols.front();
    uint8_t inv_pivot = gf256_inv(local_erasure_coefficient(code, group_index, pivot));
    for (int col = 0; col < extra; col++) {
        int symbol = erased_symbols[static_cast<std::size_t>(col + 1)];
        uint8_t scale = gf256_mul(local_erasure_coefficient(code, group_index, symbol), inv_pivot);
        for (int row = 0; row < global_rows; row++) {
            residual[idx(row, col, extra)] =
                global_erasure_coefficient(code, row, symbol) ^
                gf256_mul(scale, global_erasure_coefficient(code, row, pivot));
        }
    }
    return residual;
}

void build_residual_erasure_block_a1_small(const Code &code, int group_index,
                                           const std::vector<int> &erased_symbols, int extra,
                                           std::array<uint8_t, 16> *residual)
{
    residual->fill(0);
    int pivot = erased_symbols.front();
    uint8_t inv_pivot = gf256_inv(local_erasure_coefficient(code, group_index, pivot));
    for (int col = 0; col < extra; col++) {
        int symbol = erased_symbols[static_cast<std::size_t>(col + 1)];
        uint8_t scale = gf256_mul(local_erasure_coefficient(code, group_index, symbol), inv_pivot);
        for (int row = 0; row < code.global_parity; row++) {
            (*residual)[static_cast<std::size_t>(idx(row, col, extra))] =
                global_erasure_coefficient(code, row, symbol) ^
                gf256_mul(scale, global_erasure_coefficient(code, row, pivot));
        }
    }
}

std::vector<uint8_t> build_residual_erasure_block(const Code &code, int group_index,
                                                  const std::vector<int> &erased_symbols,
                                                  int extra)
{
    const auto &group = code.groups[static_cast<std::size_t>(group_index)];
    if (group.local_parity == 1) {
        return build_residual_erasure_block_a1(code, group_index, erased_symbols, extra);
    }

    std::vector<uint8_t> local = build_local_erasure_submatrix(code, group_index, erased_symbols);
    std::vector<uint8_t> nullspace =
        gf256_nullspace_basis(local, group.local_parity, static_cast<int>(erased_symbols.size()));
    std::vector<uint8_t> global = build_global_erasure_submatrix(code, erased_symbols);
    std::vector<uint8_t> residual(static_cast<std::size_t>(code.global_parity) *
                                      static_cast<std::size_t>(extra),
                                  0);

    for (int row = 0; row < code.global_parity; row++) {
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
    std::array<uint8_t, 16> small_matrix{};
    int small_rows = 0;
    int small_cols = 0;
    bool has_small_matrix = false;
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
            } else if (code_.groups.size() == 2 && (h_ == 3 || h_ == 4)) {
                verify_two_group_projection_fast_path();
            } else if (h_ == 3) {
                verify_by_streaming_dual(0);
            } else if (h_ == 4 && !verify_by_streaming_dual(kH4StreamingDualLimit)) {
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
        return build_local_erasure_submatrix(code_, group_index, symbols);
    }

    std::vector<uint8_t> build_residual_block(int group_index, const std::vector<int> &erased_symbols,
                                              int extra) const
    {
        return build_residual_erasure_block(code_, group_index, erased_symbols, extra);
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

    struct ProjectedBlock {
        const ResidualBlock *source = nullptr;
        std::vector<uint8_t> matrix;
    };

    bool projected_blocks_for_group(int group_index, int extra, const std::vector<int> &kept_rows,
                                    std::vector<ProjectedBlock> *projected)
    {
        projected->clear();
        if (extra == 0) {
            projected->push_back(ProjectedBlock{nullptr, {}});
            return true;
        }

        const auto &group_blocks = residual_blocks_[static_cast<std::size_t>(group_index)];
        if (extra >= static_cast<int>(group_blocks.by_extra.size())) {
            return true;
        }

        int rows = static_cast<int>(kept_rows.size());
        projected->reserve(group_blocks.by_extra[static_cast<std::size_t>(extra)].size());
        for (const auto &block : group_blocks.by_extra[static_cast<std::size_t>(extra)]) {
            std::vector<uint8_t> matrix = project_rows(block.matrix, h_, extra, kept_rows);
            result_->patterns_checked++;
            if (gf256_rank(matrix, rows, extra) < extra) {
                std::vector<std::vector<int>> group_erased(code_.groups.size());
                group_erased[static_cast<std::size_t>(group_index)] = block.erased_symbols;
                mark_failure_with_erased(group_erased, kept_rows);
                return false;
            }
            projected->push_back(ProjectedBlock{&block, std::move(matrix)});
        }
        return true;
    }

    bool check_two_group_projection_case(const std::vector<int> &kept_rows, int extra0, int extra1)
    {
        int rows = static_cast<int>(kept_rows.size());
        std::vector<ProjectedBlock> group0;
        std::vector<ProjectedBlock> group1;
        if (!projected_blocks_for_group(0, extra0, kept_rows, &group0) ||
            !projected_blocks_for_group(1, extra1, kept_rows, &group1)) {
            return false;
        }

        for (const auto &left : group0) {
            for (const auto &right : group1) {
                if (extra0 == 0 || extra1 == 0) {
                    continue;
                }

                std::vector<uint8_t> joined =
                    append_columns(left.matrix, rows, extra0, right.matrix, extra1);
                result_->patterns_checked++;
                if (gf256_rank(joined, rows, rows) < rows) {
                    std::vector<std::vector<int>> group_erased(code_.groups.size());
                    group_erased[0] = left.source->erased_symbols;
                    group_erased[1] = right.source->erased_symbols;
                    mark_failure_with_erased(group_erased, kept_rows);
                    return false;
                }
            }
        }
        return true;
    }

    bool check_two_group_rows(const std::vector<int> &kept_rows)
    {
        int residual_size = static_cast<int>(kept_rows.size());
        for (int extra0 = 0; extra0 <= residual_size; extra0++) {
            int extra1 = residual_size - extra0;
            if (!check_two_group_projection_case(kept_rows, extra0, extra1)) {
                return false;
            }
        }
        return true;
    }

    bool verify_two_group_projection_fast_path()
    {
        std::vector<int> rows;
        rows.reserve(static_cast<std::size_t>(h_));
        for (int row = 0; row < h_; row++) {
            rows.push_back(row);
        }

        for (int residual_size = 1; residual_size <= h_; residual_size++) {
            std::vector<int> kept_rows;
            if (!enumerate_combinations(rows, residual_size, 0, &kept_rows,
                                        [&](const std::vector<int> &choice) {
                                            if (should_stop()) {
                                                return false;
                                            }
                                            return check_two_group_rows(choice);
                                        })) {
                return false;
            }
        }
        return true;
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

    bool annihilates(const std::vector<uint8_t> &y, const ResidualBlock &block) const
    {
        for (int col = 0; col < block.extra; col++) {
            uint8_t value = 0;
            for (int row = 0; row < h_; row++) {
                value ^= gf256_mul(y[static_cast<std::size_t>(row)],
                                   block.matrix[idx(row, col, block.extra)]);
            }
            if (value != 0) {
                return false;
            }
        }
        return true;
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

    std::vector<DualGroupOptions> dual_options_for_vector(const std::vector<uint8_t> &y) const
    {
        std::vector<DualGroupOptions> options(code_.groups.size());
        for (auto &group_options : options) {
            group_options.by_extra.assign(static_cast<std::size_t>(h_ + 1), nullptr);
        }

        for (int group_index = 0; group_index < static_cast<int>(residual_blocks_.size()); group_index++) {
            const auto &group_blocks = residual_blocks_[static_cast<std::size_t>(group_index)];
            for (int extra = 1; extra < static_cast<int>(group_blocks.by_extra.size()) && extra <= h_;
                 extra++) {
                auto &slot = options[static_cast<std::size_t>(group_index)]
                                 .by_extra[static_cast<std::size_t>(extra)];
                if (slot != nullptr) {
                    continue;
                }
                for (const auto &block : group_blocks.by_extra[static_cast<std::size_t>(extra)]) {
                    if (annihilates(y, block)) {
                        slot = &block;
                        break;
                    }
                }
            }
        }

        return options;
    }

    bool check_streaming_dual_from_coefficients(const std::vector<uint8_t> &basis, int dim,
                                                const std::vector<uint8_t> &coefficients)
    {
        std::vector<uint8_t> y = dual_from_coefficients(basis, dim, coefficients);
        std::string key = canonical_dual_key(y);
        if (key.empty()) {
            return true;
        }
        if (checked_streaming_dual_vectors_.find(key) != checked_streaming_dual_vectors_.end()) {
            return true;
        }
        if (streaming_dual_limit_ != 0 &&
            checked_streaming_dual_vectors_.size() >= streaming_dual_limit_) {
            streaming_dual_limit_hit_ = true;
            return false;
        }

        checked_streaming_dual_vectors_.insert(key);
        return check_dual_options(y, dual_options_for_vector(y));
    }

    bool enumerate_streaming_annihilator_suffix(const std::vector<uint8_t> &basis, int dim, int pos,
                                                std::vector<uint8_t> *coefficients)
    {
        if (should_stop() || streaming_dual_limit_hit_) {
            return false;
        }
        if (pos == dim) {
            return check_streaming_dual_from_coefficients(basis, dim, *coefficients);
        }

        for (int value = 0; value < 256; value++) {
            (*coefficients)[static_cast<std::size_t>(pos)] = static_cast<uint8_t>(value);
            if (!enumerate_streaming_annihilator_suffix(basis, dim, pos + 1, coefficients)) {
                return false;
            }
        }
        return true;
    }

    bool stream_annihilator_projective_vectors(const ResidualBlock &block)
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
            if (!enumerate_streaming_annihilator_suffix(basis, dim, pivot + 1, &coefficients)) {
                return false;
            }
        }
        return true;
    }

    bool verify_by_streaming_dual(uint64_t limit)
    {
        streaming_dual_limit_ = limit;
        streaming_dual_limit_hit_ = false;
        checked_streaming_dual_vectors_.clear();

        for (int extra = h_; extra >= 1; extra--) {
            for (const auto &group_blocks : residual_blocks_) {
                if (extra >= static_cast<int>(group_blocks.by_extra.size())) {
                    continue;
                }
                for (const auto &block : group_blocks.by_extra[static_cast<std::size_t>(extra)]) {
                    if (!stream_annihilator_projective_vectors(block)) {
                        return should_stop() || !streaming_dual_limit_hit_;
                    }
                }
            }
        }

        return true;
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
            if (checked_streaming_dual_vectors_.find(entry.first) != checked_streaming_dual_vectors_.end()) {
                continue;
            }
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
    std::unordered_set<std::string> checked_streaming_dual_vectors_;
    uint64_t streaming_dual_limit_ = 0;
    bool streaming_dual_limit_hit_ = false;
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

        result_->patterns_checked++;
        if (!erased_pattern_is_recoverable(code_, erased_)) {
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

enum class StressProbeKind {
    SingleRank,
    LineCollision,
    VectorPlane,
    PlanePlane,
    HyperplaneVector,
    GenericRank,
    RandomResidual,
};

struct StressProbe {
    StressProbeKind kind = StressProbeKind::RandomResidual;
    std::array<int, 4> groups{-1, -1, -1, -1};
    std::array<int, 4> extras{0, 0, 0, 0};
    int terms = 0;
    std::vector<int> kept_rows;
    uint64_t scan_limit = 4096;
    uint64_t random_sample = 0;
};

struct PrefilterPlan {
    std::vector<StressProbe> probes;
    std::vector<std::vector<std::vector<std::vector<int>>>> erasure_choices;
};

std::vector<int> all_row_ids(int count)
{
    std::vector<int> rows;
    rows.reserve(static_cast<std::size_t>(count));
    for (int row = 0; row < count; row++) {
        rows.push_back(row);
    }
    return rows;
}

std::vector<std::vector<int>> row_subsets(int rows, int need)
{
    std::vector<std::vector<int>> subsets;
    std::vector<int> all = all_row_ids(rows);
    std::vector<int> chosen;
    enumerate_combinations(all, need, 0, &chosen, [&](const std::vector<int> &choice) {
        subsets.push_back(choice);
        return true;
    });
    return subsets;
}

StressProbe make_probe(StressProbeKind kind, const std::vector<int> &groups,
                       const std::vector<int> &extras, const std::vector<int> &kept_rows)
{
    StressProbe probe;
    probe.kind = kind;
    std::size_t terms = std::min<std::size_t>(4, std::min(groups.size(), extras.size()));
    probe.terms = static_cast<int>(terms);
    probe.kept_rows = kept_rows;
    for (std::size_t i = 0; i < terms; i++) {
        probe.groups[i] = groups[i];
        probe.extras[i] = extras[i];
    }
    return probe;
}

void add_probe_family_round_robin(const std::vector<std::vector<StressProbe>> &families,
                                  uint64_t limit, std::vector<StressProbe> *out)
{
    if (out == nullptr) {
        return;
    }

    std::vector<std::size_t> cursors(families.size(), 0);
    while (out->size() < limit) {
        bool added = false;
        for (std::size_t family = 0; family < families.size() && out->size() < limit; family++) {
            if (cursors[family] >= families[family].size()) {
                continue;
            }
            out->push_back(families[family][cursors[family]++]);
            added = true;
        }
        if (!added) {
            break;
        }
    }
}

PrefilterPlan build_prefilter_plan(const Code &code, uint64_t probe_count)
{
    PrefilterPlan plan;
    if (probe_count == 0) {
        return plan;
    }

    int g = static_cast<int>(code.groups.size());
    int h = code.global_parity;
    std::vector<std::vector<StressProbe>> families;
    std::vector<int> all_rows = all_row_ids(h);

    if (code.groups.empty() || code.groups.front().local_parity != 1 || h < 2 || h > 4 || g < 2) {
        families.push_back({});
    } else if (h == 2) {
        std::vector<StressProbe> line;
        for (int left = 0; left < g; left++) {
            for (int right = left + 1; right < g; right++) {
                line.push_back(make_probe(StressProbeKind::LineCollision, {left, right}, {1, 1},
                                          all_rows));
            }
        }

        std::vector<StressProbe> single;
        for (const auto &rows : row_subsets(h, 1)) {
            for (int group = 0; group < g; group++) {
                single.push_back(make_probe(StressProbeKind::SingleRank, {group}, {1}, rows));
            }
        }
        for (int group = 0; group < g; group++) {
            single.push_back(make_probe(StressProbeKind::SingleRank, {group}, {2}, all_rows));
        }
        families = {line, single};
    } else if (h == 3) {
        std::vector<StressProbe> vector_plane;
        for (int plane_group = 0; plane_group < g; plane_group++) {
            for (int vector_group = 0; vector_group < g; vector_group++) {
                if (plane_group == vector_group) {
                    continue;
                }
                vector_plane.push_back(make_probe(StressProbeKind::VectorPlane,
                                                  {plane_group, vector_group}, {2, 1}, all_rows));
            }
        }

        std::vector<StressProbe> three_vectors;
        if (g >= 3) {
            for (int a = 0; a < g; a++) {
                for (int b = a + 1; b < g; b++) {
                    for (int c = b + 1; c < g; c++) {
                        three_vectors.push_back(make_probe(StressProbeKind::GenericRank,
                                                           {a, b, c}, {1, 1, 1}, all_rows));
                    }
                }
            }
        }

        std::vector<StressProbe> projected_lines;
        for (const auto &rows : row_subsets(h, 2)) {
            for (int left = 0; left < g; left++) {
                for (int right = left + 1; right < g; right++) {
                    projected_lines.push_back(make_probe(StressProbeKind::LineCollision,
                                                         {left, right}, {1, 1}, rows));
                }
            }
        }

        std::vector<StressProbe> single;
        for (int extra = 1; extra <= h; extra++) {
            for (const auto &rows : row_subsets(h, extra)) {
                for (int group = 0; group < g; group++) {
                    single.push_back(make_probe(StressProbeKind::SingleRank, {group}, {extra},
                                                rows));
                }
            }
        }
        families = {vector_plane, three_vectors, projected_lines, single};
    } else {
        std::vector<StressProbe> single;
        for (int extra : std::vector<int>{2, 1, 3, 4}) {
            for (const auto &rows : row_subsets(h, extra)) {
                for (int group = 0; group < g; group++) {
                    single.push_back(make_probe(StressProbeKind::SingleRank, {group}, {extra},
                                                rows));
                }
            }
        }

        std::vector<StressProbe> plane_plane;
        for (int left = 0; left < g; left++) {
            for (int right = left + 1; right < g; right++) {
                plane_plane.push_back(make_probe(StressProbeKind::PlanePlane, {left, right},
                                                 {2, 2}, all_rows));
            }
        }

        std::vector<StressProbe> hyperplane_vector;
        for (int hyperplane_group = 0; hyperplane_group < g; hyperplane_group++) {
            for (int vector_group = 0; vector_group < g; vector_group++) {
                if (hyperplane_group == vector_group) {
                    continue;
                }
                hyperplane_vector.push_back(make_probe(StressProbeKind::HyperplaneVector,
                                                       {hyperplane_group, vector_group}, {3, 1},
                                                       all_rows));
            }
        }

        std::vector<StressProbe> plane_vectors;
        if (g >= 3) {
            for (int plane_group = 0; plane_group < g; plane_group++) {
                for (int first_vector = 0; first_vector < g; first_vector++) {
                    if (first_vector == plane_group) {
                        continue;
                    }
                    for (int second_vector = first_vector + 1; second_vector < g; second_vector++) {
                        if (second_vector == plane_group) {
                            continue;
                        }
                        plane_vectors.push_back(make_probe(StressProbeKind::GenericRank,
                                                           {plane_group, first_vector, second_vector},
                                                           {2, 1, 1}, all_rows));
                    }
                }
            }
        }

        std::vector<StressProbe> four_vectors;
        if (g >= 4) {
            for (int a = 0; a < g; a++) {
                for (int b = a + 1; b < g; b++) {
                    for (int c = b + 1; c < g; c++) {
                        for (int d = c + 1; d < g; d++) {
                            four_vectors.push_back(make_probe(StressProbeKind::GenericRank,
                                                              {a, b, c, d}, {1, 1, 1, 1},
                                                              all_rows));
                        }
                    }
                }
            }
        }

        std::vector<StressProbe> projected_vector_plane;
        for (const auto &rows : row_subsets(h, 3)) {
            for (int plane_group = 0; plane_group < g; plane_group++) {
                for (int vector_group = 0; vector_group < g; vector_group++) {
                    if (plane_group == vector_group) {
                        continue;
                    }
                    projected_vector_plane.push_back(make_probe(StressProbeKind::VectorPlane,
                                                                {plane_group, vector_group},
                                                                {2, 1}, rows));
                }
            }
        }

        std::vector<StressProbe> projected_lines;
        for (const auto &rows : row_subsets(h, 2)) {
            for (int left = 0; left < g; left++) {
                for (int right = left + 1; right < g; right++) {
                    projected_lines.push_back(make_probe(StressProbeKind::LineCollision,
                                                         {left, right}, {1, 1}, rows));
                }
            }
        }

        families = {single, plane_plane, hyperplane_vector, plane_vectors, four_vectors,
                    projected_vector_plane, projected_lines};
    }

    std::vector<StressProbe> targeted;
    add_probe_family_round_robin(families, probe_count, &targeted);
    plan.probes = std::move(targeted);

    for (uint64_t sample = 0; plan.probes.size() < probe_count; sample++) {
        StressProbe probe;
        probe.kind = StressProbeKind::RandomResidual;
        probe.random_sample = sample;
        plan.probes.push_back(std::move(probe));
    }

    plan.erasure_choices.assign(
        code.groups.size(),
        std::vector<std::vector<std::vector<int>>>(
            static_cast<std::size_t>(code.global_parity + 1)));
    for (int group_index = 0; group_index < static_cast<int>(code.groups.size()); group_index++) {
        const auto &group = code.groups[static_cast<std::size_t>(group_index)];
        std::vector<int> symbols = group_symbol_ids(code, group_index);
        int max_extra = std::min(code.global_parity,
                                 static_cast<int>(symbols.size()) - group.local_parity);
        for (int extra = 0; extra <= max_extra; extra++) {
            int need = group.local_parity + extra;
            std::vector<int> chosen;
            enumerate_combinations(symbols, need, 0, &chosen,
                                   [&](const std::vector<int> &erased_symbols) {
                                       plan.erasure_choices[static_cast<std::size_t>(group_index)]
                                                           [static_cast<std::size_t>(extra)]
                                                               .push_back(erased_symbols);
                                       return true;
                                   });
        }
    }
    return plan;
}

class RandomStressPrefilter {
public:
    RandomStressPrefilter(const Code &code, const PrefilterPlan &plan)
        : code_(code), plan_(plan),
          erased_(static_cast<std::size_t>(code.symbols), 0), extra_counts_(code.groups.size(), 0),
          extra_caps_(code.groups.size(), 0), group_erased_(code.groups.size()),
          group_symbols_(code.groups.size()), global_rows_(static_cast<std::size_t>(code.global_parity)),
          global_erased_(static_cast<std::size_t>(code.global_parity), 0),
          keep_global_mask_(static_cast<std::size_t>(code.global_parity), 0),
          group_order_(code.groups.size())
    {
        targeted_blocks_.assign(code_.groups.size(),
                                std::vector<std::vector<ResidualBlock>>(
                                    static_cast<std::size_t>(code_.global_parity + 1)));
        targeted_extra_built_.assign(
            code_.groups.size(),
            std::vector<uint8_t>(static_cast<std::size_t>(code_.global_parity + 1), 0));
        baseline_erased_.assign(code_.groups.size(), {});

        all_local_parity_one_ = true;
        std::size_t max_subset_size = static_cast<std::size_t>(code_.global_parity);
        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
            group_symbols_[static_cast<std::size_t>(group_index)] =
                group_symbol_ids(code_, group_index);
            int symbols = static_cast<int>(group_symbols_[static_cast<std::size_t>(group_index)].size());
            extra_caps_[static_cast<std::size_t>(group_index)] =
                std::max(0, symbols - group.local_parity);
            group_erased_[static_cast<std::size_t>(group_index)].reserve(
                static_cast<std::size_t>(symbols));
            baseline_erased_[static_cast<std::size_t>(group_index)].reserve(
                static_cast<std::size_t>(group.local_parity));
            group_order_[static_cast<std::size_t>(group_index)] = group_index;
            max_subset_size = std::max(max_subset_size, static_cast<std::size_t>(symbols));

            if (group.local_parity == 1) {
                baseline_erased_[static_cast<std::size_t>(group_index)].push_back(
                    group.local_row_start);
            } else {
                all_local_parity_one_ = false;
            }
        }

        for (int global = 0; global < code_.global_parity; global++) {
            global_rows_[static_cast<std::size_t>(global)] = global;
        }
        available_groups_.reserve(code_.groups.size());
        subset_scratch_.reserve(max_subset_size);
        erased_global_rows_.reserve(static_cast<std::size_t>(code_.global_parity));
        kept_global_rows_.reserve(static_cast<std::size_t>(code_.global_parity));
        for (auto &blocks : projected_blocks_) {
            blocks.reserve(static_cast<std::size_t>(kTargetedBlockLimit));
        }
        line_index_.reserve(static_cast<std::size_t>(kTargetedBlockLimit));

        targeted_layout_static_ =
            all_local_parity_one_ && code_.global_parity <= 4 && !plan_.erasure_choices.empty();
        if (targeted_layout_static_) {
            initialize_static_targeted_layout();
        }
    }

    void reset(Rng &rng, CheckResult *result, bool capture_failure)
    {
        rng_ = &rng;
        result_ = result;
        capture_failure_ = capture_failure;
        targeted_blocks_initialized_ = targeted_layout_static_;
        std::fill(erased_.begin(), erased_.end(), 0);
        std::fill(extra_counts_.begin(), extra_counts_.end(), 0);
        for (auto &erased : group_erased_) {
            erased.clear();
        }
        kept_global_rows_.clear();
        if (!targeted_layout_static_) {
            for (auto &by_extra : targeted_blocks_) {
                for (auto &blocks : by_extra) {
                    blocks.clear();
                }
            }
        }
        for (auto &built : targeted_extra_built_) {
            std::fill(built.begin(), built.end(), 0);
        }
        if (!targeted_layout_static_) {
            for (auto &baseline : baseline_erased_) {
                baseline.clear();
            }
        }
    }

    void run()
    {
        int max_extra_used = std::min(code_.global_parity, total_extra_capacity());
        for (const auto &probe : plan_.probes) {
            if (should_stop()) {
                break;
            }
            if (probe.kind == StressProbeKind::RandomResidual) {
                int extra_used = choose_extra_used(probe.random_sample, max_extra_used);
                build_erased_pattern(probe.random_sample, extra_used, max_extra_used);
                check_current_pattern();
                continue;
            }
            execute_targeted_probe(probe);
        }

        finish();
    }

private:
    static constexpr uint64_t kTargetedBlockLimit = 4096;
    using ProbeSelection = std::array<std::pair<int, const ResidualBlock *>, 4>;

    void finish()
    {
        result_->is_mr = result_->failures == 0;
        if (!result_->is_mr && capture_failure_) {
            result_->message = "candidate failed random stress prefilter";
        }
    }

    bool should_stop() const
    {
        return result_->failures > 0;
    }

    template <typename Selected>
    void mark_probe_failure_range(const Selected &selected, const std::vector<int> &kept_rows)
    {
        std::fill(erased_.begin(), erased_.end(), 0);

        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            const ResidualBlock *block = nullptr;
            for (const auto &entry : selected) {
                if (entry.first == group_index) {
                    block = entry.second;
                    break;
                }
            }

            const auto &erased_symbols =
                block == nullptr ? baseline_erased_[static_cast<std::size_t>(group_index)]
                                  : block->erased_symbols;
            mark_erased_symbols(erased_symbols);
        }

        std::fill(keep_global_mask_.begin(), keep_global_mask_.end(), 0);
        for (int row : kept_rows) {
            keep_global_mask_[static_cast<std::size_t>(row)] = 1;
        }

        int first_global = code_.data + code_.local_rows;
        for (int global = 0; global < code_.global_parity; global++) {
            if (keep_global_mask_[static_cast<std::size_t>(global)] == 0) {
                erased_[static_cast<std::size_t>(first_global + global)] = 1;
            }
        }

        result_->failures = 1;
        if (capture_failure_) {
            result_->first_failed_erased = erased_symbols_from_mask(code_, erased_);
        }
    }

    void mark_probe_failure(std::initializer_list<std::pair<int, const ResidualBlock *>> selected,
                            const std::vector<int> &kept_rows)
    {
        mark_probe_failure_range(selected, kept_rows);
    }

    void mark_probe_failure(const ProbeSelection &selected, int selected_count,
                            const std::vector<int> &kept_rows)
    {
        std::fill(erased_.begin(), erased_.end(), 0);

        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            const ResidualBlock *block = nullptr;
            for (int entry_index = 0; entry_index < selected_count; entry_index++) {
                const auto &entry = selected[static_cast<std::size_t>(entry_index)];
                if (entry.first == group_index) {
                    block = entry.second;
                    break;
                }
            }

            const auto &erased_symbols =
                block == nullptr ? baseline_erased_[static_cast<std::size_t>(group_index)]
                                  : block->erased_symbols;
            mark_erased_symbols(erased_symbols);
        }

        std::fill(keep_global_mask_.begin(), keep_global_mask_.end(), 0);
        for (int row : kept_rows) {
            keep_global_mask_[static_cast<std::size_t>(row)] = 1;
        }

        int first_global = code_.data + code_.local_rows;
        for (int global = 0; global < code_.global_parity; global++) {
            if (keep_global_mask_[static_cast<std::size_t>(global)] == 0) {
                erased_[static_cast<std::size_t>(first_global + global)] = 1;
            }
        }

        result_->failures = 1;
        if (capture_failure_) {
            result_->first_failed_erased = erased_symbols_from_mask(code_, erased_);
        }
    }

    struct ProjectedProbeBlock {
        const ResidualBlock *block = nullptr;
        std::array<uint8_t, 16> matrix{};
        int rows = 0;
        int cols = 0;
    };

    bool project_block_for_probe_small(int group_index, const ResidualBlock &block,
                                       const std::vector<int> &kept_rows,
                                       ProjectedProbeBlock *projected)
    {
        (void)group_index;
        if (projected == nullptr || (block.matrix.empty() && !block.has_small_matrix)) {
            return false;
        }

        int rows = static_cast<int>(kept_rows.size());
        int cols = block.extra;
        if (rows > 4 || cols > 4) {
            return false;
        }

        ProjectedProbeBlock result;
        result.block = &block;
        result.rows = rows;
        result.cols = cols;
        for (int row = 0; row < rows; row++) {
            int source_row = kept_rows[static_cast<std::size_t>(row)];
            for (int col = 0; col < cols; col++) {
                result.matrix[static_cast<std::size_t>(idx(row, col, cols))] =
                    block.has_small_matrix
                        ? block.small_matrix[static_cast<std::size_t>(
                              idx(source_row, col, cols))]
                        : block.matrix[static_cast<std::size_t>(idx(source_row, col, cols))];
            }
        }
        if (gf256_rank_small(result.matrix, rows, cols) < cols) {
            return false;
        }

        *projected = result;
        return true;
    }

    bool collect_projected_blocks(int group, int extra, const std::vector<int> &kept_rows,
                                  uint64_t scan_limit, uint64_t *scanned,
                                  std::vector<ProjectedProbeBlock> *projected)
    {
        if (scanned == nullptr || projected == nullptr) {
            return false;
        }

        projected->clear();
        for (const auto &block : probe_blocks(group, extra)) {
            if (*scanned >= scan_limit) {
                break;
            }
            (*scanned)++;

            ProjectedProbeBlock matrix;
            if (!project_block_for_probe_small(group, block, kept_rows, &matrix)) {
                result_->patterns_checked++;
                mark_probe_failure({{group, &block}}, kept_rows);
                return false;
            }
            projected->push_back(matrix);
        }
        return !should_stop();
    }

    bool confirm_projected_pair(int left_group, const ProjectedProbeBlock &left,
                                int right_group, const ProjectedProbeBlock &right,
                                const std::vector<int> &kept_rows)
    {
        if (left.block == nullptr || right.block == nullptr) {
            return false;
        }

        int rows = static_cast<int>(kept_rows.size());
        int left_extra = left.block->extra;
        int right_extra = right.block->extra;
        if (left.rows != rows || right.rows != rows || left.cols != left_extra ||
            right.cols != right_extra || left_extra + right_extra != rows) {
            return false;
        }

        int cols = left_extra + right_extra;
        if (rows == 0 || cols != rows) {
            return false;
        }

        std::array<uint8_t, 16> joined{};
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < left_extra; col++) {
                joined[static_cast<std::size_t>(idx(row, col, cols))] =
                    left.matrix[static_cast<std::size_t>(idx(row, col, left_extra))];
            }
            for (int col = 0; col < right_extra; col++) {
                joined[static_cast<std::size_t>(idx(row, left_extra + col, cols))] =
                    right.matrix[static_cast<std::size_t>(idx(row, col, right_extra))];
            }
        }
        if (gf256_rank_small(joined, rows, cols) == rows) {
            return false;
        }

        result_->patterns_checked++;
        mark_probe_failure({{left_group, left.block}, {right_group, right.block}}, kept_rows);
        return true;
    }

    bool confirm_projected_selection(
        const StressProbe &probe, const std::array<const ProjectedProbeBlock *, 4> &selected)
    {
        int rows = static_cast<int>(probe.kept_rows.size());
        if (rows <= 0 || rows > 4) {
            return false;
        }

        int used = 0;
        std::array<uint8_t, 16> joined{};
        ProbeSelection witness{};
        for (int term = 0; term < probe.terms; term++) {
            const ProjectedProbeBlock *block = selected[static_cast<std::size_t>(term)];
            if (block == nullptr || block->block == nullptr || block->rows != rows) {
                return false;
            }
            int extra = block->cols;
            int group = probe.groups[static_cast<std::size_t>(term)];
            witness[static_cast<std::size_t>(term)] = {group, block->block};
            for (int row = 0; row < rows; row++) {
                for (int col = 0; col < extra; col++) {
                    joined[static_cast<std::size_t>(idx(row, used + col, rows))] =
                        block->matrix[static_cast<std::size_t>(idx(row, col, extra))];
                }
            }
            used += extra;
        }
        if (used != rows) {
            return false;
        }
        if (gf256_rank_small(joined, rows, rows) == rows) {
            return false;
        }

        result_->patterns_checked++;
        mark_probe_failure(witness, probe.terms, probe.kept_rows);
        return true;
    }

    bool confirm_probe_blocks(const ProbeSelection &selected, int selected_count,
                              const std::vector<int> &kept_rows)
    {
        int rows = static_cast<int>(kept_rows.size());
        if (rows <= 0 || rows > 4) {
            return false;
        }

        std::array<uint8_t, 16> joined{};
        int used = 0;

        for (int entry_index = 0; entry_index < selected_count; entry_index++) {
            const auto &entry = selected[static_cast<std::size_t>(entry_index)];
            if (entry.second == nullptr) {
                return false;
            }

            ProjectedProbeBlock projected;
            if (!project_block_for_probe_small(entry.first, *entry.second, kept_rows, &projected)) {
                result_->patterns_checked++;
                mark_probe_failure(selected, selected_count, kept_rows);
                return true;
            }

            int extra = entry.second->extra;
            if (used + extra > rows) {
                return false;
            }
            for (int row = 0; row < rows; row++) {
                for (int col = 0; col < extra; col++) {
                    joined[static_cast<std::size_t>(idx(row, used + col, rows))] =
                        projected.matrix[static_cast<std::size_t>(idx(row, col, extra))];
                }
            }
            used += extra;
        }

        if (used != rows) {
            return false;
        }

        if (gf256_rank_small(joined, rows, rows) == rows) {
            return false;
        }

        result_->patterns_checked++;
        mark_probe_failure(selected, selected_count, kept_rows);
        return true;
    }

    const std::vector<ResidualBlock> &probe_blocks(int group, int extra)
    {
        static const std::vector<ResidualBlock> empty;
        if (!build_targeted_blocks(group, extra)) {
            return empty;
        }
        if (group < 0 || group >= static_cast<int>(targeted_blocks_.size())) {
            return empty;
        }
        const auto &by_extra = targeted_blocks_[static_cast<std::size_t>(group)];
        if (extra < 0 || extra >= static_cast<int>(by_extra.size())) {
            return empty;
        }
        return by_extra[static_cast<std::size_t>(extra)];
    }

    bool execute_single_rank_probe(const StressProbe &probe)
    {
        int group = probe.groups[0];
        int extra = probe.extras[0];
        uint64_t scanned = 0;
        for (const auto &block : probe_blocks(group, extra)) {
            if (scanned++ >= probe.scan_limit) {
                break;
            }
            ProjectedProbeBlock projected;
            if (!project_block_for_probe_small(group, block, probe.kept_rows, &projected)) {
                result_->patterns_checked++;
                mark_probe_failure({{group, &block}}, probe.kept_rows);
                return true;
            }
        }
        return false;
    }

    bool execute_line_collision_probe(const StressProbe &probe)
    {
        int left_group = probe.groups[0];
        int right_group = probe.groups[1];
        uint64_t scanned = 0;
        auto &left_blocks = projected_blocks_[0];
        auto &right_blocks = projected_blocks_[1];
        if (!collect_projected_blocks(left_group, 1, probe.kept_rows, probe.scan_limit,
                                      &scanned, &left_blocks) ||
            !collect_projected_blocks(right_group, 1, probe.kept_rows, probe.scan_limit,
                                      &scanned, &right_blocks)) {
            return should_stop();
        }

        line_index_.clear();
        line_index_.reserve(left_blocks.size());
        for (const auto &left : left_blocks) {
            uint32_t key = 0;
            if (projective_column_key_packed(left.matrix, left.rows, left.cols, 0, &key)) {
                line_index_.emplace(key, &left);
            }
        }

        for (const auto &right : right_blocks) {
            uint32_t key = 0;
            if (!projective_column_key_packed(right.matrix, right.rows, right.cols, 0, &key)) {
                continue;
            }
            auto found = line_index_.find(key);
            if (found != line_index_.end()) {
                return confirm_projected_pair(left_group, *found->second, right_group, right,
                                              probe.kept_rows);
            }
        }
        return false;
    }

    bool execute_vector_plane_direction(const StressProbe &probe, int plane_group, int vector_group)
    {
        uint64_t scanned = 0;
        auto &planes = projected_blocks_[0];
        auto &vectors = projected_blocks_[1];
        if (!collect_projected_blocks(plane_group, 2, probe.kept_rows, probe.scan_limit,
                                      &scanned, &planes) ||
            !collect_projected_blocks(vector_group, 1, probe.kept_rows, probe.scan_limit,
                                      &scanned, &vectors)) {
            return should_stop();
        }

        uint64_t pair_scanned = 0;
        for (const auto &plane : planes) {
            for (const auto &vector : vectors) {
                if (pair_scanned++ >= probe.scan_limit) {
                    return false;
                }
                if (confirm_projected_pair(plane_group, plane, vector_group, vector,
                                           probe.kept_rows)) {
                    return true;
                }
            }
        }
        return false;
    }

    bool execute_vector_plane_probe(const StressProbe &probe)
    {
        int plane_group = probe.groups[0];
        int vector_group = probe.groups[1];
        return execute_vector_plane_direction(probe, plane_group, vector_group);
    }

    bool execute_plane_plane_probe(const StressProbe &probe)
    {
        int left_group = probe.groups[0];
        int right_group = probe.groups[1];
        uint64_t scanned = 0;
        auto &left_blocks = projected_blocks_[0];
        auto &right_blocks = projected_blocks_[1];
        if (!collect_projected_blocks(left_group, 2, probe.kept_rows, probe.scan_limit,
                                      &scanned, &left_blocks) ||
            !collect_projected_blocks(right_group, 2, probe.kept_rows, probe.scan_limit,
                                      &scanned, &right_blocks)) {
            return should_stop();
        }

        uint64_t pair_scanned = 0;
        for (const auto &left : left_blocks) {
            for (const auto &right : right_blocks) {
                if (pair_scanned++ >= probe.scan_limit) {
                    return false;
                }
                if (confirm_projected_pair(left_group, left, right_group, right,
                                           probe.kept_rows)) {
                    return true;
                }
            }
        }
        return false;
    }

    bool execute_hyperplane_vector_direction(const StressProbe &probe, int hyperplane_group,
                                             int vector_group)
    {
        uint64_t scanned = 0;
        auto &hyperplanes = projected_blocks_[0];
        auto &vectors = projected_blocks_[1];
        if (!collect_projected_blocks(hyperplane_group, 3, probe.kept_rows, probe.scan_limit,
                                      &scanned, &hyperplanes) ||
            !collect_projected_blocks(vector_group, 1, probe.kept_rows, probe.scan_limit,
                                      &scanned, &vectors)) {
            return should_stop();
        }

        uint64_t pair_scanned = 0;
        for (const auto &hyperplane : hyperplanes) {
            for (const auto &vector : vectors) {
                if (pair_scanned++ >= probe.scan_limit) {
                    return false;
                }
                if (confirm_projected_pair(hyperplane_group, hyperplane, vector_group, vector,
                                           probe.kept_rows)) {
                    return true;
                }
            }
        }
        return false;
    }

    bool execute_hyperplane_vector_probe(const StressProbe &probe)
    {
        int hyperplane_group = probe.groups[0];
        int vector_group = probe.groups[1];
        return execute_hyperplane_vector_direction(probe, hyperplane_group, vector_group);
    }

    bool execute_generic_rank_probe_recursive(
        const StressProbe &probe, int term, uint64_t *scanned,
        ProbeSelection *selected)
    {
        if (scanned == nullptr || selected == nullptr || *scanned >= probe.scan_limit) {
            return false;
        }
        if (term == probe.terms) {
            (*scanned)++;
            return confirm_probe_blocks(*selected, probe.terms, probe.kept_rows);
        }
        if (term < 0 || term >= 4) {
            return false;
        }

        int group = probe.groups[static_cast<std::size_t>(term)];
        int extra = probe.extras[static_cast<std::size_t>(term)];
        for (const auto &block : probe_blocks(group, extra)) {
            if (*scanned >= probe.scan_limit) {
                break;
            }
            (*selected)[static_cast<std::size_t>(term)] = {group, &block};
            if (execute_generic_rank_probe_recursive(probe, term + 1, scanned, selected)) {
                return true;
            }
            (*selected)[static_cast<std::size_t>(term)] = {-1, nullptr};
        }
        return false;
    }

    static uint8_t projected_vector_entry(const ProjectedProbeBlock &block, int row)
    {
        return block.matrix[static_cast<std::size_t>(idx(row, 0, 1))];
    }

    static uint8_t det3(uint8_t a00, uint8_t a01, uint8_t a02,
                        uint8_t a10, uint8_t a11, uint8_t a12,
                        uint8_t a20, uint8_t a21, uint8_t a22)
    {
        return gf256_mul(a00, gf256_mul(a11, a22)) ^
               gf256_mul(a00, gf256_mul(a12, a21)) ^
               gf256_mul(a01, gf256_mul(a10, a22)) ^
               gf256_mul(a01, gf256_mul(a12, a20)) ^
               gf256_mul(a02, gf256_mul(a10, a21)) ^
               gf256_mul(a02, gf256_mul(a11, a20));
    }

    static std::array<uint8_t, 4> four_dimensional_normal(const ProjectedProbeBlock &a,
                                                          const ProjectedProbeBlock &b,
                                                          const ProjectedProbeBlock &c)
    {
        std::array<uint8_t, 4> av{
            projected_vector_entry(a, 0), projected_vector_entry(a, 1),
            projected_vector_entry(a, 2), projected_vector_entry(a, 3)};
        std::array<uint8_t, 4> bv{
            projected_vector_entry(b, 0), projected_vector_entry(b, 1),
            projected_vector_entry(b, 2), projected_vector_entry(b, 3)};
        std::array<uint8_t, 4> cv{
            projected_vector_entry(c, 0), projected_vector_entry(c, 1),
            projected_vector_entry(c, 2), projected_vector_entry(c, 3)};

        return {
            det3(av[1], av[2], av[3], bv[1], bv[2], bv[3], cv[1], cv[2], cv[3]),
            det3(av[0], av[2], av[3], bv[0], bv[2], bv[3], cv[0], cv[2], cv[3]),
            det3(av[0], av[1], av[3], bv[0], bv[1], bv[3], cv[0], cv[1], cv[3]),
            det3(av[0], av[1], av[2], bv[0], bv[1], bv[2], cv[0], cv[1], cv[2]),
        };
    }

    bool execute_three_vector_probe(const StressProbe &probe)
    {
        if (probe.terms != 3 || static_cast<int>(probe.kept_rows.size()) != 3) {
            return false;
        }
        for (int term = 0; term < 3; term++) {
            if (probe.extras[static_cast<std::size_t>(term)] != 1) {
                return false;
            }
        }

        uint64_t scanned = 0;
        for (int term = 0; term < 3; term++) {
            int group = probe.groups[static_cast<std::size_t>(term)];
            if (!collect_projected_blocks(group, 1, probe.kept_rows, probe.scan_limit,
                                          &scanned,
                                          &projected_blocks_[static_cast<std::size_t>(term)])) {
                return should_stop();
            }
        }
        if (projected_blocks_[0].empty() || projected_blocks_[1].empty() ||
            projected_blocks_[2].empty()) {
            return false;
        }

        uint64_t combinations_scanned = 0;
        std::array<const ProjectedProbeBlock *, 4> selected{};
        for (const auto &left : projected_blocks_[0]) {
            selected[0] = &left;
            uint8_t a0 = projected_vector_entry(left, 0);
            uint8_t a1 = projected_vector_entry(left, 1);
            uint8_t a2 = projected_vector_entry(left, 2);
            for (const auto &middle : projected_blocks_[1]) {
                selected[1] = &middle;
                uint8_t b0 = projected_vector_entry(middle, 0);
                uint8_t b1 = projected_vector_entry(middle, 1);
                uint8_t b2 = projected_vector_entry(middle, 2);
                std::array<uint8_t, 3> normal{
                    static_cast<uint8_t>(gf256_mul(a1, b2) ^ gf256_mul(a2, b1)),
                    static_cast<uint8_t>(gf256_mul(a0, b2) ^ gf256_mul(a2, b0)),
                    static_cast<uint8_t>(gf256_mul(a0, b1) ^ gf256_mul(a1, b0)),
                };

                for (const auto &right : projected_blocks_[2]) {
                    if (combinations_scanned++ >= probe.scan_limit) {
                        return false;
                    }
                    selected[2] = &right;
                    uint8_t dot =
                        gf256_mul(normal[0], projected_vector_entry(right, 0)) ^
                        gf256_mul(normal[1], projected_vector_entry(right, 1)) ^
                        gf256_mul(normal[2], projected_vector_entry(right, 2));
                    if (dot == 0 && confirm_projected_selection(probe, selected)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    bool execute_four_vector_probe(const StressProbe &probe)
    {
        if (probe.terms != 4 || static_cast<int>(probe.kept_rows.size()) != 4) {
            return false;
        }
        for (int term = 0; term < 4; term++) {
            if (probe.extras[static_cast<std::size_t>(term)] != 1) {
                return false;
            }
        }

        uint64_t scanned = 0;
        for (int term = 0; term < 4; term++) {
            int group = probe.groups[static_cast<std::size_t>(term)];
            if (!collect_projected_blocks(group, 1, probe.kept_rows, probe.scan_limit,
                                          &scanned,
                                          &projected_blocks_[static_cast<std::size_t>(term)])) {
                return should_stop();
            }
        }
        for (int term = 0; term < 4; term++) {
            if (projected_blocks_[static_cast<std::size_t>(term)].empty()) {
                return false;
            }
        }

        uint64_t combinations_scanned = 0;
        std::array<const ProjectedProbeBlock *, 4> selected{};
        for (const auto &first : projected_blocks_[0]) {
            selected[0] = &first;
            for (const auto &second : projected_blocks_[1]) {
                selected[1] = &second;
                for (const auto &third : projected_blocks_[2]) {
                    selected[2] = &third;
                    std::array<uint8_t, 4> normal =
                        four_dimensional_normal(first, second, third);
                    for (const auto &fourth : projected_blocks_[3]) {
                        if (combinations_scanned++ >= probe.scan_limit) {
                            return false;
                        }
                        selected[3] = &fourth;
                        uint8_t dot =
                            gf256_mul(normal[0], projected_vector_entry(fourth, 0)) ^
                            gf256_mul(normal[1], projected_vector_entry(fourth, 1)) ^
                            gf256_mul(normal[2], projected_vector_entry(fourth, 2)) ^
                            gf256_mul(normal[3], projected_vector_entry(fourth, 3));
                        if (dot == 0 && confirm_projected_selection(probe, selected)) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    bool search_projected_selection(const StressProbe &probe, int term,
                                    std::array<const ProjectedProbeBlock *, 4> *selected,
                                    uint64_t *combinations_scanned)
    {
        if (selected == nullptr || combinations_scanned == nullptr ||
            *combinations_scanned >= probe.scan_limit) {
            return false;
        }
        if (term == probe.terms) {
            (*combinations_scanned)++;
            return confirm_projected_selection(probe, *selected);
        }

        for (const auto &block : projected_blocks_[static_cast<std::size_t>(term)]) {
            (*selected)[static_cast<std::size_t>(term)] = &block;
            if (search_projected_selection(probe, term + 1, selected, combinations_scanned)) {
                return true;
            }
            if (*combinations_scanned >= probe.scan_limit) {
                return false;
            }
        }
        (*selected)[static_cast<std::size_t>(term)] = nullptr;
        return false;
    }

    bool execute_generic_rank_probe(const StressProbe &probe)
    {
        bool three_vector_probe =
            probe.terms == 3 && static_cast<int>(probe.kept_rows.size()) == 3;
        for (int term = 0; term < probe.terms && three_vector_probe; term++) {
            three_vector_probe = probe.extras[static_cast<std::size_t>(term)] == 1;
        }
        if (three_vector_probe) {
            return execute_three_vector_probe(probe);
        }
        bool four_vector_probe =
            probe.terms == 4 && static_cast<int>(probe.kept_rows.size()) == 4;
        for (int term = 0; term < probe.terms && four_vector_probe; term++) {
            four_vector_probe = probe.extras[static_cast<std::size_t>(term)] == 1;
        }
        if (four_vector_probe) {
            return execute_four_vector_probe(probe);
        }

        int rows = static_cast<int>(probe.kept_rows.size());
        int total_extra = 0;
        for (int term = 0; term < probe.terms; term++) {
            total_extra += probe.extras[static_cast<std::size_t>(term)];
        }
        if (probe.terms > 0 && probe.terms <= 4 && rows > 0 && rows <= 4 &&
            total_extra == rows) {
            uint64_t scanned = 0;
            for (int term = 0; term < probe.terms; term++) {
                int group = probe.groups[static_cast<std::size_t>(term)];
                int extra = probe.extras[static_cast<std::size_t>(term)];
                if (!collect_projected_blocks(group, extra, probe.kept_rows, probe.scan_limit,
                                              &scanned,
                                              &projected_blocks_[static_cast<std::size_t>(term)])) {
                    return should_stop();
                }
            }

            uint64_t combinations_scanned = 0;
            std::array<const ProjectedProbeBlock *, 4> selected{};
            return search_projected_selection(probe, 0, &selected, &combinations_scanned);
        }

        uint64_t scanned = 0;
        ProbeSelection selected{};
        return execute_generic_rank_probe_recursive(probe, 0, &scanned, &selected);
    }

    void execute_targeted_probe(const StressProbe &probe)
    {
        switch (probe.kind) {
        case StressProbeKind::SingleRank:
            (void)execute_single_rank_probe(probe);
            return;
        case StressProbeKind::LineCollision:
            (void)execute_line_collision_probe(probe);
            return;
        case StressProbeKind::VectorPlane:
            (void)execute_vector_plane_probe(probe);
            return;
        case StressProbeKind::PlanePlane:
            (void)execute_plane_plane_probe(probe);
            return;
        case StressProbeKind::HyperplaneVector:
            (void)execute_hyperplane_vector_probe(probe);
            return;
        case StressProbeKind::GenericRank:
            (void)execute_generic_rank_probe(probe);
            return;
        case StressProbeKind::RandomResidual:
            return;
        }
    }

    void initialize_static_targeted_layout()
    {
        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            for (int extra = 1; extra <= code_.global_parity; extra++) {
                const auto &choices =
                    plan_.erasure_choices[static_cast<std::size_t>(group_index)]
                                         [static_cast<std::size_t>(extra)];
                if (choices.empty() || choices.size() > kTargetedBlockLimit) {
                    continue;
                }

                auto &blocks = targeted_blocks_[static_cast<std::size_t>(group_index)]
                                              [static_cast<std::size_t>(extra)];
                blocks.reserve(choices.size());
                for (const auto &erased_symbols : choices) {
                    ResidualBlock block;
                    block.extra = extra;
                    block.erased_symbols = erased_symbols;
                    block.small_rows = code_.global_parity;
                    block.small_cols = extra;
                    block.has_small_matrix = true;
                    blocks.push_back(std::move(block));
                }
            }
        }
    }

    bool initialize_targeted_storage()
    {
        if (targeted_blocks_initialized_) {
            return !targeted_blocks_.empty();
        }
        targeted_blocks_initialized_ = true;

        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
            if (group.local_parity == 1) {
                if (baseline_erased_[static_cast<std::size_t>(group_index)].empty()) {
                    baseline_erased_[static_cast<std::size_t>(group_index)].push_back(
                        group.local_row_start);
                }
                continue;
            }

            const auto &choices = plan_.erasure_choices[static_cast<std::size_t>(group_index)][0];
            for (const auto &erased_symbols : choices) {
                std::vector<uint8_t> local =
                    build_local_erasure_submatrix(code_, group_index, erased_symbols);
                if (gf256_rank(local, group.local_parity,
                               static_cast<int>(erased_symbols.size())) == group.local_parity) {
                    baseline_erased_[static_cast<std::size_t>(group_index)] = erased_symbols;
                    break;
                }
            }
            if (baseline_erased_[static_cast<std::size_t>(group_index)].empty()) {
                const auto &symbols = group_symbols_[static_cast<std::size_t>(group_index)];
                auto &baseline = baseline_erased_[static_cast<std::size_t>(group_index)];
                for (int i = 0; i < group.local_parity && i < static_cast<int>(symbols.size()); i++) {
                    baseline.push_back(symbols[static_cast<std::size_t>(i)]);
                }
            }
        }

        return true;
    }

    bool build_targeted_blocks(int group_index, int extra)
    {
        if (!initialize_targeted_storage()) {
            return false;
        }
        if (group_index < 0 || group_index >= static_cast<int>(targeted_blocks_.size()) ||
            extra <= 0 || extra > code_.global_parity) {
            return false;
        }

        auto &built = targeted_extra_built_[static_cast<std::size_t>(group_index)]
                                           [static_cast<std::size_t>(extra)];
        if (built != 0) {
            return true;
        }
        built = 1;

        const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
        const auto &choices = plan_.erasure_choices[static_cast<std::size_t>(group_index)]
                                                   [static_cast<std::size_t>(extra)];
        if (choices.empty()) {
            return true;
        }

        if (choices.size() > kTargetedBlockLimit) {
            return false;
        }
        auto &blocks = targeted_blocks_[static_cast<std::size_t>(group_index)]
                                      [static_cast<std::size_t>(extra)];
        if (!targeted_layout_static_) {
            blocks.reserve(choices.size());
        }
        bool local_rank_is_automatic = group.local_parity == 1;
        bool use_small_matrix = local_rank_is_automatic && code_.global_parity <= 4 && extra <= 4;
        if (targeted_layout_static_) {
            for (auto &block : blocks) {
                build_residual_erasure_block_a1_small(code_, group_index, block.erased_symbols,
                                                      extra, &block.small_matrix);
            }
            return true;
        }

        for (const auto &erased_symbols : choices) {
            if (!local_rank_is_automatic) {
                std::vector<uint8_t> local =
                    build_local_erasure_submatrix(code_, group_index, erased_symbols);
                if (gf256_rank(local, group.local_parity,
                               static_cast<int>(erased_symbols.size())) <
                    group.local_parity) {
                    ResidualBlock block;
                    block.extra = extra;
                    block.erased_symbols = erased_symbols;
                    blocks.push_back(std::move(block));
                    continue;
                }
            }

            ResidualBlock block;
            block.extra = extra;
            block.erased_symbols = erased_symbols;
            if (use_small_matrix) {
                block.small_rows = code_.global_parity;
                block.small_cols = extra;
                block.has_small_matrix = true;
                build_residual_erasure_block_a1_small(code_, group_index, erased_symbols,
                                                      extra, &block.small_matrix);
            } else {
                block.matrix = build_residual_erasure_block(code_, group_index,
                                                            erased_symbols, extra);
            }
            blocks.push_back(std::move(block));
        }

        return true;
    }


    int total_extra_capacity() const
    {
        int total = 0;
        for (int cap : extra_caps_) {
            total += cap;
        }
        return total;
    }

    int choose_extra_used(uint64_t sample, int max_extra_used) const
    {
        if (max_extra_used <= 0) {
            return 0;
        }
        return 1 + static_cast<int>(sample % static_cast<uint64_t>(max_extra_used));
    }

    void shuffle_group_order()
    {
        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            group_order_[static_cast<std::size_t>(group_index)] = group_index;
        }
        shuffle(*rng_, &group_order_);
    }

    bool fill_balanced_extra_counts(int extra_used)
    {
        std::fill(extra_counts_.begin(), extra_counts_.end(), 0);
        shuffle_group_order();
        int remaining = extra_used;
        while (remaining > 0) {
            bool placed = false;
            for (int group_index : group_order_) {
                auto pos = static_cast<std::size_t>(group_index);
                if (extra_counts_[pos] >= extra_caps_[pos]) {
                    continue;
                }
                extra_counts_[pos]++;
                remaining--;
                placed = true;
                if (remaining == 0) {
                    break;
                }
            }
            if (!placed) {
                return false;
            }
        }
        return true;
    }

    bool fill_random_extra_counts(int extra_used)
    {
        std::fill(extra_counts_.begin(), extra_counts_.end(), 0);
        for (int used = 0; used < extra_used; used++) {
            available_groups_.clear();
            for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
                auto pos = static_cast<std::size_t>(group_index);
                if (extra_counts_[pos] < extra_caps_[pos]) {
                    available_groups_.push_back(group_index);
                }
            }
            if (available_groups_.empty()) {
                return false;
            }
            int picked = available_groups_[static_cast<std::size_t>(
                rng_->range(available_groups_.size()))];
            extra_counts_[static_cast<std::size_t>(picked)]++;
        }
        return true;
    }

    void choose_extra_counts(uint64_t sample, int extra_used, int max_extra_used)
    {
        uint64_t sweep = max_extra_used <= 0 ? 0 : sample / static_cast<uint64_t>(max_extra_used);
        bool use_balanced_shape = sweep % 4 != 3;
        bool ok = use_balanced_shape ? fill_balanced_extra_counts(extra_used)
                                     : fill_random_extra_counts(extra_used);
        if (!ok) {
            ok = use_balanced_shape ? fill_random_extra_counts(extra_used)
                                    : fill_balanced_extra_counts(extra_used);
        }
        if (!ok) {
            throw ValidationError("could not build a feasible prefilter erasure shape");
        }
    }

    void choose_random_subset_into(const std::vector<int> &symbols, int need,
                                   std::vector<int> *chosen)
    {
        if (chosen == nullptr) {
            return;
        }
        chosen->clear();
        if (need <= 0) {
            return;
        }
        if (need > static_cast<int>(symbols.size())) {
            throw ValidationError("prefilter erasure shape asks for too many symbols");
        }

        subset_scratch_.assign(symbols.begin(), symbols.end());
        shuffle(*rng_, &subset_scratch_);
        for (int i = 0; i < need; i++) {
            chosen->push_back(subset_scratch_[static_cast<std::size_t>(i)]);
        }
    }

    void mark_erased_symbols(const std::vector<int> &symbols)
    {
        for (int symbol : symbols) {
            erased_[static_cast<std::size_t>(symbol)] = 1;
        }
    }

    void build_erased_pattern(uint64_t sample, int extra_used, int max_extra_used)
    {
        std::fill(erased_.begin(), erased_.end(), 0);
        for (auto &erased : group_erased_) {
            erased.clear();
        }
        kept_global_rows_.clear();
        choose_extra_counts(sample, extra_used, max_extra_used);

        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
            int need = group.local_parity + extra_counts_[static_cast<std::size_t>(group_index)];
            choose_random_subset_into(group_symbols_[static_cast<std::size_t>(group_index)], need,
                                      &group_erased_[static_cast<std::size_t>(group_index)]);
            mark_erased_symbols(group_erased_[static_cast<std::size_t>(group_index)]);
        }

        int first_global = code_.data + code_.local_rows;
        choose_random_subset_into(global_rows_, code_.global_parity - extra_used,
                                  &erased_global_rows_);
        std::fill(global_erased_.begin(), global_erased_.end(), 0);
        for (int global : erased_global_rows_) {
            global_erased_[static_cast<std::size_t>(global)] = 1;
            erased_[static_cast<std::size_t>(first_global + global)] = 1;
        }
        for (int global = 0; global < code_.global_parity; global++) {
            if (global_erased_[static_cast<std::size_t>(global)] == 0) {
                kept_global_rows_.push_back(global);
            }
        }
    }

    bool current_pattern_is_recoverable_by_residual_small_a1()
    {
        int residual_size = static_cast<int>(kept_global_rows_.size());
        std::array<uint8_t, 16> joined{};
        int used = 0;

        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            const auto &erased_symbols = group_erased_[static_cast<std::size_t>(group_index)];
            int extra = static_cast<int>(erased_symbols.size()) - 1;
            if (extra < 0 || extra > 4) {
                return false;
            }
            if (extra == 0) {
                continue;
            }

            std::array<uint8_t, 16> block{};
            build_residual_erasure_block_a1_small(code_, group_index, erased_symbols, extra, &block);

            std::array<uint8_t, 16> projected{};
            for (int row = 0; row < residual_size; row++) {
                int source_row = kept_global_rows_[static_cast<std::size_t>(row)];
                for (int col = 0; col < extra; col++) {
                    projected[static_cast<std::size_t>(idx(row, col, extra))] =
                        block[static_cast<std::size_t>(idx(source_row, col, extra))];
                }
            }
            if (gf256_rank_small(projected, residual_size, extra) < extra) {
                return false;
            }
            if (used + extra > residual_size) {
                return false;
            }

            for (int row = 0; row < residual_size; row++) {
                for (int col = 0; col < extra; col++) {
                    joined[static_cast<std::size_t>(idx(row, used + col, residual_size))] =
                        projected[static_cast<std::size_t>(idx(row, col, extra))];
                }
            }
            used += extra;
            if (gf256_rank_small(joined, residual_size, used) < used) {
                return false;
            }
        }

        return used == residual_size &&
               (residual_size == 0 ||
                gf256_rank_small(joined, residual_size, residual_size) == residual_size);
    }

    bool current_pattern_is_recoverable_by_residual()
    {
        if (all_local_parity_one_ && code_.global_parity <= 4) {
            return current_pattern_is_recoverable_by_residual_small_a1();
        }

        int residual_size = static_cast<int>(kept_global_rows_.size());
        std::vector<uint8_t> joined;
        int used = 0;

        for (int group_index = 0; group_index < static_cast<int>(code_.groups.size()); group_index++) {
            const auto &group = code_.groups[static_cast<std::size_t>(group_index)];
            const auto &erased_symbols = group_erased_[static_cast<std::size_t>(group_index)];
            int extra = static_cast<int>(erased_symbols.size()) - group.local_parity;
            if (extra < 0) {
                return false;
            }

            std::vector<uint8_t> local =
                build_local_erasure_submatrix(code_, group_index, erased_symbols);
            if (gf256_rank(local, group.local_parity, static_cast<int>(erased_symbols.size())) <
                group.local_parity) {
                return false;
            }
            if (extra == 0) {
                continue;
            }

            std::vector<uint8_t> block =
                build_residual_erasure_block(code_, group_index, erased_symbols, extra);
            std::vector<uint8_t> projected =
                project_rows(block, code_.global_parity, extra, kept_global_rows_);
            if (gf256_rank(projected, residual_size, extra) < extra) {
                return false;
            }

            joined = append_columns(joined, residual_size, used, projected, extra);
            used += extra;
            if (gf256_rank(joined, residual_size, used) < used) {
                return false;
            }
        }

        return used == residual_size &&
               (residual_size == 0 || gf256_rank(joined, residual_size, residual_size) == residual_size);
    }

    void check_current_pattern()
    {
        if (should_stop()) {
            return;
        }

        result_->patterns_checked++;
        if (!current_pattern_is_recoverable_by_residual()) {
            result_->failures = 1;
            if (capture_failure_) {
                result_->first_failed_erased = erased_symbols_from_mask(code_, erased_);
            }
        }
    }

    const Code &code_;
    const PrefilterPlan &plan_;
    Rng *rng_ = nullptr;
    CheckResult *result_ = nullptr;
    bool capture_failure_ = false;
    std::vector<uint8_t> erased_;
    std::vector<int> extra_counts_;
    std::vector<int> extra_caps_;
    std::vector<std::vector<int>> group_erased_;
    std::vector<int> kept_global_rows_;
    std::vector<std::vector<int>> group_symbols_;
    std::vector<int> global_rows_;
    std::vector<int> erased_global_rows_;
    std::vector<uint8_t> global_erased_;
    std::vector<uint8_t> keep_global_mask_;
    std::vector<int> group_order_;
    std::vector<int> available_groups_;
    std::vector<int> subset_scratch_;
    std::array<std::vector<ProjectedProbeBlock>, 4> projected_blocks_;
    ankerl::unordered_dense::map<uint32_t, const ProjectedProbeBlock *> line_index_;
    bool all_local_parity_one_ = false;
    bool targeted_layout_static_ = false;
    bool targeted_blocks_initialized_ = false;
    std::vector<std::vector<std::vector<ResidualBlock>>> targeted_blocks_;
    std::vector<std::vector<uint8_t>> targeted_extra_built_;
    std::vector<std::vector<int>> baseline_erased_;
};

} // namespace

const char *family_name(MatrixFamily family)
{
    switch (family) {
    case MatrixFamily::Cauchy:
        return "cauchy";
    case MatrixFamily::ColumnMultiplierCauchy:
        return "column_multiplier_cauchy";
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
    if (name == "column_multiplier_cauchy" || name == "column-multiplier-cauchy" ||
        name == "cm_cauchy" || name == "cm-cauchy") {
        *family = MatrixFamily::ColumnMultiplierCauchy;
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
        result.strict_complete = false;
        result.failures = gate.failures;
        result.first_failed_erased = std::move(gate.first_failed_erased);
        result.message = "residual checker passed, but the full maximal erasure gate failed";
        return result;
    }

    result.strict_complete = true;
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
        result.cauchy_dedup_enabled = use_cauchy_dedup(params);
        result.prefilter_enabled = params.prefilter_count != 0;
        if (result.cauchy_dedup_enabled) {
            result.cauchy_dedup_key_bytes = cauchy_dedup_key_bytes(base_code);
        }
        PrefilterPlan prefilter_plan = build_prefilter_plan(base_code, params.prefilter_count);

        struct WorkerStats {
            uint64_t attempts_done = 0;
            uint64_t checked_candidates = 0;
            uint64_t strict_candidates = 0;
            uint64_t total_patterns_checked = 0;
            uint64_t failed_candidates = 0;
            uint64_t duplicate_candidates = 0;
            uint64_t prefilter_candidates = 0;
            uint64_t prefilter_rejected_candidates = 0;
            uint64_t prefilter_patterns_checked = 0;
        };

        struct alignas(64) PublishedWorkerStats {
            std::atomic<uint64_t> attempts_done{0};
            std::atomic<uint64_t> checked_candidates{0};
            std::atomic<uint64_t> strict_candidates{0};
            std::atomic<uint64_t> total_patterns_checked{0};
            std::atomic<uint64_t> failed_candidates{0};
            std::atomic<uint64_t> duplicate_candidates{0};
            std::atomic<uint64_t> prefilter_candidates{0};
            std::atomic<uint64_t> prefilter_rejected_candidates{0};
            std::atomic<uint64_t> prefilter_patterns_checked{0};
        };

        constexpr uint64_t kStatsPublishInterval = 64;
        uint64_t worker_count = std::min(params.thread_count, params.random_limit);
        bool progress_enabled = params.progress_callback && params.step_time != 0;
        std::vector<WorkerStats> final_worker_stats(static_cast<std::size_t>(worker_count));
        std::unique_ptr<PublishedWorkerStats[]> published_stats(
            progress_enabled && worker_count != 0 ? new PublishedWorkerStats[worker_count]
                                                  : nullptr);
        std::atomic<bool> stop{false};
        std::atomic<bool> progress_done{false};
        std::atomic<bool> first_failure_recorded{false};
        std::mutex result_mutex;
        std::mutex dedup_mutex;
        std::mutex progress_mutex;
        std::condition_variable progress_cv;
        std::unique_ptr<CauchyDeduper> cauchy_deduper;
        if (result.cauchy_dedup_enabled) {
            cauchy_deduper = make_cauchy_deduper(result.cauchy_dedup_key_bytes);
        }
        std::string worker_error;
        std::string progress_error;

        auto accumulate_stats = [](WorkerStats *total, const WorkerStats &stats) {
            total->attempts_done += stats.attempts_done;
            total->checked_candidates += stats.checked_candidates;
            total->strict_candidates += stats.strict_candidates;
            total->total_patterns_checked += stats.total_patterns_checked;
            total->failed_candidates += stats.failed_candidates;
            total->duplicate_candidates += stats.duplicate_candidates;
            total->prefilter_candidates += stats.prefilter_candidates;
            total->prefilter_rejected_candidates += stats.prefilter_rejected_candidates;
            total->prefilter_patterns_checked += stats.prefilter_patterns_checked;
        };

        auto publish_stats = [&](std::size_t worker_index, const WorkerStats &stats) {
            if (!progress_enabled) {
                return;
            }
            PublishedWorkerStats &published = published_stats[worker_index];
            published.attempts_done.store(stats.attempts_done, std::memory_order_relaxed);
            published.checked_candidates.store(stats.checked_candidates, std::memory_order_relaxed);
            published.strict_candidates.store(stats.strict_candidates, std::memory_order_relaxed);
            published.total_patterns_checked.store(stats.total_patterns_checked,
                                                   std::memory_order_relaxed);
            published.failed_candidates.store(stats.failed_candidates, std::memory_order_relaxed);
            published.duplicate_candidates.store(stats.duplicate_candidates,
                                                 std::memory_order_relaxed);
            published.prefilter_candidates.store(stats.prefilter_candidates,
                                                 std::memory_order_relaxed);
            published.prefilter_rejected_candidates.store(
                stats.prefilter_rejected_candidates, std::memory_order_relaxed);
            published.prefilter_patterns_checked.store(stats.prefilter_patterns_checked,
                                                       std::memory_order_relaxed);
        };

        auto sum_stats = [&]() {
            WorkerStats total;
            for (uint64_t worker = 0; worker < worker_count; worker++) {
                PublishedWorkerStats &published = published_stats[static_cast<std::size_t>(worker)];
                total.attempts_done +=
                    published.attempts_done.load(std::memory_order_relaxed);
                total.checked_candidates +=
                    published.checked_candidates.load(std::memory_order_relaxed);
                total.strict_candidates +=
                    published.strict_candidates.load(std::memory_order_relaxed);
                total.total_patterns_checked +=
                    published.total_patterns_checked.load(std::memory_order_relaxed);
                total.failed_candidates +=
                    published.failed_candidates.load(std::memory_order_relaxed);
                total.duplicate_candidates +=
                    published.duplicate_candidates.load(std::memory_order_relaxed);
                total.prefilter_candidates +=
                    published.prefilter_candidates.load(std::memory_order_relaxed);
                total.prefilter_rejected_candidates +=
                    published.prefilter_rejected_candidates.load(std::memory_order_relaxed);
                total.prefilter_patterns_checked +=
                    published.prefilter_patterns_checked.load(std::memory_order_relaxed);
            }
            return total;
        };

        auto sum_final_stats = [&]() {
            WorkerStats total;
            for (const WorkerStats &stats : final_worker_stats) {
                accumulate_stats(&total, stats);
            }
            return total;
        };

        auto finish_stats = [&]() {
            WorkerStats total = sum_final_stats();
            result.attempts_done = total.attempts_done;
            result.unique_candidates_checked = total.checked_candidates;
            result.strict_candidates_checked = total.strict_candidates;
            result.duplicate_candidates_skipped = total.duplicate_candidates;
            result.prefilter_candidates_checked = total.prefilter_candidates;
            result.prefilter_candidates_rejected = total.prefilter_rejected_candidates;
            result.prefilter_patterns_checked = total.prefilter_patterns_checked;
            return total;
        };

        std::thread progress_worker;
        if (params.progress_callback && params.step_time != 0) {
            progress_worker = std::thread([&]() {
                std::unique_lock<std::mutex> lock(progress_mutex);
                uint64_t progress_step = 0;
                for (;;) {
                    if (progress_cv.wait_for(lock, std::chrono::seconds(params.step_time),
                                             [&]() { return progress_done.load(std::memory_order_relaxed); })) {
                        break;
                    }

                    progress_step++;
                    WorkerStats progress_stats = sum_stats();
                    uint64_t searched = progress_stats.attempts_done;
                    uint64_t strict_checked = progress_stats.strict_candidates;
                    lock.unlock();
                    try {
                        params.progress_callback(progress_step, searched, strict_checked);
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

        auto run_attempt = [&](uint64_t attempt, WorkerStats *stats, Code *code,
                               RandomStressPrefilter *stress_prefilter) {
            code->attempt = attempt;
            uint64_t attempt_seed = seed_for_attempt(params.seed, attempt);
            Rng rng(attempt_seed);
            CauchyKeyBytes cauchy_key;
            CauchyKeyBytes *cauchy_key_ptr = result.cauchy_dedup_enabled ? &cauchy_key : nullptr;
            build_candidate(code, rng, cauchy_key_ptr);
            stats->attempts_done++;

            if (cauchy_key_ptr != nullptr) {
                std::lock_guard<std::mutex> lock(dedup_mutex);
                if (!cauchy_deduper->insert(cauchy_key)) {
                    stats->duplicate_candidates++;
                    return;
                }
            }

            stats->checked_candidates++;
            uint64_t candidate_prefilter_patterns = 0;
            if (params.prefilter_count != 0) {
                CheckResult prefilter;
                Rng prefilter_rng(
                    splitmix64_value(attempt_seed ^ UINT64_C(0xa4d1f3c9b27e596d)));
                bool capture_prefilter_failure =
                    !first_failure_recorded.load(std::memory_order_relaxed);
                stress_prefilter->reset(prefilter_rng, &prefilter, capture_prefilter_failure);
                stress_prefilter->run();
                candidate_prefilter_patterns = prefilter.patterns_checked;
                stats->prefilter_candidates++;
                stats->prefilter_patterns_checked += candidate_prefilter_patterns;

                if (!prefilter.is_mr) {
                    stats->total_patterns_checked += candidate_prefilter_patterns;
                    stats->prefilter_rejected_candidates++;
                    stats->failed_candidates++;

                    if (!first_failure_recorded.load(std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lock(result_mutex);
                        if (!first_failure_recorded.load(std::memory_order_relaxed)) {
                            code->patterns_checked = prefilter.patterns_checked;
                            result.code = *code;
                            result.check = std::move(prefilter);
                            first_failure_recorded.store(true, std::memory_order_relaxed);
                        }
                    }
                    return;
                }
            }

            CheckResult check = check_mr(*code);
            check.patterns_checked += candidate_prefilter_patterns;
            stats->total_patterns_checked += check.patterns_checked;
            if (check.strict_complete) {
                stats->strict_candidates++;
            }

            if (check.is_mr) {
                bool expected = false;
                if (stop.compare_exchange_strong(expected, true)) {
                    code->patterns_checked = check.patterns_checked;
                    std::lock_guard<std::mutex> lock(result_mutex);
                    result.status = 0;
                    result.message = "MR-LRC matrix found on attempt " + std::to_string(attempt);
                    result.code = *code;
                    result.check = std::move(check);
                }
                return;
            }

            stats->failed_candidates++;
            if (!first_failure_recorded.load(std::memory_order_relaxed)) {
                std::lock_guard<std::mutex> lock(result_mutex);
                if (!first_failure_recorded.load(std::memory_order_relaxed)) {
                    code->patterns_checked = check.patterns_checked;
                    result.code = *code;
                    result.check = std::move(check);
                    first_failure_recorded.store(true, std::memory_order_relaxed);
                }
            }
        };

        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(worker_count));

        for (uint64_t worker = 0; worker < worker_count; worker++) {
            workers.emplace_back([&, worker]() {
                WorkerStats local_stats;
                Code worker_code = base_code;
                std::unique_ptr<RandomStressPrefilter> stress_prefilter;
                if (params.prefilter_count != 0) {
                    stress_prefilter =
                        std::make_unique<RandomStressPrefilter>(worker_code, prefilter_plan);
                }
                uint64_t attempts_since_publish = 0;
                auto publish_local_stats = [&]() {
                    publish_stats(static_cast<std::size_t>(worker), local_stats);
                };
                auto finish_local_stats = [&]() {
                    publish_local_stats();
                    final_worker_stats[static_cast<std::size_t>(worker)] = local_stats;
                };
                try {
                    for (uint64_t attempt = worker + 1; attempt <= params.random_limit;) {
                        if (stop.load(std::memory_order_relaxed)) {
                            break;
                        }
                        run_attempt(attempt, &local_stats, &worker_code, stress_prefilter.get());
                        attempts_since_publish++;
                        if (attempts_since_publish >= kStatsPublishInterval) {
                            publish_local_stats();
                            attempts_since_publish = 0;
                        }
                        if (attempt > params.random_limit - worker_count) {
                            break;
                        }
                        attempt += worker_count;
                    }
                    finish_local_stats();
                } catch (const std::exception &ex) {
                    finish_local_stats();
                    {
                        std::lock_guard<std::mutex> lock(result_mutex);
                        if (worker_error.empty()) {
                            worker_error = ex.what();
                        }
                    }
                    stop.store(true, std::memory_order_relaxed);
                } catch (...) {
                    finish_local_stats();
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
            (void)finish_stats();
            result.status = -1;
            result.message = progress_error;
            return result;
        }

        if (!worker_error.empty()) {
            (void)finish_stats();
            result.status = -1;
            result.message = worker_error;
            return result;
        }

        WorkerStats final_stats = finish_stats();
        if (result.status == 0) {
            return result;
        }

        result.status = 2;
        result.message = "no MR-LRC matrix found after " +
                         std::to_string(final_stats.attempts_done) + " attempts";
        result.check.patterns_checked = final_stats.total_patterns_checked;
        result.check.failures = final_stats.failed_candidates;
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
