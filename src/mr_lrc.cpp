#include "mr_lrc.h"

#include "gf256.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
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

class Enumerator {
public:
    explicit Enumerator(const Code &code, CheckResult *result)
        : code_(code), result_(result), erased_(static_cast<std::size_t>(code.symbols), 0),
          choose_counts_(code.groups.size(), 0)
    {
    }

    void run()
    {
        enumerate_extra_counts(0, 0);
        result_->is_mr = result_->failures == 0;
        result_->message = result_->is_mr
                               ? "candidate is recoverable for every maximal erasure pattern"
                               : "candidate is singular for at least one maximal erasure pattern";
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
    Enumerator enumerator(code, &result);
    enumerator.run();
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
        std::mutex result_mutex;
        std::string worker_error;

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
