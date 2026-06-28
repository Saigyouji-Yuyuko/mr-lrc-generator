#ifndef MR_LRC_GENERATOR_H
#define MR_LRC_GENERATOR_H

#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace mrlrc {

enum class MatrixFamily {
    Cauchy,
    ColumnMultiplierCauchy,
    Vandermonde,
    Random,
};

struct Params {
    int data = 0;
    int groups = 0;
    int local_parity = 0;
    int global_parity = 0;
    uint64_t seed = 0;
    uint64_t random_limit = std::numeric_limits<uint64_t>::max();
    uint64_t thread_count = 1;
    MatrixFamily local_family = MatrixFamily::Cauchy;
    MatrixFamily global_family = MatrixFamily::Cauchy;
    bool construction = true;
    bool cauchy_dedup = false;
    uint64_t step_time = 0;
    std::function<void(uint64_t)> progress_callback;
};

struct LocalGroup {
    std::vector<int> data;
    int local_parity = 0;
    int local_row_start = 0;
};

struct Code {
    int data = 0;
    int total_parity = 0;
    int global_parity = 0;
    int local_rows = 0;
    int symbols = 0;
    uint64_t seed = 0;
    uint64_t attempt = 0;
    uint64_t patterns_checked = 0;
    MatrixFamily local_family = MatrixFamily::Cauchy;
    MatrixFamily global_family = MatrixFamily::Cauchy;
    bool construction = true;
    std::string candidate_source = "search";
    std::vector<LocalGroup> groups;
    std::vector<uint8_t> matrix;
};

struct CheckResult {
    bool is_mr = false;
    uint64_t patterns_checked = 0;
    uint64_t failures = 0;
    std::string message;
    std::vector<int> first_failed_erased;
};

struct GenerateResult {
    int status = 0;
    std::string message;
    Code code;
    CheckResult check;
    uint64_t attempts_done = 0;
    uint64_t unique_candidates_checked = 0;
    uint64_t duplicate_candidates_skipped = 0;
    uint64_t cauchy_dedup_key_bytes = 0;
    bool cauchy_dedup_enabled = false;
};

const char *family_name(MatrixFamily family);
bool parse_family(const std::string &name, MatrixFamily *family);

GenerateResult generate(const Params &params);
CheckResult check_mr(const Code &code);
std::string symbol_label(const Code &code, int symbol);
std::vector<int> erased_symbols_from_mask(const Code &code, const std::vector<uint8_t> &erased);

}

#endif
