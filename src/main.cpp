#include "gf256.h"
#include "mr_lrc.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gflags/gflags.h>

DEFINE_int32(data, 0, "Number of systematic data symbols.");
DEFINE_int32(groups, 0, "Number of local groups; data is split nearly evenly.");
DEFINE_int32(local_parity, 0, "Local parity symbols per group.");
DEFINE_int32(global_parity, 0, "Global parity symbols.");
DEFINE_uint64(seed, 0, "Random seed. If omitted, a seed is generated and printed.");
DEFINE_string(local_method, "cauchy", "Local parity construction: cauchy, vandermonde, random.");
DEFINE_string(global_method, "cauchy",
              "Global parity construction: cauchy, column_multiplier_cauchy, vandermonde, random.");
DEFINE_string(method, "",
              "Alias that sets both local_method and global_method; global-only methods are rejected.");
DEFINE_string(construction, "true",
              "Enable registered data-local constructions: true/false, on/off, 1/0, yes/no.");
DEFINE_uint64(random_limit, std::numeric_limits<uint64_t>::max(), "Maximum random candidate attempts.");
DEFINE_uint64(thread_count, 1, "Parallel search worker count, max 256.");
DEFINE_uint64(step_time, 30, "Print search progress every N seconds; 0 disables progress.");
DEFINE_string(json, "", "Write the found matrix as pretty JSON to this file.");
DEFINE_bool(cauchy_dedup, false,
            "Skip duplicate all-Cauchy candidates using canonical Cauchy parameter keys.");

namespace {

class CliError : public std::runtime_error {
public:
    explicit CliError(const std::string &message) : std::runtime_error(message) {}
};

std::string usage_message()
{
    return "Generate and exactly verify a GF(256) MR-LRC generator matrix.\n\n"
           "Required parameters:\n"
           "  -k, --data K\n"
           "  -g, --groups G\n"
           "  -r, --local-parity R\n"
           "  -p, --global-parity M\n\n"
           "Optional aliases are accepted for compatibility:\n"
           "  -s, --seed S\n"
           "  -t, --thread-count N\n"
           "  -m, --method M\n"
           "  --local-method M, --global-method M, --random-limit N\n"
           "  --json FILE\n";
}

bool parse_bool(const std::string &text)
{
    if (text == "1" || text == "true" || text == "on" || text == "yes") {
        return true;
    }
    if (text == "0" || text == "false" || text == "off" || text == "no") {
        return false;
    }
    throw CliError("invalid boolean: " + text);
}

bool flag_was_set(const char *name)
{
    gflags::CommandLineFlagInfo info;
    return gflags::GetCommandLineFlagInfo(name, &info) && !info.is_default;
}

uint64_t splitmix64_value(uint64_t value)
{
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

uint64_t generate_seed()
{
    std::random_device random;
    uint64_t value = (static_cast<uint64_t>(random()) << 32U) ^ static_cast<uint64_t>(random());
    value ^= static_cast<uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    value ^= static_cast<uint64_t>(reinterpret_cast<std::uintptr_t>(&value));
    return splitmix64_value(value);
}

std::string normalize_flag_name(const std::string &arg)
{
    static const std::vector<std::pair<std::string, std::string>> aliases = {
        {"--local-parity", "--local_parity"},
        {"--global-parity", "--global_parity"},
        {"--local-method", "--local_method"},
        {"--global-method", "--global_method"},
        {"--random-limit", "--random_limit"},
        {"--thread-count", "--thread_count"},
        {"--step-time", "--step_time"},
        {"--step-times", "--step_time"},
        {"--matrix-json", "--json"},
        {"--cauchy-dedup", "--cauchy_dedup"},
        {"-k", "--data"},
        {"-g", "--groups"},
        {"-r", "--local_parity"},
        {"-p", "--global_parity"},
        {"-s", "--seed"},
        {"-t", "--thread_count"},
        {"-m", "--method"},
        {"-h", "--help"},
        {"--family", "--method"},
    };

    for (const auto &alias : aliases) {
        const auto &from = alias.first;
        const auto &to = alias.second;
        if (arg == from) {
            return to;
        }
        if (arg.rfind(from + "=", 0) == 0) {
            return to + arg.substr(from.size());
        }
    }

    if (arg == "--noconstruction") {
        return "--construction=false";
    }
    if (arg == "--no-cauchy-dedup") {
        return "--nocauchy_dedup";
    }
    return arg;
}

std::vector<std::string> normalize_args(int argc, char **argv)
{
    std::vector<std::string> normalized;
    normalized.reserve(static_cast<std::size_t>(argc));
    for (int i = 0; i < argc; i++) {
        normalized.push_back(normalize_flag_name(argv[i]));
    }
    return normalized;
}

struct CliOptions {
    mrlrc::Params params;
    std::string json_path;
};

CliOptions parse_args(int argc, char **argv)
{
    std::vector<std::string> normalized = normalize_args(argc, argv);
    std::vector<char *> normalized_argv;
    normalized_argv.reserve(normalized.size());
    for (auto &arg : normalized) {
        normalized_argv.push_back(arg.data());
    }

    int normalized_argc = static_cast<int>(normalized_argv.size());
    char **normalized_argv_data = normalized_argv.data();
    gflags::SetUsageMessage(usage_message());
    gflags::ParseCommandLineFlags(&normalized_argc, &normalized_argv_data, true);

    if (normalized_argc > 1) {
        throw CliError(std::string("unexpected positional argument: ") + normalized_argv_data[1]);
    }

    std::vector<std::string> missing;
    if (!flag_was_set("data")) {
        missing.push_back("--data");
    }
    if (!flag_was_set("groups")) {
        missing.push_back("--groups");
    }
    if (!flag_was_set("local_parity")) {
        missing.push_back("--local-parity");
    }
    if (!flag_was_set("global_parity")) {
        missing.push_back("--global-parity");
    }
    if (!missing.empty()) {
        std::string message = "missing required parameter(s): ";
        for (std::size_t i = 0; i < missing.size(); i++) {
            if (i != 0) {
                message += ", ";
            }
            message += missing[i];
        }
        message += "; run with --help for usage";
        throw CliError(message);
    }

    CliOptions opts;
    opts.params.data = FLAGS_data;
    opts.params.groups = FLAGS_groups;
    opts.params.local_parity = FLAGS_local_parity;
    opts.params.global_parity = FLAGS_global_parity;
    opts.params.seed = flag_was_set("seed") ? FLAGS_seed : generate_seed();
    opts.params.random_limit = FLAGS_random_limit;
    opts.params.thread_count = FLAGS_thread_count;
    opts.params.construction = parse_bool(FLAGS_construction);
    opts.params.cauchy_dedup = FLAGS_cauchy_dedup;
    opts.params.step_time = FLAGS_step_time;
    if (opts.params.step_time != 0) {
        opts.params.progress_callback = [](uint64_t step, uint64_t searched, uint64_t strict_checked) {
            std::time_t now = std::time(nullptr);
            std::tm local_time{};
            localtime_r(&now, &local_time);
            std::cerr << "step: " << step << " "
                      << std::put_time(&local_time, "%F %T")
                      << " search: " << searched
                      << ", strict search:" << strict_checked << ";\n";
        };
    }
    opts.json_path = FLAGS_json;

    if (!FLAGS_method.empty()) {
        mrlrc::MatrixFamily parsed;
        if (!mrlrc::parse_family(FLAGS_method, &parsed)) {
            throw CliError("invalid method: " + FLAGS_method);
        }
        opts.params.local_family = parsed;
        opts.params.global_family = parsed;
    }
    if (FLAGS_method.empty() || flag_was_set("local_method")) {
        if (!mrlrc::parse_family(FLAGS_local_method, &opts.params.local_family)) {
            throw CliError("invalid local method: " + FLAGS_local_method);
        }
    }
    if (FLAGS_method.empty() || flag_was_set("global_method")) {
        if (!mrlrc::parse_family(FLAGS_global_method, &opts.params.global_family)) {
            throw CliError("invalid global method: " + FLAGS_global_method);
        }
    }
    return opts;
}

void print_layout(const mrlrc::Code &code)
{
    std::cout << "layout:\n";
    for (std::size_t group_index = 0; group_index < code.groups.size(); group_index++) {
        const auto &group = code.groups[group_index];
        std::cout << "  group " << group_index << ": data";
        for (int data : group.data) {
            std::cout << " D" << data;
        }
        std::cout << ", local";
        for (int local = 0; local < group.local_parity; local++) {
            std::cout << " " << mrlrc::symbol_label(code, group.local_row_start + local);
        }
        std::cout << "\n";
    }
}

void print_matrix(const mrlrc::Code &code)
{
    std::cout << "matrix_hex (" << code.symbols << " x " << code.data << "):\n";
    for (int row = 0; row < code.symbols; row++) {
        std::cout << "  " << std::left << std::setw(8) << mrlrc::symbol_label(code, row) << std::right;
        for (int col = 0; col < code.data; col++) {
            auto pos = static_cast<std::size_t>(row) * static_cast<std::size_t>(code.data) +
                       static_cast<std::size_t>(col);
            std::cout << " " << std::hex << std::setw(2) << std::setfill('0')
                      << static_cast<unsigned>(code.matrix[pos]) << std::dec << std::setfill(' ');
        }
        std::cout << "\n";
    }
}

void print_matrix_json(std::ostream &out, const mrlrc::Code &code, int local_parity)
{
    out << "{\n";
    out << "  \"data_cnt\": " << code.data << ",\n";
    out << "  \"group_cnt\": " << code.groups.size() << ",\n";
    out << "  \"local_parity_cnt\": " << local_parity << ",\n";
    out << "  \"global_parity_cnt\": " << code.global_parity << ",\n";
    out << "  \"matrix\": [\n";
    for (int row = 0; row < code.symbols; row++) {
        out << "    [";
        for (int col = 0; col < code.data; col++) {
            if (col != 0) {
                out << ", ";
            }
            auto pos = static_cast<std::size_t>(row) * static_cast<std::size_t>(code.data) +
                       static_cast<std::size_t>(col);
            out << static_cast<unsigned>(code.matrix[pos]);
        }
        out << "]";
        if (row + 1 != code.symbols) {
            out << ",";
        }
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

void write_matrix_json_file(const std::string &path, const mrlrc::Code &code, int local_parity)
{
    std::ofstream out(path);
    if (!out) {
        throw CliError("could not open JSON output file: " + path);
    }
    print_matrix_json(out, code, local_parity);
    if (!out) {
        throw CliError("could not write JSON output file: " + path);
    }
}

void print_failed_pattern(const mrlrc::Code &code, const std::vector<int> &erased)
{
    if (erased.empty()) {
        return;
    }

    std::cout << "first_unrecoverable_erasure:";
    for (int symbol : erased) {
        std::cout << " " << mrlrc::symbol_label(code, symbol);
    }
    std::cout << "\n";
}

template <typename T>
void print_attribute(const char *name, const T &value)
{
    std::cout << name << "=" << value << "\n";
}

} // namespace

int main(int argc, char **argv)
{
    try {
        CliOptions opts = parse_args(argc, argv);
        mrlrc::GenerateResult result = mrlrc::generate(opts.params);
        const auto &code = result.code;

        if (result.status != 0) {
            std::cout << "MR-LRC matrix not found\n";
            print_attribute("message", result.message);
            if (!code.matrix.empty()) {
                print_attribute("patterns_checked", result.check.patterns_checked);
                print_attribute("failures", result.check.failures);
                if (result.cauchy_dedup_enabled) {
                    print_attribute("cauchy_dedup", "on");
                    print_attribute("cauchy_dedup_key_bytes", result.cauchy_dedup_key_bytes);
                    print_attribute("attempts_done", result.attempts_done);
                    print_attribute("unique_candidates_checked", result.unique_candidates_checked);
                    print_attribute("strict_candidates_checked", result.strict_candidates_checked);
                    print_attribute("duplicates_skipped", result.duplicate_candidates_skipped);
                }
                print_failed_pattern(code, result.check.first_failed_erased);
            }
            gflags::ShutDownCommandLineFlags();
            return result.status == 2 ? 2 : 1;
        }

        if (!opts.json_path.empty()) {
            write_matrix_json_file(opts.json_path, code, opts.params.local_parity);
        }

        std::cout << "MR-LRC matrix found\n";
        print_attribute("data", code.data);
        print_attribute("groups", code.groups.size());
        print_attribute("local_parity", opts.params.local_parity);
        print_attribute("local_rows", code.local_rows);
        print_attribute("global_parity", code.global_parity);
        print_attribute("total_parity", code.total_parity);
        print_attribute("symbols", code.symbols);
        print_attribute("construction", code.construction ? "on" : "off");
        print_attribute("candidate_source", code.candidate_source);
        print_attribute("local_method", mrlrc::family_name(code.local_family));
        print_attribute("global_method", mrlrc::family_name(code.global_family));
        print_attribute("seed", code.seed);
        print_attribute("random_limit", opts.params.random_limit);
        print_attribute("thread_count", opts.params.thread_count);
        print_attribute("attempt", code.attempt);
        print_attribute("patterns_checked", result.check.patterns_checked);
        if (result.cauchy_dedup_enabled) {
            print_attribute("cauchy_dedup", "on");
            print_attribute("cauchy_dedup_key_bytes", result.cauchy_dedup_key_bytes);
            print_attribute("attempts_done", result.attempts_done);
            print_attribute("unique_candidates_checked", result.unique_candidates_checked);
            print_attribute("strict_candidates_checked", result.strict_candidates_checked);
            print_attribute("duplicates_skipped", result.duplicate_candidates_skipped);
        }
        print_attribute("gf256_backend", mrlrc::gf256_backend());
        print_layout(code);
        print_matrix(code);
        gflags::ShutDownCommandLineFlags();
        return 0;
    } catch (const std::exception &ex) {
        std::cerr << "error: " << ex.what() << "\n";
        gflags::ShutDownCommandLineFlags();
        return 1;
    }
}
