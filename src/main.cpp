#include "gf256.h"
#include "mr_lrc.h"

#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

class CliError : public std::runtime_error {
public:
    explicit CliError(const std::string &message) : std::runtime_error(message) {}
};

void print_usage(std::ostream &out, const char *program)
{
    out << "Usage: " << program
        << " --data K --groups G --local-parity R --global-parity M --seed S "
           "[--construction BOOL] [--local-method METHOD --global-method METHOD] "
           "[--random-limit N] [--thread-count N]\n\n"
        << "Generate and exhaustively verify a GF(256) MR-LRC generator matrix.\n\n"
        << "Required parameters:\n"
        << "  -k, --data K             Number of systematic data symbols.\n"
        << "  -g, --groups G           Number of local groups; data is split nearly evenly.\n"
        << "  -r, --local-parity R     Local parity symbols per group.\n"
        << "  -p, --global-parity M    Global parity symbols.\n"
        << "  -s, --seed S             Random seed.\n"
        << "      --local-method M     Local parity construction: cauchy, vandermonde, random.\n"
        << "      --global-method M    Global parity construction: cauchy, vandermonde, random.\n"
        << "      --construction BOOL  Enable registered data-local constructions: true/false, default true.\n"
        << "      --random-limit N     Maximum random candidate attempts, default 1000.\n"
        << "  -t, --thread-count N     Parallel search worker count, default 1, max 256.\n"
        << "  -m, --method M           Alias: set both local and global methods.\n"
        << "  -h, --help               Show this help.\n";
}

int parse_i32(const std::string &text)
{
    char *end = nullptr;
    errno = 0;
    long parsed = std::strtol(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0' || parsed < INT_MIN || parsed > INT_MAX) {
        throw CliError("invalid integer: " + text);
    }
    return static_cast<int>(parsed);
}

uint64_t parse_u64(const std::string &text)
{
    char *end = nullptr;
    errno = 0;
    unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0') {
        throw CliError("invalid unsigned integer: " + text);
    }
    return static_cast<uint64_t>(parsed);
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

const char *require_value(int argc, char **argv, int *i)
{
    if (*i + 1 >= argc) {
        throw CliError(std::string("missing value for ") + argv[*i]);
    }
    (*i)++;
    return argv[*i];
}

struct CliOptions {
    mrlrc::Params params;
    bool have_data = false;
    bool have_groups = false;
    bool have_local_parity = false;
    bool have_global_parity = false;
    bool have_seed = false;
    bool have_local_method = false;
    bool have_global_method = false;
};

CliOptions parse_args(int argc, char **argv)
{
    CliOptions opts;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_usage(std::cout, argv[0]);
            std::exit(0);
        } else if (arg == "-k" || arg == "--data") {
            opts.params.data = parse_i32(require_value(argc, argv, &i));
            opts.have_data = true;
        } else if (arg == "-g" || arg == "--groups") {
            opts.params.groups = parse_i32(require_value(argc, argv, &i));
            opts.have_groups = true;
        } else if (arg == "-r" || arg == "--local-parity") {
            opts.params.local_parity = parse_i32(require_value(argc, argv, &i));
            opts.have_local_parity = true;
        } else if (arg == "-p" || arg == "--global-parity") {
            opts.params.global_parity = parse_i32(require_value(argc, argv, &i));
            opts.have_global_parity = true;
        } else if (arg == "-s" || arg == "--seed") {
            opts.params.seed = parse_u64(require_value(argc, argv, &i));
            opts.have_seed = true;
        } else if (arg == "--random-limit") {
            opts.params.random_limit = parse_u64(require_value(argc, argv, &i));
        } else if (arg == "-t" || arg == "--thread-count") {
            opts.params.thread_count = parse_u64(require_value(argc, argv, &i));
        } else if (arg == "--local-method") {
            std::string method = require_value(argc, argv, &i);
            if (!mrlrc::parse_family(method, &opts.params.local_family)) {
                throw CliError("invalid local method: " + method);
            }
            opts.have_local_method = true;
        } else if (arg == "--global-method") {
            std::string method = require_value(argc, argv, &i);
            if (!mrlrc::parse_family(method, &opts.params.global_family)) {
                throw CliError("invalid global method: " + method);
            }
            opts.have_global_method = true;
        } else if (arg == "--construction") {
            opts.params.construction = parse_bool(require_value(argc, argv, &i));
        } else if (arg == "-m" || arg == "--method" || arg == "--family") {
            std::string method = require_value(argc, argv, &i);
            mrlrc::MatrixFamily parsed;
            if (!mrlrc::parse_family(method, &parsed)) {
                throw CliError("invalid method: " + method);
            }
            opts.params.local_family = parsed;
            opts.params.global_family = parsed;
            opts.have_local_method = true;
            opts.have_global_method = true;
        } else {
            throw CliError("unknown argument: " + arg);
        }
    }

    if (!opts.have_data || !opts.have_groups || !opts.have_local_parity ||
        !opts.have_global_parity || !opts.have_seed) {
        throw CliError("missing required parameter; run with --help for usage");
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

} // namespace

int main(int argc, char **argv)
{
    try {
        CliOptions opts = parse_args(argc, argv);
        mrlrc::GenerateResult result = mrlrc::generate(opts.params);
        const auto &code = result.code;

        if (result.status != 0) {
            std::cout << "MR-LRC matrix not found\n";
            std::cout << "message=" << result.message << "\n";
            if (!code.matrix.empty()) {
                std::cout << "patterns_checked=" << result.check.patterns_checked
                          << " failures=" << result.check.failures << "\n";
                print_failed_pattern(code, result.check.first_failed_erased);
            }
            return result.status == 2 ? 2 : 1;
        }

        std::cout << "MR-LRC matrix found\n";
        std::cout << "data=" << code.data << " groups=" << code.groups.size()
                  << " local_parity=" << opts.params.local_parity
                  << " local_rows=" << code.local_rows
                  << " global_parity=" << code.global_parity
                  << " total_parity=" << code.total_parity
                  << " symbols=" << code.symbols << "\n";
        std::cout << "construction=" << (code.construction ? "on" : "off")
                  << " candidate_source=" << code.candidate_source
                  << " local_method=" << mrlrc::family_name(code.local_family)
                  << " global_method=" << mrlrc::family_name(code.global_family)
                  << " seed=" << code.seed
                  << " random_limit=" << opts.params.random_limit
                  << " thread_count=" << opts.params.thread_count
                  << " attempt=" << code.attempt
                  << " patterns_checked=" << result.check.patterns_checked
                  << " gf256_backend=" << mrlrc::gf256_backend() << "\n";
        print_layout(code);
        print_matrix(code);
        return 0;
    } catch (const std::exception &ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 1;
    }
}
