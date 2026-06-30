#include "gf256.h"
#include "mr_lrc.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
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
DEFINE_uint64(construction, 0,
              "Difference-pack construction attempts before random search; 0 disables.");
DEFINE_uint64(random_limit, std::numeric_limits<uint64_t>::max(), "Maximum random candidate attempts.");
DEFINE_uint64(thread_count, 1, "Parallel search worker count, max 256.");
DEFINE_uint64(step_time, 30, "Print search progress every N seconds; 0 disables progress.");
DEFINE_string(json, "", "Write the found matrix as pretty JSON to this file.");
DEFINE_string(check_json, "", "Read a matrix JSON file and exactly verify whether it is MR-LRC.");
DEFINE_bool(cauchy_dedup, false,
            "Skip duplicate all-Cauchy candidates using canonical Cauchy parameter keys.");
DEFINE_uint64(prefilter_count, 0,
              "Stress-test patterns per candidate before exact verification; 0 disables.");

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
           "  --prefilter-count N\n"
           "  --construction N\n"
           "  --json FILE\n"
           "  --check-json FILE\n";
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
        {"--check-matrix-json", "--check_json"},
        {"--input-json", "--check_json"},
        {"--verify-json", "--check_json"},
        {"--cauchy-dedup", "--cauchy_dedup"},
        {"--prefilter-count", "--prefilter_count"},
        {"--random-prefilter", "--prefilter_count"},
        {"--stress-prefilter", "--prefilter_count"},
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

    auto normalize_construction_value = [](const std::string &value) {
        if (value == "true" || value == "on" || value == "yes") {
            return std::string("1");
        }
        if (value == "false" || value == "off" || value == "no") {
            return std::string("0");
        }
        return value;
    };

    if (arg.rfind("--construction=", 0) == 0) {
        return "--construction=" + normalize_construction_value(arg.substr(15));
    }
    if (arg == "--noconstruction") {
        return "--construction=0";
    }
    if (arg == "--no-cauchy-dedup") {
        return "--nocauchy_dedup";
    }
    return arg;
}

std::vector<std::string> normalize_args(int argc, char **argv)
{
    auto normalize_construction_value = [](const std::string &value) {
        if (value == "true" || value == "on" || value == "yes") {
            return std::string("1");
        }
        if (value == "false" || value == "off" || value == "no") {
            return std::string("0");
        }
        return value;
    };

    std::vector<std::string> normalized;
    normalized.reserve(static_cast<std::size_t>(argc));
    for (int i = 0; i < argc; i++) {
        std::string arg = normalize_flag_name(argv[i]);
        if (!normalized.empty() && normalized.back() == "--construction") {
            arg = normalize_construction_value(arg);
        }
        normalized.push_back(std::move(arg));
    }
    return normalized;
}

struct CliOptions {
    mrlrc::Params params;
    std::string json_path;
    std::string check_json_path;
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

    CliOptions opts;
    opts.json_path = FLAGS_json;
    opts.check_json_path = FLAGS_check_json;
    if (!opts.check_json_path.empty()) {
        return opts;
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

    opts.params.data = FLAGS_data;
    opts.params.groups = FLAGS_groups;
    opts.params.local_parity = FLAGS_local_parity;
    opts.params.global_parity = FLAGS_global_parity;
    opts.params.seed = flag_was_set("seed") ? FLAGS_seed : generate_seed();
    opts.params.random_limit = FLAGS_random_limit;
    opts.params.thread_count = FLAGS_thread_count;
    opts.params.construction = FLAGS_construction;
    opts.params.cauchy_dedup = FLAGS_cauchy_dedup;
    opts.params.prefilter_count = FLAGS_prefilter_count;
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

void print_json_string(std::ostream &out, const std::string &value)
{
    out << '"';
    for (unsigned char ch : value) {
        switch (ch) {
        case '"':
            out << "\\\"";
            break;
        case '\\':
            out << "\\\\";
            break;
        case '\b':
            out << "\\b";
            break;
        case '\f':
            out << "\\f";
            break;
        case '\n':
            out << "\\n";
            break;
        case '\r':
            out << "\\r";
            break;
        case '\t':
            out << "\\t";
            break;
        default:
            if (ch < 0x20) {
                auto flags = out.flags();
                auto fill = out.fill();
                out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                    << static_cast<unsigned>(ch);
                out.flags(flags);
                out.fill(fill);
            } else {
                out << static_cast<char>(ch);
            }
            break;
        }
    }
    out << '"';
}

void print_json_array(std::ostream &out, const std::vector<int> &values)
{
    out << "[";
    for (std::size_t i = 0; i < values.size(); i++) {
        if (i != 0) {
            out << ", ";
        }
        out << values[i];
    }
    out << "]";
}

void print_json_hex_byte(std::ostream &out, uint8_t value)
{
    auto flags = out.flags();
    auto fill = out.fill();
    out << '"' << std::hex << std::setw(2) << std::setfill('0')
        << static_cast<unsigned>(value) << '"';
    out.flags(flags);
    out.fill(fill);
}

void print_groups_json(std::ostream &out, const mrlrc::Code &code)
{
    out << "  \"groups\": [\n";
    for (std::size_t group_index = 0; group_index < code.groups.size(); group_index++) {
        const auto &group = code.groups[group_index];
        std::vector<int> local;
        local.reserve(static_cast<std::size_t>(group.local_parity));
        for (int local_index = 0; local_index < group.local_parity; local_index++) {
            local.push_back(group.local_row_start + local_index);
        }

        out << "    { \"data\": ";
        print_json_array(out, group.data);
        out << ", \"local\": ";
        print_json_array(out, local);
        out << " }";
        if (group_index + 1 != code.groups.size()) {
            out << ",";
        }
        out << "\n";
    }
    out << "  ],\n";
}

void print_row_labels_json(std::ostream &out, const mrlrc::Code &code)
{
    out << "  \"row_labels\": [";
    for (int row = 0; row < code.symbols; row++) {
        if (row != 0) {
            out << ", ";
        }
        print_json_string(out, mrlrc::symbol_label(code, row));
    }
    out << "],\n";
}

void print_matrix_decimal_json(std::ostream &out, const mrlrc::Code &code)
{
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
    out << "  ],\n";
}

void print_matrix_hex_json(std::ostream &out, const mrlrc::Code &code)
{
    out << "  \"matrix_hex\": [\n";
    for (int row = 0; row < code.symbols; row++) {
        out << "    [";
        for (int col = 0; col < code.data; col++) {
            if (col != 0) {
                out << ", ";
            }
            auto pos = static_cast<std::size_t>(row) * static_cast<std::size_t>(code.data) +
                       static_cast<std::size_t>(col);
            print_json_hex_byte(out, code.matrix[pos]);
        }
        out << "]";
        if (row + 1 != code.symbols) {
            out << ",";
        }
        out << "\n";
    }
    out << "  ]\n";
}

void print_matrix_json(std::ostream &out, const mrlrc::GenerateResult &result,
                       const mrlrc::Params &params)
{
    const auto &code = result.code;
    out << "{\n";
    out << "  \"status\": \"found\",\n";
    out << "  \"data\": " << code.data << ",\n";
    out << "  \"data_cnt\": " << code.data << ",\n";
    out << "  \"group_cnt\": " << code.groups.size() << ",\n";
    out << "  \"local_parity\": " << params.local_parity << ",\n";
    out << "  \"local_parity_cnt\": " << params.local_parity << ",\n";
    out << "  \"local_rows\": " << code.local_rows << ",\n";
    out << "  \"global_parity\": " << code.global_parity << ",\n";
    out << "  \"global_parity_cnt\": " << code.global_parity << ",\n";
    out << "  \"total_parity\": " << code.total_parity << ",\n";
    out << "  \"symbols\": " << code.symbols << ",\n";
    out << "  \"construction\": " << code.construction << ",\n";
    out << "  \"candidate_source\": ";
    print_json_string(out, code.candidate_source);
    out << ",\n";
    out << "  \"local_method\": ";
    print_json_string(out, mrlrc::family_name(code.local_family));
    out << ",\n";
    out << "  \"global_method\": ";
    print_json_string(out, mrlrc::family_name(code.global_family));
    out << ",\n";
    out << "  \"seed\": " << code.seed << ",\n";
    out << "  \"random_limit\": " << params.random_limit << ",\n";
    out << "  \"thread_count\": " << params.thread_count << ",\n";
    out << "  \"attempt\": " << code.attempt << ",\n";
    out << "  \"patterns_checked\": " << result.check.patterns_checked << ",\n";
    out << "  \"strict_complete\": "
        << (result.check.strict_complete ? "true" : "false") << ",\n";
    if (result.prefilter_enabled) {
        out << "  \"prefilter_count\": " << params.prefilter_count << ",\n";
        out << "  \"prefilter_candidates_checked\": "
            << result.prefilter_candidates_checked << ",\n";
        out << "  \"prefilter_candidates_rejected\": "
            << result.prefilter_candidates_rejected << ",\n";
        out << "  \"prefilter_patterns_checked\": "
            << result.prefilter_patterns_checked << ",\n";
    }
    if (result.cauchy_dedup_enabled) {
        out << "  \"cauchy_dedup\": true,\n";
        out << "  \"cauchy_dedup_key_bytes\": " << result.cauchy_dedup_key_bytes << ",\n";
        out << "  \"attempts_done\": " << result.attempts_done << ",\n";
        out << "  \"unique_candidates_checked\": " << result.unique_candidates_checked << ",\n";
        out << "  \"strict_candidates_checked\": " << result.strict_candidates_checked << ",\n";
        out << "  \"duplicates_skipped\": " << result.duplicate_candidates_skipped << ",\n";
    }
    out << "  \"gf256_backend\": ";
    print_json_string(out, mrlrc::gf256_backend());
    out << ",\n";
    print_groups_json(out, code);
    print_row_labels_json(out, code);
    print_matrix_decimal_json(out, code);
    print_matrix_hex_json(out, code);
    out << "}\n";
}

void write_matrix_json_file(const std::string &path, const mrlrc::GenerateResult &result,
                            const mrlrc::Params &params)
{
    std::ofstream out(path);
    if (!out) {
        throw CliError("could not open JSON output file: " + path);
    }
    print_matrix_json(out, result, params);
    if (!out) {
        throw CliError("could not write JSON output file: " + path);
    }
}

enum class JsonType {
    Null,
    Bool,
    Number,
    String,
    Array,
    Object,
};

struct JsonValue {
    JsonType type = JsonType::Null;
    bool bool_value = false;
    std::string text;
    std::vector<JsonValue> array;
    std::map<std::string, JsonValue> object;
};

class JsonParser {
public:
    explicit JsonParser(std::string text) : text_(std::move(text)) {}

    JsonValue parse()
    {
        JsonValue value = parse_value();
        skip_ws();
        if (pos_ != text_.size()) {
            throw CliError("unexpected trailing content in JSON");
        }
        return value;
    }

private:
    char peek() const
    {
        return pos_ < text_.size() ? text_[pos_] : '\0';
    }

    char take()
    {
        if (pos_ >= text_.size()) {
            throw CliError("unexpected end of JSON");
        }
        return text_[pos_++];
    }

    void skip_ws()
    {
        while (pos_ < text_.size() &&
               std::isspace(static_cast<unsigned char>(text_[pos_])) != 0) {
            pos_++;
        }
    }

    void expect(char expected)
    {
        char actual = take();
        if (actual != expected) {
            throw CliError(std::string("expected '") + expected + "' in JSON");
        }
    }

    bool consume_literal(const char *literal)
    {
        std::size_t start = pos_;
        for (const char *p = literal; *p != '\0'; p++) {
            if (pos_ >= text_.size() || text_[pos_] != *p) {
                pos_ = start;
                return false;
            }
            pos_++;
        }
        return true;
    }

    JsonValue parse_value()
    {
        skip_ws();
        char ch = peek();
        if (ch == '{') {
            return parse_object();
        }
        if (ch == '[') {
            return parse_array();
        }
        if (ch == '"') {
            JsonValue value;
            value.type = JsonType::String;
            value.text = parse_string();
            return value;
        }
        if (ch == '-' || std::isdigit(static_cast<unsigned char>(ch)) != 0) {
            JsonValue value;
            value.type = JsonType::Number;
            value.text = parse_number();
            return value;
        }
        if (consume_literal("true")) {
            JsonValue value;
            value.type = JsonType::Bool;
            value.bool_value = true;
            return value;
        }
        if (consume_literal("false")) {
            JsonValue value;
            value.type = JsonType::Bool;
            value.bool_value = false;
            return value;
        }
        if (consume_literal("null")) {
            JsonValue value;
            value.type = JsonType::Null;
            return value;
        }
        throw CliError("invalid JSON value");
    }

    JsonValue parse_object()
    {
        JsonValue value;
        value.type = JsonType::Object;
        expect('{');
        skip_ws();
        if (peek() == '}') {
            take();
            return value;
        }

        for (;;) {
            skip_ws();
            if (peek() != '"') {
                throw CliError("expected object key in JSON");
            }
            std::string key = parse_string();
            skip_ws();
            expect(':');
            value.object.emplace(std::move(key), parse_value());
            skip_ws();
            char ch = take();
            if (ch == '}') {
                break;
            }
            if (ch != ',') {
                throw CliError("expected ',' or '}' in JSON object");
            }
        }
        return value;
    }

    JsonValue parse_array()
    {
        JsonValue value;
        value.type = JsonType::Array;
        expect('[');
        skip_ws();
        if (peek() == ']') {
            take();
            return value;
        }

        for (;;) {
            value.array.push_back(parse_value());
            skip_ws();
            char ch = take();
            if (ch == ']') {
                break;
            }
            if (ch != ',') {
                throw CliError("expected ',' or ']' in JSON array");
            }
        }
        return value;
    }

    std::string parse_string()
    {
        expect('"');
        std::string value;
        while (pos_ < text_.size()) {
            char ch = take();
            if (ch == '"') {
                return value;
            }
            if (static_cast<unsigned char>(ch) < 0x20) {
                throw CliError("control character in JSON string");
            }
            if (ch != '\\') {
                value.push_back(ch);
                continue;
            }

            char escaped = take();
            switch (escaped) {
            case '"':
            case '\\':
            case '/':
                value.push_back(escaped);
                break;
            case 'b':
                value.push_back('\b');
                break;
            case 'f':
                value.push_back('\f');
                break;
            case 'n':
                value.push_back('\n');
                break;
            case 'r':
                value.push_back('\r');
                break;
            case 't':
                value.push_back('\t');
                break;
            case 'u':
                parse_ascii_unicode_escape(&value);
                break;
            default:
                throw CliError("invalid JSON string escape");
            }
        }
        throw CliError("unterminated JSON string");
    }

    int hex_digit(char ch) const
    {
        if (ch >= '0' && ch <= '9') {
            return ch - '0';
        }
        if (ch >= 'a' && ch <= 'f') {
            return 10 + ch - 'a';
        }
        if (ch >= 'A' && ch <= 'F') {
            return 10 + ch - 'A';
        }
        return -1;
    }

    void parse_ascii_unicode_escape(std::string *out)
    {
        int value = 0;
        for (int i = 0; i < 4; i++) {
            int digit = hex_digit(take());
            if (digit < 0) {
                throw CliError("invalid JSON unicode escape");
            }
            value = (value << 4) | digit;
        }
        if (value > 0x7f) {
            throw CliError("only ASCII JSON unicode escapes are supported");
        }
        out->push_back(static_cast<char>(value));
    }

    std::string parse_number()
    {
        std::size_t start = pos_;
        if (peek() == '-') {
            pos_++;
        }
        if (std::isdigit(static_cast<unsigned char>(peek())) == 0) {
            throw CliError("invalid JSON number");
        }
        while (std::isdigit(static_cast<unsigned char>(peek())) != 0) {
            pos_++;
        }
        return text_.substr(start, pos_ - start);
    }

    std::string text_;
    std::size_t pos_ = 0;
};

std::string read_text_file(const std::string &path)
{
    std::ifstream in(path);
    if (!in) {
        throw CliError("could not open JSON input file: " + path);
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

const JsonValue *json_field(const JsonValue &object, const std::string &name)
{
    if (object.type != JsonType::Object) {
        return nullptr;
    }
    auto found = object.object.find(name);
    return found == object.object.end() ? nullptr : &found->second;
}

const JsonValue *json_any_field(const JsonValue &object, std::initializer_list<const char *> names)
{
    for (const char *name : names) {
        const JsonValue *value = json_field(object, name);
        if (value != nullptr) {
            return value;
        }
    }
    return nullptr;
}

int json_int_value(const JsonValue &value, const std::string &name)
{
    if (value.type != JsonType::Number) {
        throw CliError("JSON field '" + name + "' must be an integer");
    }
    std::size_t used = 0;
    long long parsed = 0;
    try {
        parsed = std::stoll(value.text, &used, 10);
    } catch (const std::exception &) {
        throw CliError("JSON field '" + name + "' is out of integer range");
    }
    if (used != value.text.size()) {
        throw CliError("JSON field '" + name + "' must be an integer");
    }
    if (parsed < std::numeric_limits<int>::min() || parsed > std::numeric_limits<int>::max()) {
        throw CliError("JSON field '" + name + "' is out of int range");
    }
    return static_cast<int>(parsed);
}

int required_json_int(const JsonValue &root, std::initializer_list<const char *> names)
{
    const JsonValue *value = json_any_field(root, names);
    if (value == nullptr) {
        return std::numeric_limits<int>::min();
    }
    return json_int_value(*value, *names.begin());
}

int optional_json_int(const JsonValue &root, std::initializer_list<const char *> names, int fallback)
{
    const JsonValue *value = json_any_field(root, names);
    return value == nullptr ? fallback : json_int_value(*value, *names.begin());
}

uint64_t optional_json_uint64(const JsonValue &root, const std::string &name, uint64_t fallback)
{
    const JsonValue *value = json_field(root, name);
    if (value == nullptr) {
        return fallback;
    }
    if (value->type == JsonType::Bool) {
        return value->bool_value ? 1 : 0;
    }
    if (value->type == JsonType::Number || value->type == JsonType::String) {
        std::string text = value->text;
        if (text == "true" || text == "on" || text == "yes") {
            return 1;
        }
        if (text == "false" || text == "off" || text == "no") {
            return 0;
        }
        if (text.empty() || text.front() == '-') {
            throw CliError("JSON field '" + name + "' must be a non-negative integer");
        }
        std::size_t used = 0;
        unsigned long long parsed = 0;
        try {
            parsed = std::stoull(text, &used, 10);
        } catch (const std::exception &) {
            throw CliError("JSON field '" + name + "' is out of uint64 range");
        }
        if (used != text.size()) {
            throw CliError("JSON field '" + name + "' must be a non-negative integer");
        }
        return static_cast<uint64_t>(parsed);
    }
    throw CliError("JSON field '" + name + "' must be a non-negative integer");
}

std::string optional_json_string(const JsonValue &root, const std::string &name,
                                 const std::string &fallback)
{
    const JsonValue *value = json_field(root, name);
    if (value == nullptr) {
        return fallback;
    }
    if (value->type != JsonType::String) {
        throw CliError("JSON field '" + name + "' must be a string");
    }
    return value->text;
}

int json_array_int(const JsonValue &value, const std::string &name)
{
    return json_int_value(value, name);
}

bool is_hex_text(const std::string &text)
{
    if (text.rfind("0x", 0) == 0 || text.rfind("0X", 0) == 0) {
        return true;
    }
    for (char ch : text) {
        if ((ch >= 'a' && ch <= 'f') || (ch >= 'A' && ch <= 'F')) {
            return true;
        }
    }
    return false;
}

uint8_t parse_json_byte(const JsonValue &value, const std::string &name, bool hex_context)
{
    int parsed = 0;
    if (value.type == JsonType::Number) {
        parsed = json_int_value(value, name);
    } else if (value.type == JsonType::String) {
        int base = hex_context || is_hex_text(value.text) ? 16 : 10;
        std::string text = value.text;
        if (base == 16 && (text.rfind("0x", 0) == 0 || text.rfind("0X", 0) == 0)) {
            text = text.substr(2);
        }
        if (text.empty()) {
            throw CliError("empty byte value in JSON matrix");
        }
        std::size_t used = 0;
        try {
            parsed = std::stoi(text, &used, base);
        } catch (const std::exception &) {
            throw CliError("invalid byte value in JSON matrix");
        }
        if (used != text.size()) {
            throw CliError("invalid byte value in JSON matrix");
        }
    } else {
        throw CliError("JSON matrix entries must be integers or hex strings");
    }
    if (parsed < 0 || parsed > 255) {
        throw CliError("JSON matrix byte value must be in 0..255");
    }
    return static_cast<uint8_t>(parsed);
}

std::vector<int> parse_int_array(const JsonValue &value, const std::string &name)
{
    if (value.type != JsonType::Array) {
        throw CliError("JSON field '" + name + "' must be an array");
    }
    std::vector<int> result;
    result.reserve(value.array.size());
    for (const JsonValue &entry : value.array) {
        result.push_back(json_array_int(entry, name));
    }
    return result;
}

void apply_json_groups(const JsonValue &groups_json, mrlrc::Code *code)
{
    if (groups_json.type != JsonType::Array) {
        throw CliError("JSON field 'groups' must be an array");
    }
    if (groups_json.array.size() != code->groups.size()) {
        throw CliError("JSON group count does not match group_cnt");
    }

    std::vector<uint8_t> seen_data(static_cast<std::size_t>(code->data), 0);
    int next_local_row = code->data;
    for (std::size_t group_index = 0; group_index < groups_json.array.size(); group_index++) {
        const JsonValue &group_json = groups_json.array[group_index];
        if (group_json.type != JsonType::Object) {
            throw CliError("each JSON group must be an object");
        }
        const JsonValue *data = json_field(group_json, "data");
        if (data == nullptr) {
            throw CliError("each JSON group needs a data array");
        }
        std::vector<int> group_data = parse_int_array(*data, "groups.data");
        for (int symbol : group_data) {
            if (symbol < 0 || symbol >= code->data) {
                throw CliError("group data symbol is out of range");
            }
            if (seen_data[static_cast<std::size_t>(symbol)] != 0) {
                throw CliError("group data symbols must be disjoint");
            }
            seen_data[static_cast<std::size_t>(symbol)] = 1;
        }

        const JsonValue *local = json_field(group_json, "local");
        if (local != nullptr) {
            std::vector<int> local_rows = parse_int_array(*local, "groups.local");
            if (local_rows.size() !=
                static_cast<std::size_t>(code->groups[group_index].local_parity)) {
                throw CliError("group local row count does not match local_parity");
            }
            for (std::size_t local_index = 0; local_index < local_rows.size(); local_index++) {
                if (local_rows[local_index] != next_local_row + static_cast<int>(local_index)) {
                    throw CliError("JSON local rows must be contiguous after data rows");
                }
            }
        }

        code->groups[group_index].data = std::move(group_data);
        code->groups[group_index].local_row_start = next_local_row;
        next_local_row += code->groups[group_index].local_parity;
    }

    for (uint8_t seen : seen_data) {
        if (seen == 0) {
            throw CliError("JSON groups must cover every data symbol exactly once");
        }
    }
}

std::vector<uint8_t> parse_json_matrix(const JsonValue &root, int expected_rows, int expected_cols)
{
    bool hex_context = false;
    const JsonValue *matrix = json_field(root, "matrix");
    if (matrix == nullptr) {
        matrix = json_field(root, "matrix_hex");
        hex_context = true;
    }
    if (matrix == nullptr) {
        throw CliError("JSON input needs a matrix or matrix_hex field");
    }
    if (matrix->type != JsonType::Array) {
        throw CliError("JSON matrix must be an array of rows");
    }
    if (matrix->array.size() != static_cast<std::size_t>(expected_rows)) {
        throw CliError("JSON matrix row count does not match data/local/global parameters");
    }

    std::vector<uint8_t> values;
    values.reserve(static_cast<std::size_t>(expected_rows) * static_cast<std::size_t>(expected_cols));
    for (const JsonValue &row : matrix->array) {
        if (row.type != JsonType::Array) {
            throw CliError("JSON matrix rows must be arrays");
        }
        if (row.array.size() != static_cast<std::size_t>(expected_cols)) {
            throw CliError("JSON matrix column count does not match data_cnt");
        }
        for (const JsonValue &entry : row.array) {
            values.push_back(parse_json_byte(entry, "matrix", hex_context));
        }
    }
    return values;
}

void validate_systematic_data_local_matrix(const mrlrc::Code &code)
{
    for (int row = 0; row < code.data; row++) {
        for (int col = 0; col < code.data; col++) {
            uint8_t expected = row == col ? 1 : 0;
            if (code.matrix[static_cast<std::size_t>(row) *
                                static_cast<std::size_t>(code.data) +
                            static_cast<std::size_t>(col)] != expected) {
                throw CliError("JSON matrix must have identity data rows");
            }
        }
    }

    std::vector<uint8_t> in_group(static_cast<std::size_t>(code.data), 0);
    for (const auto &group : code.groups) {
        std::fill(in_group.begin(), in_group.end(), 0);
        for (int data : group.data) {
            in_group[static_cast<std::size_t>(data)] = 1;
        }
        for (int local = 0; local < group.local_parity; local++) {
            int row = group.local_row_start + local;
            for (int col = 0; col < code.data; col++) {
                if (in_group[static_cast<std::size_t>(col)] == 0 &&
                    code.matrix[static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(code.data) +
                                static_cast<std::size_t>(col)] != 0) {
                    throw CliError("JSON local rows must only cover data symbols in their group");
                }
            }
        }
    }
}

mrlrc::Code read_code_json(const std::string &path)
{
    JsonParser parser(read_text_file(path));
    JsonValue root = parser.parse();
    if (root.type != JsonType::Object) {
        throw CliError("JSON input root must be an object");
    }

    const JsonValue *groups_json = json_field(root, "groups");
    int data = required_json_int(root, {"data", "data_cnt"});
    int groups = optional_json_int(root, {"group_cnt"}, std::numeric_limits<int>::min());
    if (groups == std::numeric_limits<int>::min() && groups_json != nullptr &&
        groups_json->type == JsonType::Array) {
        groups = static_cast<int>(groups_json->array.size());
    }
    int local_parity =
        optional_json_int(root, {"local_parity", "local_parity_cnt"}, std::numeric_limits<int>::min());
    if (local_parity == std::numeric_limits<int>::min() && groups_json != nullptr &&
        groups_json->type == JsonType::Array && !groups_json->array.empty()) {
        const JsonValue *local = json_field(groups_json->array.front(), "local");
        if (local != nullptr && local->type == JsonType::Array) {
            local_parity = static_cast<int>(local->array.size());
        }
    }
    int global_parity = required_json_int(root, {"global_parity", "global_parity_cnt"});
    if (data == std::numeric_limits<int>::min()) {
        throw CliError("JSON input needs data or data_cnt");
    }
    if (groups == std::numeric_limits<int>::min()) {
        throw CliError("JSON input needs group_cnt or groups");
    }
    if (local_parity == std::numeric_limits<int>::min()) {
        throw CliError("JSON input needs local_parity/local_parity_cnt or groups.local");
    }
    if (global_parity == std::numeric_limits<int>::min()) {
        throw CliError("JSON input needs global_parity or global_parity_cnt");
    }

    mrlrc::Params params;
    params.data = data;
    params.groups = groups;
    params.local_parity = local_parity;
    params.global_parity = global_parity;
    params.local_family = mrlrc::MatrixFamily::Random;
    params.global_family = mrlrc::MatrixFamily::Random;
    params.construction = optional_json_uint64(root, "construction", 0);

    mrlrc::Code code = mrlrc::make_code_layout(params);
    code.candidate_source = "json";
    std::string local_method = optional_json_string(root, "local_method", "");
    if (!local_method.empty()) {
        if (!mrlrc::parse_family(local_method, &code.local_family)) {
            throw CliError("invalid local_method in JSON: " + local_method);
        }
    }
    std::string global_method = optional_json_string(root, "global_method", "");
    if (!global_method.empty()) {
        if (!mrlrc::parse_family(global_method, &code.global_family)) {
            throw CliError("invalid global_method in JSON: " + global_method);
        }
    }
    if (groups_json != nullptr) {
        apply_json_groups(*groups_json, &code);
    }

    code.matrix = parse_json_matrix(root, code.symbols, code.data);
    validate_systematic_data_local_matrix(code);
    return code;
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

void print_prefilter_stats(const mrlrc::GenerateResult &result, uint64_t prefilter_count)
{
    if (!result.prefilter_enabled) {
        return;
    }
    print_attribute("prefilter_count", prefilter_count);
    print_attribute("prefilter_candidates_checked", result.prefilter_candidates_checked);
    print_attribute("prefilter_candidates_rejected", result.prefilter_candidates_rejected);
    print_attribute("prefilter_patterns_checked", result.prefilter_patterns_checked);
}

int run_check_json(const std::string &path)
{
    mrlrc::Code code = read_code_json(path);
    mrlrc::CheckResult check = mrlrc::check_mr(code);
    code.patterns_checked = check.patterns_checked;

    if (check.is_mr) {
        std::cout << "MR-LRC matrix verified\n";
        print_attribute("data", code.data);
        print_attribute("groups", code.groups.size());
        print_attribute("local_rows", code.local_rows);
        print_attribute("global_parity", code.global_parity);
        print_attribute("total_parity", code.total_parity);
        print_attribute("symbols", code.symbols);
        print_attribute("candidate_source", code.candidate_source);
        print_attribute("patterns_checked", check.patterns_checked);
        print_attribute("strict_complete", check.strict_complete ? "true" : "false");
        print_attribute("gf256_backend", mrlrc::gf256_backend());
        return 0;
    }

    std::cout << "MR-LRC matrix rejected\n";
    print_attribute("message", check.message);
    print_attribute("patterns_checked", check.patterns_checked);
    print_attribute("failures", check.failures);
    print_attribute("strict_complete", check.strict_complete ? "true" : "false");
    print_failed_pattern(code, check.first_failed_erased);
    return 2;
}

} // namespace

int main(int argc, char **argv)
{
    try {
        CliOptions opts = parse_args(argc, argv);
        if (!opts.check_json_path.empty()) {
            int status = run_check_json(opts.check_json_path);
            gflags::ShutDownCommandLineFlags();
            return status;
        }

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
                print_prefilter_stats(result, opts.params.prefilter_count);
                print_failed_pattern(code, result.check.first_failed_erased);
            }
            gflags::ShutDownCommandLineFlags();
            return result.status == 2 ? 2 : 1;
        }

        if (!opts.json_path.empty()) {
            write_matrix_json_file(opts.json_path, result, opts.params);
        }

        std::cout << "MR-LRC matrix found\n";
        print_attribute("data", code.data);
        print_attribute("groups", code.groups.size());
        print_attribute("local_parity", opts.params.local_parity);
        print_attribute("local_rows", code.local_rows);
        print_attribute("global_parity", code.global_parity);
        print_attribute("total_parity", code.total_parity);
        print_attribute("symbols", code.symbols);
        print_attribute("construction", code.construction);
        print_attribute("candidate_source", code.candidate_source);
        print_attribute("local_method", mrlrc::family_name(code.local_family));
        print_attribute("global_method", mrlrc::family_name(code.global_family));
        print_attribute("seed", code.seed);
        print_attribute("random_limit", opts.params.random_limit);
        print_attribute("thread_count", opts.params.thread_count);
        print_attribute("attempt", code.attempt);
        print_attribute("patterns_checked", result.check.patterns_checked);
        print_attribute("strict_complete",
                        result.check.strict_complete ? "true" : "false");
        print_prefilter_stats(result, opts.params.prefilter_count);
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
