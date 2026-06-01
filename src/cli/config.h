// ======================================================================
// PureCLIP — TOML configuration file loader
// ======================================================================
// Loads a TOML config into AppOptions.  CLI flags override TOML values.
// Unknown keys are silently ignored for forward compatibility.
// ======================================================================

#ifndef PURECLIP_CONFIG_H_
#define PURECLIP_CONFIG_H_

#include "toml.hpp"
#include "util.h"

#include <string>
#include <iostream>
#include <fstream>

using namespace seqan;

// Helper: set a scalar if the key exists in the TOML table
template<typename T>
void set_if_present(const toml::table& tbl, const char* key, T& dest) {
    auto node = tbl.get(key);
    if (node && node->is<T>()) {
        dest = node->value<T>().value_or(dest);
    }
}

// Overload for unsigned (TOML uses int64_t)
inline void set_if_present_uint(const toml::table& tbl, const char* key, unsigned& dest) {
    auto node = tbl.get(key);
    if (node && node->is_integer()) {
        auto val = node->value<int64_t>();
        if (val && *val >= 0) dest = static_cast<unsigned>(*val);
    }
}

// Helper: report an error and return false
inline bool config_error(const std::string& msg) {
    std::cerr << "ERROR: config file: " << msg << std::endl;
    return false;
}

// =============================================================================
// Load TOML config into AppOptions.  Returns true on success.
// Does NOT override values already set via CLI — call before parseCommandLine
// and the CLI parser will naturally override.
// =============================================================================
inline bool load_toml_config(const std::string& path, AppOptions& options) {
    // Parse the file
    toml::table tbl;
    try {
        tbl = toml::parse_file(path);
    } catch (const toml::parse_error& e) {
        return config_error(e.what());
    }

    // ── Required: input data ────────────────────────────────────────────
    if (auto bam = tbl.get("bam")) {
        if (auto arr = bam->as_array()) {
            clear(options.bamFileNames);
            for (const auto& v : *arr) {
                if (auto s = v.value<std::string>())
                    appendValue(options.bamFileNames, s->c_str());
            }
        }
    }
    if (auto bai = tbl.get("bai")) {
        if (auto arr = bai->as_array()) {
            clear(options.baiFileNames);
            for (const auto& v : *arr) {
                if (auto s = v.value<std::string>())
                    appendValue(options.baiFileNames, s->c_str());
            }
        }
    }
    if (auto g = tbl.get("genome")) {
        if (auto s = g->value<std::string>()) options.refFileName = s->c_str();
    }

    // ── Output prefix ───────────────────────────────────────────────────
    if (auto o = tbl.get("output_prefix")) {
        if (auto s = o->value<std::string>()) {
            std::string prefix = *s;
            if (options.outFileName == "") {
                options.outFileName = (prefix + "_sites.bed").c_str();
            }
            if (options.outRegionsFileName == "") {
                options.outRegionsFileName = (prefix + "_regions.bed").c_str();
            }
            if (options.parFileName == "") {
                options.parFileName = (prefix + "_params.txt").c_str();
            }
        }
    }

    // ── Genomic subsets ─────────────────────────────────────────────────
    if (auto l = tbl.get("learn_on")) {
        if (auto s = l->value<std::string>()) options.intervals_str = s->c_str();
    }
    if (auto a = tbl.get("apply_on")) {
        if (auto s = a->value<std::string>()) options.applyChr_str = s->c_str();
    }

    // ── Performance ─────────────────────────────────────────────────────
    set_if_present_uint(tbl, "threads", options.numThreads);
    set_if_present_uint(tbl, "threads_apply", options.numThreadsA);

    // ── Bandwidths ──────────────────────────────────────────────────────
    set_if_present_uint(tbl, "bandwidth", options.bandwidth);
    set_if_present_uint(tbl, "bandwidth_n", options.bandwidthN);

    // ── Binding mode ────────────────────────────────────────────────────
    if (auto bc = tbl.get("binding_characteristics")) {
        if (auto val = bc->value<int64_t>()) {
            if (*val == 1) {
                options.bandwidthN = 100;
                options.get_nThreshold = true;
                options.p1 = 0.01;
                options.p2 = 0.1;
            }
        }
    }

    // ── Initial values ──────────────────────────────────────────────────
    set_if_present(tbl, "bin1_p_init", options.p1);
    set_if_present(tbl, "bin2_p_init", options.p2);

    // ── Thresholds ──────────────────────────────────────────────────────
    set_if_present(tbl, "min_trans_prob_crosslink", options.minTransProbCS);
    set_if_present_uint(tbl, "prior_enrichment_threshold", options.prior_enrichmentThreshold);

    // ── Convergence ─────────────────────────────────────────────────────
    set_if_present_uint(tbl, "max_iter_baumwelch", options.maxIter_bw);
    set_if_present_uint(tbl, "max_iter_simplex", options.maxIter_simplex);
    set_if_present(tbl, "gamma_k_convergence", options.gamma_k_conv);
    set_if_present(tbl, "gamma_theta_convergence", options.gamma_theta_conv);
    set_if_present(tbl, "bin_p_convergence", options.bin_p_conv);

    // ── Gamma constraints ───────────────────────────────────────────────
    set_if_present(tbl, "g1_k_min", options.g1_kMin);
    set_if_present(tbl, "g1_k_max", options.g1_kMax);
    set_if_present(tbl, "g2_k_min", options.g2_kMin);
    set_if_present(tbl, "g2_k_max", options.g2_kMax);

    // ── Scoring / output ────────────────────────────────────────────────
    set_if_present_uint(tbl, "scoring_scheme", options.score_type);
    set_if_present_uint(tbl, "merge_distance", options.distMerge);

    // ── Debug ───────────────────────────────────────────────────────────
    if (auto v = tbl.get("verbosity")) {
        if (auto val = v->value<int64_t>()) options.verbosity = static_cast<int>(*val);
    }
    if (auto hp = tbl.get("high_precision")) {
        if (auto val = hp->value<bool>()) options.useHighPrecision = *val;
    }

    // ── Covariates ──────────────────────────────────────────────────────
    if (auto is_ = tbl.get("input_signal_file")) {
        if (auto s = is_->value<std::string>()) options.rpkmFileName = s->c_str();
    }
    if (auto ib = tbl.get("input_bam_file")) {
        if (auto s = ib->value<std::string>()) options.inputBamFileName = s->c_str();
    }
    if (auto ibai = tbl.get("input_bai_file")) {
        if (auto s = ibai->value<std::string>()) options.inputBaiFileName = s->c_str();
    }
    if (auto fi = tbl.get("fimo_file")) {
        if (auto s = fi->value<std::string>()) options.fimoFileName = s->c_str();
    }

    // Validate: if input_bam is set, input_bai must also be set
    if (options.inputBamFileName != "" && options.inputBaiFileName == "") {
        return config_error("input_bam_file set but input_bai_file not set");
    }
    if (options.inputBamFileName == "" && options.inputBaiFileName != "") {
        return config_error("input_bai_file set but input_bam_file not set");
    }

    if (options.verbosity >= 1) {
        std::cout << "Loaded configuration from " << path << std::endl;
    }

    return true;
}

#endif // PURECLIP_CONFIG_H_
