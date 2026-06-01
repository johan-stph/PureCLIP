// ======================================================================
// PureCLIP — minimal unit tests for GAMMA / ZTBIN / LogSumExp
//
// Compiled as a standalone executable sharing PureCLIP's build infra.
// Includes the same headers as pureclip.cpp to get all SeqAn types.
// ======================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

// Include PureCLIP headers in the same order as pureclip.cpp
#include "version.h"
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/misc/name_store_cache.h>
#include <seqan/arg_parse.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

#include "util.h"
#include "call_sites.h"
#include "call_sites_replicates.h"

using namespace seqan;

#include <cmath>
#include <cassert>

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) \
    do { std::cout << "  " << #name << "... "; } while(0)

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::cerr << "FAIL at line " << __LINE__ << ": " << #cond << std::endl; \
            tests_failed++; \
            return; \
        } \
    } while(0)

#define CHECK_CLOSE(a, b, tol) \
    do { \
        if (std::fabs((double)(a) - (double)(b)) > (double)(tol)) { \
            std::cerr << "FAIL: " << #a << "=" << (a) << " vs " << #b << "=" << (b) \
                      << " (tol=" << tol << ")" << std::endl; \
            tests_failed++; \
            return; \
        } \
    } while(0)

#define PASS() \
    do { std::cout << "OK" << std::endl; tests_passed++; } while(0)

// =============================================================================
// GAMMA::getDensity tests
// =============================================================================

void test_gamma_density_at_mode() {
    TEST(gamma_density_at_mode);
    GAMMA g;
    g.k = 2.0;
    g.b0 = log(2.0);  // theta = exp(b0)/k = 1.0, mean = k*theta = 2.0
    g.tp = 0.0;

    // For k=2, theta=1, gamma PDF: x * exp(-x)
    // Mode at (k-1)*theta = 1.0, PDF = exp(-1) ≈ 0.367879
    Float d = g.getDensity(1.0);
    Float expected = exp(static_cast<Float>(-1.0));
    CHECK_CLOSE(d, expected, 1e-6);
    PASS();
}

void test_gamma_density_below_truncation() {
    TEST(gamma_density_below_truncation);
    GAMMA g;
    g.k = 2.0;
    g.b0 = log(2.0);
    g.tp = 0.5;

    Float d = g.getDensity(0.3);
    CHECK(d == 0.0);
    PASS();
}

void test_gamma_density_above_truncation() {
    TEST(gamma_density_above_truncation);
    GAMMA g;
    g.k = 2.0;
    g.b0 = log(2.0);
    g.tp = 0.5;

    Float d = g.getDensity(1.0);
    CHECK(d > 0.0);
    PASS();
}

void test_gamma_cache_invalidation() {
    TEST(gamma_cache_invalidation);
    GAMMA g;
    g.k = 2.0;
    g.b0 = log(2.0);
    g.tp = 0.0;

    Float d1 = g.getDensity(1.0);
    g.k = 3.0;
    g.invalidateCache();
    Float d2 = g.getDensity(1.0);
    CHECK(std::fabs(d1 - d2) > 1e-10);
    PASS();
}

void test_gamma_small_values() {
    TEST(gamma_small_values);
    GAMMA g;
    g.k = 1.5;
    g.b0 = log(0.01);
    g.tp = 0.0001;

    Float d = g.getDensity(0.0002);
    CHECK(!std::isnan(d));
    CHECK(!std::isinf(d));
    CHECK(d >= 0.0);
    PASS();
}

// =============================================================================
// ZTBIN::getDensity tests
// =============================================================================

void test_ztbin_zero_k() {
    TEST(ztbin_zero_k);
    ZTBIN bin;
    bin.p = 0.15;
    AppOptions opts;
    Float d = bin.getDensity(0, 10, opts);
    CHECK(d == 0.0);
    PASS();
}

void test_ztbin_k_equals_n() {
    TEST(ztbin_k_equals_n);
    ZTBIN bin;
    bin.p = 0.15;
    AppOptions opts;
    Float d = bin.getDensity(5, 5, opts);
    CHECK(d > 0.0);
    CHECK(d <= 1.0);
    PASS();
}

void test_ztbin_probability_sums() {
    TEST(ztbin_probability_sums);
    ZTBIN bin;
    bin.p = 0.3;
    AppOptions opts;
    unsigned n = 10;
    Float sum = 0.0;
    for (unsigned k = 1; k <= n; ++k) {
        sum += bin.getDensity(k, n, opts);
    }
    CHECK_CLOSE(sum, 1.0, 1e-6);
    PASS();
}

// =============================================================================
// LogSumExp tests
// =============================================================================

void test_logsumexp_identity() {
    TEST(logsumexp_identity);
    LogSumExp_lookupTable lut(1000, -20.0);

    Float a = 0.5;
    Float b = 0.5;
    Float result = lut.logSumExp_add(a, b);
    Float expected = 0.5 + log(2.0);
    CHECK_CLOSE(result, expected, 1e-6);
    PASS();
}

void test_logsumexp_large_difference() {
    TEST(logsumexp_large_difference);
    LogSumExp_lookupTable lut(1000, -20.0);

    Float result = lut.logSumExp_add(10.0, -10.0);
    // log(exp(10) + exp(-10)) = 10 + log(1+exp(-20)) ≈ 10 + 2.06e-9
    // Lookup table interpolation adds ~step/2 error
    CHECK_CLOSE(result, 10.0, 1e-7);
    PASS();
}

void test_logsumexp_small_difference() {
    TEST(logsumexp_small_difference);
    LogSumExp_lookupTable lut(1000, -20.0);

    Float a = 5.0;
    Float b = 0.0;
    Float result = lut.logSumExp_add(a, b);
    Float expected = log(exp(5.0) + exp(0.0));
    CHECK_CLOSE(result, expected, 1e-6);
    PASS();
}

void test_logsumexp_below_precision() {
    TEST(logsumexp_below_precision);
    LogSumExp_lookupTable lut(1000, -20.0);

    Float result = lut.logSumExp_add(5.0, -30.0);
    CHECK_CLOSE(result, 5.0, 1e-15);
    PASS();
}

// =============================================================================
// myLog / myExp tests
// =============================================================================

void test_mylog_zero() {
    TEST(mylog_zero);
    Float result = myLog(0.0);
    CHECK(std::isnan(result));
    PASS();
}

void test_myexp_nan() {
    TEST(myexp_nan);
    Float result = myExp(std::numeric_limits<Float>::quiet_NaN());
    CHECK(result == 0.0);
    PASS();
}

// =============================================================================

int main() {
    std::cout << "=== Unit tests: GAMMA ===" << std::endl;
    test_gamma_density_at_mode();
    test_gamma_density_below_truncation();
    test_gamma_density_above_truncation();
    test_gamma_cache_invalidation();
    test_gamma_small_values();

    std::cout << "=== Unit tests: ZTBIN ===" << std::endl;
    test_ztbin_zero_k();
    test_ztbin_k_equals_n();
    test_ztbin_probability_sums();

    std::cout << "=== Unit tests: LogSumExp ===" << std::endl;
    test_logsumexp_identity();
    test_logsumexp_large_difference();
    test_logsumexp_small_difference();
    test_logsumexp_below_precision();

    std::cout << "=== Unit tests: myLog/myExp ===" << std::endl;
    test_mylog_zero();
    test_myexp_nan();

    std::cout << std::endl;
    std::cout << "Results: " << tests_passed << " passed, " 
              << tests_failed << " failed" << std::endl;
    return tests_failed > 0 ? 1 : 0;
}
