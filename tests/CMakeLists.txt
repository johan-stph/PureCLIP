# =============================================================================
# tests/CMakeLists.txt
#
# CTest-based regression test suite for PureCLIP.
#
# Labels
# ------
#   setup      – data-generation helpers (not real tests)
#   tier1      – fast synthetic tests  (~1-5 s each, <256 MB)
#   tier2      – real chrM tests       (~30-60 s each, <512 MB)
#   regression – exact output comparison against golden files
#   mode       – smoke-only mode/flag variants
#   determinism– threading consistency checks
#
# Typical usage:
#   ctest -L tier1                   # quick sanity check
#   ctest -L "tier1|tier2"           # full regression
#   ctest -L tier1 --output-on-failure
#
# Golden files:
#   Missing golden files cause the test to SKIP (exit 77), not FAIL.
#   Generate them once with:  tests/generate_golden.sh <build>/pureclip
#   Then commit tests/golden/.
# =============================================================================

set(T           "${CMAKE_CURRENT_SOURCE_DIR}")   # tests/ root
set(SYN         "${T}/data/synthetic")
set(REAL        "${T}/data/real")
set(GOLDEN_SYN  "${T}/golden/synthetic")
set(GOLDEN_CHR  "${T}/golden/chrM")
set(OUT         "${CMAKE_BINARY_DIR}/test_output")
set(PY          "python3")
set(CMP         "${T}/compare.py")
set(BIN         "$<TARGET_FILE:pureclip>")

file(MAKE_DIRECTORY "${OUT}/synthetic")
file(MAKE_DIRECTORY "${OUT}/chrM")

# =============================================================================
# Tier 1 – Synthetic
# =============================================================================

# T1.0  generate synthetic data (setup step, always runs first) ──────────────
add_test(NAME t1_0_make_data
    COMMAND ${PY} "${T}/make_synthetic.py")
set_tests_properties(t1_0_make_data PROPERTIES
    LABELS   "tier1;setup"
    TIMEOUT  60)

# Helper macro: standard synthetic pureclip invocation ───────────────────────
# Usage: add_syn_run(<test_name> <out_prefix> <extra_args...>)
macro(add_syn_run TEST_NAME OUT_PFX)
    add_test(NAME ${TEST_NAME}
        COMMAND ${BIN}
            -i   "${SYN}/sample.bam"
            -bai "${SYN}/sample.bam.bai"
            -g   "${SYN}/ref.fa"
            -o   "${OUT}/synthetic/${OUT_PFX}_sites.bed"
            -or  "${OUT}/synthetic/${OUT_PFX}_regions.bed"
            -p   "${OUT}/synthetic/${OUT_PFX}_params.txt"
            -nt  1
            ${ARGN})
    set_tests_properties(${TEST_NAME} PROPERTIES
        LABELS   "tier1"
        TIMEOUT  120
        DEPENDS  t1_0_make_data)
endmacro()

# T1.1  smoke – exit 0, output files are non-empty ───────────────────────────
add_syn_run(t1_1_smoke smoke)
# Verify files are non-empty in a follow-up step
add_test(NAME t1_1_smoke_check
    COMMAND ${PY} "${CMP}" count "${OUT}/synthetic/smoke_sites.bed" 1)
set_tests_properties(t1_1_smoke_check PROPERTIES
    LABELS   "tier1"
    TIMEOUT  10
    DEPENDS  t1_1_smoke)

# T1.2  sites regression – exact BED match ────────────────────────────────────
add_syn_run(t1_2_reg_run reg)
add_test(NAME t1_2_regression_sites
    COMMAND ${PY} "${CMP}" bed
        "${OUT}/synthetic/reg_sites.bed"
        "${GOLDEN_SYN}/sites.bed")
set_tests_properties(t1_2_regression_sites PROPERTIES
    LABELS          "tier1;regression"
    TIMEOUT         15
    DEPENDS         t1_2_reg_run
    SKIP_RETURN_CODE 77)

# T1.3  regions regression ────────────────────────────────────────────────────
add_test(NAME t1_3_regression_regions
    COMMAND ${PY} "${CMP}" bed
        "${OUT}/synthetic/reg_regions.bed"
        "${GOLDEN_SYN}/regions.bed")
set_tests_properties(t1_3_regression_regions PROPERTIES
    LABELS          "tier1;regression"
    TIMEOUT         15
    DEPENDS         t1_2_reg_run
    SKIP_RETURN_CODE 77)

# T1.4  params regression – float-tolerant ────────────────────────────────────
add_test(NAME t1_4_regression_params
    COMMAND ${PY} "${CMP}" params
        "${OUT}/synthetic/reg_params.txt"
        "${GOLDEN_SYN}/params.txt"
        --tol 1e-4)
set_tests_properties(t1_4_regression_params PROPERTIES
    LABELS          "tier1;regression"
    TIMEOUT         15
    DEPENDS         t1_2_reg_run
    SKIP_RETURN_CODE 77)

# T1.5  scoring types 1-3 (code-path coverage, smoke only) ───────────────────
foreach(ST 1 2 3)
    add_syn_run(t1_5_st${ST} st${ST} -st ${ST})
endforeach()

# T1.6  binding characteristics mode 1 ───────────────────────────────────────
add_syn_run(t1_6_bc1 bc1 -bc 1)

# T1.7  output-all mode ───────────────────────────────────────────────────────
add_syn_run(t1_7_outputall oa -oa)
add_test(NAME t1_7_outputall_check
    COMMAND ${PY} "${CMP}" count "${OUT}/synthetic/oa_sites.bed" 1)
set_tests_properties(t1_7_outputall_check PROPERTIES
    LABELS   "tier1;mode"
    TIMEOUT  10
    DEPENDS  t1_7_outputall)

# T1.8  thread determinism – nt=4 must produce identical BED ─────────────────
add_syn_run(t1_8_nt4 nt4 -nt 4)
add_test(NAME t1_8_thread_determinism
    COMMAND ${PY} "${CMP}" bed
        "${OUT}/synthetic/nt4_sites.bed"
        "${GOLDEN_SYN}/sites.bed")
set_tests_properties(t1_8_thread_determinism PROPERTIES
    LABELS          "tier1;determinism"
    TIMEOUT         15
    DEPENDS         t1_8_nt4
    SKIP_RETURN_CODE 77)


# =============================================================================
# Tier 2 – Real chrM data
# Conditionally registered: skipped at cmake-configure time if data is absent.
# =============================================================================

if(EXISTS "${REAL}/chrM.bam" AND EXISTS "${REAL}/chrM.fa")

    # Helper macro: standard chrM pureclip invocation ────────────────────────
    macro(add_chrm_run TEST_NAME OUT_PFX)
        add_test(NAME ${TEST_NAME}
            COMMAND ${BIN}
                -i   "${REAL}/chrM.bam"
                -bai "${REAL}/chrM.bam.bai"
                -g   "${REAL}/chrM.fa"
                -o   "${OUT}/chrM/${OUT_PFX}_sites.bed"
                -or  "${OUT}/chrM/${OUT_PFX}_regions.bed"
                -p   "${OUT}/chrM/${OUT_PFX}_params.txt"
                -nt  4
                ${ARGN})
        set_tests_properties(${TEST_NAME} PROPERTIES
            LABELS   "tier2"
            TIMEOUT  300)
    endmacro()

    # T2.1  smoke ─────────────────────────────────────────────────────────────
    add_chrm_run(t2_1_smoke_chrM smoke)
    add_test(NAME t2_1_smoke_chrM_check
        COMMAND ${PY} "${CMP}" count "${OUT}/chrM/smoke_sites.bed" 1)
    set_tests_properties(t2_1_smoke_chrM_check PROPERTIES
        LABELS   "tier2"
        TIMEOUT  15
        DEPENDS  t2_1_smoke_chrM)

    # T2.2  sites regression ──────────────────────────────────────────────────
    add_chrm_run(t2_2_reg_run reg)
    add_test(NAME t2_2_regression_chrM_sites
        COMMAND ${PY} "${CMP}" bed
            "${OUT}/chrM/reg_sites.bed"
            "${GOLDEN_CHR}/sites.bed")
    set_tests_properties(t2_2_regression_chrM_sites PROPERTIES
        LABELS          "tier2;regression"
        TIMEOUT         15
        DEPENDS         t2_2_reg_run
        SKIP_RETURN_CODE 77)

    # T2.3  regions regression ────────────────────────────────────────────────
    add_test(NAME t2_3_regression_chrM_regions
        COMMAND ${PY} "${CMP}" bed
            "${OUT}/chrM/reg_regions.bed"
            "${GOLDEN_CHR}/regions.bed")
    set_tests_properties(t2_3_regression_chrM_regions PROPERTIES
        LABELS          "tier2;regression"
        TIMEOUT         15
        DEPENDS         t2_2_reg_run
        SKIP_RETURN_CODE 77)

    # T2.4  params regression ─────────────────────────────────────────────────
    add_test(NAME t2_4_regression_chrM_params
        COMMAND ${PY} "${CMP}" params
            "${OUT}/chrM/reg_params.txt"
            "${GOLDEN_CHR}/params.txt"
            --tol 1e-4)
    set_tests_properties(t2_4_regression_chrM_params PROPERTIES
        LABELS          "tier2;regression"
        TIMEOUT         15
        DEPENDS         t2_2_reg_run
        SKIP_RETURN_CODE 77)

    # T2.5  bc=1 mode on real data ────────────────────────────────────────────
    add_chrm_run(t2_5_chrM_bc1 bc1 -bc 1)

else()
    message(STATUS
        "Tier 2 tests DISABLED: chrM data not found.\n"
        "       Run:  tests/prepare_real.sh <path/to/sample_run>")
endif()
