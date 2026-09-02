# =============================================================================
# tests/CMakeLists.txt — correctness regression suite for PureCLIP
#
# Labels
#   tier1   fast synthetic fixture   (~1 s total)   – run on every change
#   tier2   real chrM data           (~15 s total)  – run before pushing
#   truth   scored against implanted crosslinks, not against previous output
#
# Usage
#   ctest --test-dir build -L tier1 --output-on-failure
#   ctest --test-dir build --output-on-failure
#
# A missing golden file makes a test SKIP (not fail), so the suite is usable
# on a platform whose goldens have not been generated yet. See TESTING.md.
# =============================================================================

set(T          "${CMAKE_CURRENT_SOURCE_DIR}")
set(SYN        "${T}/data/synthetic")
set(CHRM       "${T}/data/chrM")
set(GOLD_SYN   "${T}/golden/synthetic/${CMAKE_SYSTEM_NAME}")
set(GOLD_CHRM  "${T}/golden/chrM/${CMAKE_SYSTEM_NAME}")
set(GOLD_COV   "${T}/golden/chrM_cov/${CMAKE_SYSTEM_NAME}")
set(OUT        "${CMAKE_BINARY_DIR}/test_output")
set(CMP        "${T}/compare.py")
set(PURECLIP   "$<TARGET_FILE:pureclip>")

file(MAKE_DIRECTORY "${OUT}")

# --- helpers ----------------------------------------------------------------
# Register a PureCLIP run producing ${OUT}/<tag>.{sites,regions,params}
function(add_run tag tier data bam fa)
    add_test(NAME ${tier}_run_${tag}
             COMMAND ${PURECLIP}
                     -i   "${data}/${bam}"
                     -bai "${data}/${bam}.bai"
                     -g   "${data}/${fa}"
                     -o   "${OUT}/${tag}.sites.bed"
                     -or  "${OUT}/${tag}.regions.bed"
                     -p   "${OUT}/${tag}.params.txt"
                     -nt  4)
    set_tests_properties(${tier}_run_${tag} PROPERTIES
        LABELS "${tier};run" FIXTURES_SETUP ${tag}_output TIMEOUT 900)
endfunction()

# Register a compare.py check that consumes a run's output
function(add_check name tier tag)
    add_test(NAME ${name} COMMAND python3 "${CMP}" ${ARGN})
    set_tests_properties(${name} PROPERTIES
        LABELS "${tier}" FIXTURES_REQUIRED ${tag}_output
        SKIP_RETURN_CODE 77 TIMEOUT 120)
endfunction()

# --- tier 1: synthetic ------------------------------------------------------
add_test(NAME tier1_make_data COMMAND python3 "${T}/make_synthetic.py")
set_tests_properties(tier1_make_data PROPERTIES
    LABELS "tier1;setup" FIXTURES_SETUP synthetic_data TIMEOUT 300)

add_run(syn tier1 "${SYN}" sample.bam ref.fa)
set_tests_properties(tier1_run_syn PROPERTIES FIXTURES_REQUIRED synthetic_data)

add_check(tier1_smoke        tier1 syn count  "${OUT}/syn.sites.bed" 1)
add_check(tier1_sites        tier1 syn bed    "${OUT}/syn.sites.bed"   "${GOLD_SYN}/sites.bed")
add_check(tier1_regions      tier1 syn bed    "${OUT}/syn.regions.bed" "${GOLD_SYN}/regions.bed")
add_check(tier1_params       tier1 syn params "${OUT}/syn.params.txt"  "${GOLD_SYN}/params.txt" --tol 1e-4)
add_check(tier1_agreement    tier1 syn agree  "${OUT}/syn.sites.bed"   "${GOLD_SYN}/sites.bed" --min-jaccard 0.95)

# Ground truth: does PureCLIP find the crosslinks we implanted? This is the one
# check that stays meaningful if the golden files themselves are wrong.
add_check(tier1_ground_truth tier1 syn recall "${OUT}/syn.sites.bed"
          "${SYN}/truth.bed" --window 2 --min-recall 1.0)
set_tests_properties(tier1_ground_truth PROPERTIES LABELS "tier1;truth")

# Output must not depend on the thread count.
add_test(NAME tier1_determinism
         COMMAND ${CMAKE_COMMAND}
                 -DPURECLIP=${PURECLIP} -DSYN=${SYN} -DOUT=${OUT} -DCMP=${CMP}
                 -P "${T}/determinism.cmake")
set_tests_properties(tier1_determinism PROPERTIES
    LABELS "tier1;determinism" FIXTURES_REQUIRED synthetic_data TIMEOUT 600)

# --- tier 2: chrM (real data) -----------------------------------------------
if(EXISTS "${CHRM}/chrM.bam")
    add_run(chrm tier2 "${CHRM}" chrM.bam chrM.fa)
    add_check(tier2_smoke     tier2 chrm count  "${OUT}/chrm.sites.bed" 100)
    add_check(tier2_sites     tier2 chrm bed    "${OUT}/chrm.sites.bed"   "${GOLD_CHRM}/sites.bed")
    add_check(tier2_regions   tier2 chrm bed    "${OUT}/chrm.regions.bed" "${GOLD_CHRM}/regions.bed")
    add_check(tier2_params    tier2 chrm params "${OUT}/chrm.params.txt"  "${GOLD_CHRM}/params.txt" --tol 1e-4)
    add_check(tier2_agreement tier2 chrm agree  "${OUT}/chrm.sites.bed"   "${GOLD_CHRM}/sites.bed" --min-jaccard 0.95)
endif()

# --- tier 2: covariate model (-ibam) ----------------------------------------
# Exercises GAMMA_REG and the seven-argument parse_bamRegion, neither of which
# the default invocation touches. That gap let a crash ship: the +1000 jump
# lookahead ran past the contig end and SeqAn 2.5.3 aborted the run.
if(EXISTS "${CHRM}/chrM.bam")
    add_test(NAME tier2_run_cov
             COMMAND ${PURECLIP}
                     -i    "${CHRM}/chrM.bam"  -bai  "${CHRM}/chrM.bam.bai"
                     -ibam "${CHRM}/chrM.bam"  -ibai "${CHRM}/chrM.bam.bai"
                     -g    "${CHRM}/chrM.fa"
                     -o    "${OUT}/cov.sites.bed"
                     -or   "${OUT}/cov.regions.bed"
                     -p    "${OUT}/cov.params.txt"
                     -nt   4)
    set_tests_properties(tier2_run_cov PROPERTIES
        LABELS "tier2;run;covariates" FIXTURES_SETUP cov_output TIMEOUT 900)

    add_check(tier2_cov_smoke   tier2 cov count  "${OUT}/cov.sites.bed" 10)

    # Positions, not bytes. The covariate model fits a GLM, and that fit
    # amplifies small libm differences: between two arm64 Macs the same input
    # gives identical site positions but scores that differ by 1e-4 to 1e-3
    # relative, and fitted params by up to 4e-4. An exact comparison is not
    # portable here even within one platform, so assert what is stable.
    add_check(tier2_cov_sites   tier2 cov agree  "${OUT}/cov.sites.bed"
              "${GOLD_COV}/sites.bed" --min-jaccard 0.98)
    add_check(tier2_cov_params  tier2 cov params "${OUT}/cov.params.txt"
              "${GOLD_COV}/params.txt" --tol 1e-2)
endif()

# --- tier 1: TOML config file -----------------------------------------------
add_test(NAME tier1_config
         COMMAND ${CMAKE_COMMAND}
                 -DPURECLIP=${PURECLIP} -DSYN=${SYN} -DOUT=${OUT}
                 -P "${T}/config.cmake")
set_tests_properties(tier1_config PROPERTIES
    LABELS "tier1;config" FIXTURES_REQUIRED synthetic_data
    SKIP_RETURN_CODE 77 TIMEOUT 600)
