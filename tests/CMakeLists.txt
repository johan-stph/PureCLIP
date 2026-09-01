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
