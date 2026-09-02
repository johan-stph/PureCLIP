# Exercises the TOML config file: that it is read at all, and that a flag on
# the command line still wins over it.
#
# Deliberately compares against runs made in this same invocation rather than
# against a golden file. The comparison is "does a config-driven run match the
# equivalent command line", which is self-contained and therefore works on a
# platform with no recorded reference — where a golden comparison would skip,
# and this check would silently stop testing anything.

# base -> writes <base>_sites.bed, <base>_regions.bed, <base>_params.txt
function(run_pureclip base)
    execute_process(COMMAND ${PURECLIP}
                    -i ${SYN}/sample.bam -bai ${SYN}/sample.bam.bai -g ${SYN}/ref.fa
                    -o ${base}_sites.bed -or ${base}_regions.bed -p ${base}_params.txt
                    ${ARGN}
                    OUTPUT_QUIET RESULT_VARIABLE rc)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "pureclip failed (exit ${rc}) with extra args: ${ARGN}")
    endif()
endfunction()

function(expect_same a b msg)
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${a}" "${b}"
                    RESULT_VARIABLE rc)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "${msg}")
    endif()
endfunction()

# Reference: a plain, fully explicit command line.
run_pureclip("${OUT}/cfgref")

# --- 1. a config supplying everything; no flags on the command line ---------
file(WRITE "${OUT}/only.toml"
"bam = [\"${SYN}/sample.bam\"]\n"
"bai = [\"${SYN}/sample.bam.bai\"]\n"
"genome = \"${SYN}/ref.fa\"\n"
"output_prefix = \"${OUT}/cfgonly\"\n")

execute_process(COMMAND ${PURECLIP} -c "${OUT}/only.toml"
                OUTPUT_QUIET RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "config-only invocation failed (exit ${rc})")
endif()
foreach(f cfgonly_sites.bed cfgonly_regions.bed cfgonly_params.txt)
    if(NOT EXISTS "${OUT}/${f}")
        message(FATAL_ERROR "output_prefix did not produce ${f}")
    endif()
endforeach()
expect_same("${OUT}/cfgonly_sites.bed" "${OUT}/cfgref_sites.bed"
            "a config-only run does not match the equivalent command line")

# --- 2. the config must actually take effect -------------------------------
# bandwidth 100 changes the result, so a config asking for it must differ from
# the reference. Without this the override check below would pass vacuously if
# the config were never read at all.
file(WRITE "${OUT}/bw.toml" "bandwidth = 100\n")
run_pureclip("${OUT}/cfgbw" -c "${OUT}/bw.toml")
execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
                "${OUT}/cfgbw_sites.bed" "${OUT}/cfgref_sites.bed"
                RESULT_VARIABLE rc)
if(rc EQUAL 0)
    message(FATAL_ERROR "bandwidth=100 in the config changed nothing — config not applied")
endif()

# --- 3. and a flag on the command line must override it --------------------
run_pureclip("${OUT}/cfgovr" -c "${OUT}/bw.toml" -bw 50)
expect_same("${OUT}/cfgovr_sites.bed" "${OUT}/cfgref_sites.bed"
            "-bw 50 on the command line did not override bandwidth=100 in the config")
