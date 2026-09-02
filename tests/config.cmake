# Exercises the TOML config file: that it is read at all, and that a flag on
# the command line still wins over it. Invoked by the tier1_config test.
#
# Precedence is command line > config file > built-in defaults. It works
# because defaults live in the AppOptions constructor and SeqAn's
# getOptionValue() leaves its destination untouched for a flag the user did
# not pass, so the config can be loaded in between.

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

execute_process(COMMAND python3 ${CMP} bed "${OUT}/cfgonly_sites.bed" "${GOLD}/sites.bed"
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "config-only run does not reproduce the reference output")
endif()

# --- 2. the command line must override the config --------------------------
# bandwidth 100 changes the result; asking for it in the config and then
# passing the default on the command line must land back on the default.
file(WRITE "${OUT}/bw.toml" "bandwidth = 100\n")
execute_process(
    COMMAND ${PURECLIP} -i ${SYN}/sample.bam -bai ${SYN}/sample.bam.bai -g ${SYN}/ref.fa
            -o ${OUT}/ovr.sites.bed -or ${OUT}/ovr.regions.bed -p ${OUT}/ovr.params.txt
            -c "${OUT}/bw.toml" -bw 50
    OUTPUT_QUIET RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "config + overriding flag failed (exit ${rc})")
endif()
execute_process(COMMAND python3 ${CMP} bed "${OUT}/ovr.sites.bed" "${GOLD}/sites.bed"
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "-bw 50 on the command line did not override bandwidth=100 in the config")
endif()
