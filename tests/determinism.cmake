# Run the synthetic fixture at two thread counts and require identical sites.
# Invoked by the tier1_determinism test; not meant to be run directly.
foreach(NT 1 4)
    execute_process(
        COMMAND ${PURECLIP} -i ${SYN}/sample.bam -bai ${SYN}/sample.bam.bai
                -g ${SYN}/ref.fa -o ${OUT}/det_nt${NT}.sites.bed
                -or ${OUT}/det_nt${NT}.regions.bed -p ${OUT}/det_nt${NT}.params.txt
                -nt ${NT}
        OUTPUT_QUIET RESULT_VARIABLE rc)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "pureclip failed with -nt ${NT} (exit ${rc})")
    endif()
endforeach()

execute_process(
    COMMAND python3 ${CMP} bed ${OUT}/det_nt4.sites.bed ${OUT}/det_nt1.sites.bed
    RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "output differs between -nt 1 and -nt 4")
endif()
