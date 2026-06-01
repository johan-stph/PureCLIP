// ======================================================================
// PureCLIP — command-line parsing and top-level dispatch
// ======================================================================

#ifndef PURECLIP_CLI_H_
#define PURECLIP_CLI_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "version.h"
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/misc/name_store_cache.h>
#include <seqan/arg_parse.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

#include "util.h"
#include "cli/config.h"
#include "call_sites.h"
#include "call_sites_replicates.h"

using namespace seqan;

// ── Argument parsing ────────────────────────────────────────────────────────

inline ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser("pureclip");
    setShortDescription(parser, "Protein-RNA interaction site detection ");
    setVersion(parser, PURE_CLIP_VERSION);
    setDate(parser, PURE_CLIP_BUILD_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <-i \\fIBAM FILE\\fP> <-bai \\fIBAI FILE\\fP> <-g \\fIGENOME FILE\\fP> [-o \\fIOUTPUT BED FILE\\fP] ");
    addDescription(parser, "Protein-RNA interaction site detection using a non-homogeneous HMM.");

    // Config file
    addOption(parser, ArgParseOption("c", "config", "TOML configuration file.  CLI flags override config values.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "config", ".toml");

    // Required files
    addOption(parser, ArgParseOption("i", "in", "Target bam files.", ArgParseArgument::INPUT_FILE, "BAM", true));
    setValidValues(parser, "in", ".bam");

    addOption(parser, ArgParseOption("bai", "bai", "Target bam index files.", ArgParseArgument::INPUT_FILE, "BAI", true));
    setValidValues(parser, "bai", ".bai");

    addOption(parser, ArgParseOption("g", "genome", "Genome reference file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "genome", ".fa .fasta .fa.gz .fasta.gz");

    addOption(parser, ArgParseOption("o", "out", "Output file to write crosslink sites.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "out", ".bed");
    addOption(parser, ArgParseOption("or", "or", "Output file to write binding regions.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "or", ".bed");
    addOption(parser, ArgParseOption("p", "par", "Output file to write learned parameters.", ArgParseArgument::OUTPUT_FILE));

    addSection(parser, "Options");

    addOption(parser, ArgParseOption("ctr", "ctr", "Assign crosslink sites to read start positions. Default: upstream of read starts."));
    addOption(parser, ArgParseOption("st", "st", "Scoring scheme. Default: 0 -> score_UC (log posterior probability ratio).", ArgParseArgument::INTEGER));
    setMinValue(parser, "st", "0");
    setMaxValue(parser, "st", "3");

    addOption(parser, ArgParseOption("iv", "inter", "Genomic intervals for learning, e.g. 'chr1;chr2;chr3'. Default: all contigs from reference.", ArgParseArgument::STRING));
    addOption(parser, ArgParseOption("chr", "chr", "Contigs to apply HMM, e.g. 'chr1;chr2;chr3;'.", ArgParseArgument::STRING));

    addOption(parser, ArgParseOption("bc", "bc", "Binding characteristics: 0=short motifs, 1=larger clusters.", ArgParseArgument::INTEGER));
    setMinValue(parser, "bc", "0");
    setMaxValue(parser, "bc", "1");

    addOption(parser, ArgParseOption("bw", "bdw", "KDE bandwidth for enrichment. Default: 50.", ArgParseArgument::INTEGER));
    setMinValue(parser, "bdw", "1");
    setMaxValue(parser, "bdw", "500");

    addOption(parser, ArgParseOption("bwn", "bdwn", "KDE bandwidth for N estimation. Default: same as bdw.", ArgParseArgument::INTEGER));
    setMinValue(parser, "bdwn", "1");
    setMaxValue(parser, "bdwn", "500");

    addOption(parser, ArgParseOption("kgw", "kgw", "Kernel gap width", ArgParseArgument::INTEGER));
    setMinValue(parser, "kgw", "0");
    setMaxValue(parser, "kgw", "20");
    hideOption(parser, "kgw");

    addOption(parser, ArgParseOption("dm", "dm", "Distance to merge crosslink sites into regions. Default: 8", ArgParseArgument::INTEGER));

    addOption(parser, ArgParseOption("ld", "ld", "Use long double precision (slower, for stability)."));
    addOption(parser, ArgParseOption("ts", "ts", "Log-sum-exp lookup table size. Default: 600000", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("tmv", "tmv", "Log-sum-exp lookup table min value. Default: -2000", ArgParseArgument::DOUBLE));

    addOption(parser, ArgParseOption("ur", "ur", "Select read: 1=R1, 2=R2. eCLIP→R2, iCLIP→R1.", ArgParseArgument::INTEGER));
    setMinValue(parser, "ur", "1");
    setMaxValue(parser, "ur", "2");

    addSection(parser, "Options for incorporating covariates");

    addOption(parser, ArgParseOption("is", "is", "Covariates file: position-wise values (e.g. smoothed input KDEs).", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "is", ".bed");
    addOption(parser, ArgParseOption("ibam", "ibam", "Control experiment BAM file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "ibam", ".bam");
    addOption(parser, ArgParseOption("ibai", "ibai", "Control experiment BAI file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "ibai", ".bai");

    addOption(parser, ArgParseOption("fis", "fis", "Fimo input motif score covariates file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "fis", ".bed");
    addOption(parser, ArgParseOption("nim", "nim", "Max. motif ID to use. Default: 1.", ArgParseArgument::INTEGER));

    addSection(parser, "Advanced user options");

    addOption(parser, ArgParseOption("upe", "upe", "Use pseudo emission probabilities for crosslink state."));
    addOption(parser, ArgParseOption("m", "mibr", "Max iterations within BRENT algorithm.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mibr", "1");
    setMaxValue(parser, "mibr", "1000");
    addOption(parser, ArgParseOption("w", "mibw", "Max iterations within Baum-Welch algorithm.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mibw", "0");
    setMaxValue(parser, "mibw", "500");
    addOption(parser, ArgParseOption("g1kmin", "g1kmin", "Min shape k of non-enriched gamma.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("g1kmax", "g1kmax", "Max shape k of non-enriched gamma.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "g1kmin", "1.5");
    addOption(parser, ArgParseOption("g2kmin", "g2kmin", "Min shape k of enriched gamma.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("g2kmax", "g2kmax", "Max shape k of enriched gamma.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("fk", "fk", "Don't constrain g1.k <= g2.k."));

    addOption(parser, ArgParseOption("mkn", "mkn", "Max. k/N ratio for learning truncation probs. Default: 1.0", ArgParseArgument::DOUBLE));
    setMinValue(parser, "mkn", "0.5");
    setMaxValue(parser, "mkn", "1.5");
    addOption(parser, ArgParseOption("b1p", "b1p", "Initial binomial p for non-crosslink state. Default: 0.01.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("b2p", "b2p", "Initial binomial p for crosslink state. Default: 0.15.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("mtp", "mtp", "Min. transition prob 2→3. Default: 0.0001.", ArgParseArgument::DOUBLE));

    addOption(parser, ArgParseOption("mk", "mkde", "Min KDE for gamma fitting. Default: singleton read start.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("ntp", "ntp", "Min N for learning binomial p. Default: 10", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("ntp2", "ntp2", "Min N for learning transition 2→2/3. Default: 0", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("antp", "antp", "Auto-choose N thresholds from expected read start counts."));

    addOption(parser, ArgParseOption("pa", "pat", "Poly-X stretch length threshold.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("ea1", "epal", "Exclude poly-A from learning."));
    addOption(parser, ArgParseOption("ea2", "epaa", "Exclude poly-A from analysis."));
    addOption(parser, ArgParseOption("et1", "eptl", "Exclude poly-U from learning."));
    addOption(parser, ArgParseOption("et2", "epta", "Exclude poly-U from analysis."));

    addOption(parser, ArgParseOption("mrtf", "mrtf", "Min covariate value for gamma k fitting.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("mtc", "mtc", "Max read starts at one position for learning. Default: 500.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mtc", "50");
    setMaxValue(parser, "mtc", "50000");
    addOption(parser, ArgParseOption("mtc2", "mtc2", "Max read starts stored (truncated above). Default: 65000.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mtc2", "5000");
    setMaxValue(parser, "mtc2", "65000");

    addOption(parser, ArgParseOption("pet", "pet", "Prior enrichment threshold (KDE units). Default: 7", ArgParseArgument::INTEGER));
    setMinValue(parser, "pet", "2");
    setMaxValue(parser, "pet", "50");

    addSection(parser, "General user options");
    addOption(parser, ArgParseOption("nt", "nt", "Number of threads for learning.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("nta", "nta", "Number of threads for applying (default: min(nt, #chromosomes)).", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("oa", "oa", "Output all sites with at least one read start."));
    addOption(parser, ArgParseOption("oe", "oe", "Output all enriched sites."));
    hideOption(parser, "oe");

    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addTextSection(parser, "Parameter settings for proteins with different binding characteristics");
    addText(parser, "By default, parameters are optimized for proteins binding to short defined regions.");
    addListItem(parser, "\\fB0\\fP", "\\fBShort defined\\fP. Default. Equivalent to: \\fB-bdwn 50 -ntp 10 -ntp2 0 -b1p 0.01 -b2p 0.15\\fP.");
    addListItem(parser, "\\fB1\\fP", "\\fBLarger clusters\\fP. Equivalent to: \\fB-bdwn 100 -antp -b2p 0.01 -b2p 0.1\\fP. ");
    addText(parser, "");

    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBpureclip\\fP \\fB-i target.bam\\fP \\fB-bai target.bai\\fP \\fB-g ref.fasta\\fP \\fB-o sites.bed\\fP \\fB-nt 10\\fP  \\fB-iv '1;2;3;'\\fP", "Learn on chr 1-3, 10 threads.");
    addListItem(parser, "\\fBpureclip\\fP \\fB-i r1.bam\\fP \\fB-bai r1.bai\\fP \\fB-i r2.bam\\fP \\fB-bai r2.bai\\fP \\fB-g ref.fa\\fP \\fB-o sites.bed\\fP \\fB-nt 10\\fP", "Two replicates.");
    addListItem(parser, "\\fBpureclip\\fP \\fB-c config.toml\\fP", "All parameters from TOML config.");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    unsigned repNo = getOptionValueCount(parser, "in");
    if (repNo != getOptionValueCount(parser, "bai"))
    {
        std::cout << "ERROR: number of BAI files must equal number of BAM files." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    resize(options.bamFileNames, repNo);
    resize(options.baiFileNames, repNo);
    for (unsigned i = 0; i < repNo; ++i)
    {
        getOptionValue(options.bamFileNames[i], parser, "in", i);
        getOptionValue(options.baiFileNames[i], parser, "bai", i);
    }

    getOptionValue(options.refFileName, parser, "genome");
    getOptionValue(options.outFileName, parser, "out");
    getOptionValue(options.outRegionsFileName, parser, "or");
    getOptionValue(options.parFileName, parser, "par");
    getOptionValue(options.rpkmFileName, parser, "is");
    getOptionValue(options.inputBamFileName, parser, "ibam");
    getOptionValue(options.inputBaiFileName, parser, "ibai");
    if ((options.rpkmFileName != "" && options.inputBamFileName != "") ||
            (options.rpkmFileName != "" && options.inputBaiFileName != "") ||
            (options.inputBamFileName != "" && options.inputBaiFileName == "") ||
            (options.inputBamFileName == "" && options.inputBaiFileName != "") )
    {
        std::cout << "ERROR: If using background covariates, either -is or -ibam+ibai must be given." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    if (options.rpkmFileName != "" || options.inputBamFileName != "")
        options.useCov_RPKM = true;
    getOptionValue(options.fimoFileName, parser, "fis");
    if (options.fimoFileName != "")
        options.useFimoScore = true;

    if (isSet(parser, "ctr"))
        options.crosslinkAtTruncSite = true;
    getOptionValue(options.score_type, parser, "st");
    getOptionValue(options.intervals_str, parser, "inter");
    if (isSet(parser, "upe"))
        options.use_pseudoEProb = true;
    getOptionValue(options.maxIter_brent, parser, "mibr");
    getOptionValue(options.maxIter_bw, parser, "mibw");
    getOptionValue(options.g1_kMin, parser, "g1kmin");
    getOptionValue(options.g1_kMax, parser, "g1kmax");
    getOptionValue(options.g2_kMin, parser, "g2kmin");
    getOptionValue(options.g2_kMax, parser, "g2kmax");
    if (isSet(parser, "fk"))
        options.g1_k_le_g2_k = false;

    unsigned bc = 0;
    getOptionValue(bc, parser, "bc");
    if (bc == 1)
    {
        options.bandwidthN = 100;
        options.get_nThreshold = true;
        options.p1 = 0.01;
        options.p2 = 0.1;
    }

    getOptionValue(options.bandwidth, parser, "bdw");
    getOptionValue(options.bandwidthN, parser, "bdwn");
    if (options.bandwidthN == 0)
        options.bandwidthN = options.bandwidth;
    getOptionValue(options.nKernelGap, parser, "kgw");

    getOptionValue(options.useKdeThreshold, parser, "mkde");
    getOptionValue(options.nThresholdForP, parser, "ntp");
    getOptionValue(options.nThresholdForTransP, parser, "ntp2");
    if (isSet(parser, "antp"))
        options.get_nThreshold = true;
    getOptionValue(options.minTransProbCS, parser, "mtp");
    getOptionValue(options.maxkNratio, parser, "mkn");
    getOptionValue(options.p1, parser, "b1p");
    getOptionValue(options.p2, parser, "b2p");

    getOptionValue(options.distMerge, parser, "dm");
    if (isSet(parser, "ld"))
        options.useHighPrecision = true;
    getOptionValue(options.lookupTable_size, parser, "ts");
    getOptionValue(options.lookupTable_minValue, parser, "tmv");
    getOptionValue(options.selectRead, parser, "ur");

    getOptionValue(options.polyAThreshold, parser, "pat");
    if (isSet(parser, "epal")) options.excludePolyAFromLearning = true;
    if (isSet(parser, "epaa")) options.excludePolyA = true;
    if (isSet(parser, "eptl")) options.excludePolyTFromLearning = true;
    if (isSet(parser, "epta")) options.excludePolyT = true;

    getOptionValue(options.minRPKMtoFit, parser, "mrtf");
    if (isSet(parser, "mrtf")) options.mrtf_kdeSglt = false;

    getOptionValue(options.maxTruncCount, parser, "mtc");
    getOptionValue(options.maxTruncCount2, parser, "mtc2");
    getOptionValue(options.nInputMotifs, parser, "nim");
    getOptionValue(options.prior_enrichmentThreshold, parser, "pet");
    getOptionValue(options.numThreads, parser, "nt");
    getOptionValue(options.numThreadsA, parser, "nta");

    if (isSet(parser, "oa")) options.outputAll = true;
    if (isSet(parser, "quiet")) options.verbosity = 0;
    if (isSet(parser, "verbose")) options.verbosity = 2;
    if (isSet(parser, "very-verbose")) options.verbosity = 3;

    getOptionValue(options.applyChr_str, parser, "chr");

    return ArgumentParser::PARSE_OK;
}


// ── Top-level dispatch ──────────────────────────────────────────────────────

template <typename TOptions>
inline bool doIt(TOptions &options)
{
    unsigned repNo = length(options.baiFileNames);

    LogSumExp_lookupTable lookUp(options.lookupTable_size, options.lookupTable_minValue);
    options.lookUp = lookUp;

    if (options.useCov_RPKM)
    {
        if (options.useFimoScore)
        {
            ModelParams<GAMMA_REG, ZTBIN_REG> modelParams;
            modelParams.gamma1.tp = options.useKdeThreshold;
            modelParams.gamma2.tp = options.useKdeThreshold;
            modelParams.bin1.b0 = log(options.p1/(1.0 - options.p1));
            modelParams.bin2.b0 = log(options.p2/(1.0 - options.p2));
            resize(modelParams.bin1.regCoeffs, options.nInputMotifs, 0.0, Exact());
            resize(modelParams.bin2.regCoeffs, options.nInputMotifs, 0.0, Exact());
            String<ModelParams<GAMMA_REG, ZTBIN_REG> > modelParams_reps;
            resize(modelParams_reps, repNo, modelParams);
            return doIt(modelParams_reps, options);
        }
        else
        {
            ModelParams<GAMMA_REG, ZTBIN> modelParams;
            modelParams.gamma1.tp = options.useKdeThreshold;
            modelParams.gamma2.tp = options.useKdeThreshold;
            modelParams.bin1.p = options.p1;
            modelParams.bin2.p = options.p2;
            String<ModelParams<GAMMA_REG, ZTBIN> > modelParams_reps;
            resize(modelParams_reps, repNo, modelParams);
            return doIt(modelParams_reps, options);
        }
    }
    else
    {
        options.g1_kMax = 1.0;
        if (options.verbosity > 1) std::cout << "Note: set max. value of g1.k to 1.0." << std::endl;

        if (options.useFimoScore)
        {
            ModelParams<GAMMA, ZTBIN_REG> modelParams;
            modelParams.gamma1.tp = options.useKdeThreshold;
            modelParams.gamma2.tp = options.useKdeThreshold;
            modelParams.bin1.b0 = log(options.p1/(1.0 - options.p1));
            modelParams.bin2.b0 = log(options.p2/(1.0 - options.p2));
            resize(modelParams.bin1.regCoeffs, options.nInputMotifs, 0.0, Exact());
            resize(modelParams.bin2.regCoeffs, options.nInputMotifs, 0.0, Exact());
            String<ModelParams<GAMMA, ZTBIN_REG> > modelParams_reps;
            resize(modelParams_reps, repNo, modelParams);
            return doIt(modelParams_reps, options);
        }
        else
        {
            ModelParams<GAMMA, ZTBIN> modelParams;
            modelParams.gamma1.tp = options.useKdeThreshold;
            modelParams.gamma2.tp = options.useKdeThreshold;
            modelParams.bin1.p = options.p1;
            modelParams.bin2.p = options.p2;
            String<ModelParams<GAMMA, ZTBIN> > modelParams_reps;
            resize(modelParams_reps, repNo, modelParams);
            return doIt(modelParams_reps, options);
        }
    }
}

#endif // PURECLIP_CLI_H_
