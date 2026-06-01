// ======================================================================
// PureCLIP: capturing target-specific protein-RNA interaction footprints
// ======================================================================
// Copyright (C) 2017  Sabrina Krakau, Max Planck Institute for Molecular 
// Genetics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Sabrina Krakau <krakau@molgen.mpg.de>
// ==========================================================================

#define HMM_PROFILE

	
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
#include "config.h"
#include "call_sites.h"
#include "call_sites_replicates.h"


using namespace seqan;

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("pureclip");
    // Set short description, version, and date.
    setShortDescription(parser, "Protein-RNA interaction site detection ");
    setVersion(parser, PURE_CLIP_VERSION);
    setDate(parser, PURE_CLIP_BUILD_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <-i \\fIBAM FILE\\fP> <-bai \\fIBAI FILE\\fP> <-g \\fIGENOME FILE\\fP> [-o \\fIOUTPUT BED FILE\\fP] ");
    addDescription(parser, "Protein-RNA interaction site detection using a non-homogeneous HMM.");

    // Config file (loaded before CLI flags — CLI overrides TOML values)
    addOption(parser, ArgParseOption("c", "config", "TOML configuration file.  CLI flags override config values.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "config", ".toml");

    // rep1 [rep2]
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
    //setRequired(parser, "par", true);
    

    addSection(parser, "Options");

    addOption(parser, ArgParseOption("ctr", "ctr", "Assign crosslink sites to read start positions. Note: depends on RT enzyme, buffer conditions and likely on protein. Default: assign crosslink sites to positions upstream of read starts."));
    addOption(parser, ArgParseOption("st", "st", "Scoring scheme. Default: 0 -> score_UC (log posterior probability ratio of most likely and second most likely state).", ArgParseArgument::INTEGER));
    setMinValue(parser, "st", "0");
    setMaxValue(parser, "st", "3");

    addOption(parser, ArgParseOption("iv", "inter", "Genomic chromosomes to learn HMM parameters, e.g. 'chr1;chr2;chr3'. Contigs have to be in the same order as in BAM file. Useful to reduce runtime and memory consumption. Default: all contigs from reference file are used (useful when applying to transcript-wise alignments or poor data).", ArgParseArgument::STRING));
    addOption(parser, ArgParseOption("chr", "chr", "Contigs to apply HMM, e.g. 'chr1;chr2;chr3;'. Contigs have to be in the same order as in BAM file.", ArgParseArgument::STRING));

    addOption(parser, ArgParseOption("bc", "bc", "Flag to set parameters according to binding characteristics of protein: see description in section below.", ArgParseArgument::INTEGER));
    //setValidValues(parser, "bc", "0 1");
    setMinValue(parser, "bc", "0");
    setMaxValue(parser, "bc", "1");

    addOption(parser, ArgParseOption("bw", "bdw", "Bandwidth for kernel density estimation used to access enrichment. NOTE: Increasing the bandwidth increases runtime and memory consumption. Default: 50.", ArgParseArgument::INTEGER));
    setMinValue(parser, "bdw", "1");
    setMaxValue(parser, "bdw", "500");

    addOption(parser, ArgParseOption("bwn", "bdwn", "Bandwidth for kernel density estimation used to estimate n for binomial distributions. For proteins that are rather sliding along the RNA or showing long crosslink clusters this should be increased, e.g. to 100 (should be <= 4*bdw). Default: same as bdw.", ArgParseArgument::INTEGER));
    setMinValue(parser, "bdwn", "1");
    setMaxValue(parser, "bdwn", "500"); 
    addOption(parser, ArgParseOption("kgw", "kgw", "Kernel gap width", ArgParseArgument::INTEGER));
    setMinValue(parser, "kgw", "0");
    setMaxValue(parser, "kgw", "20");
    hideOption(parser, "kgw");

    addOption(parser, ArgParseOption("dm", "dm", "Distance used to merge individual crosslink sites to binding regions. Default: 8", ArgParseArgument::INTEGER));

    addOption(parser, ArgParseOption("ld", "ld", "Use higher precision to store emission probabilities, state poster posterior probabilities etc. (i.e. long double). Should not be necessary anymore, due to computations in log-space. Note: increases memory consumption. Default: double."));
    addOption(parser, ArgParseOption("ts", "ts", "Size of look-up table for log-sum-exp values. Default: 600000", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("tmv", "tmv", "Minimum value in look-up table for log-sum-exp values. Default: -2000", ArgParseArgument::DOUBLE));

    addOption(parser, ArgParseOption("ur", "ur", "Flag to define which read should be selected for the analysis: 1->R1, 2->R2. Note: PureCLIP uses read starts corresponding to 3' cDNA ends. Thus if providing paired-end data, only the corresponding read should be selected (e.g. eCLIP->R2, iCLIP->R1). If applicable, used for input BAM file as well. Default: uses read starts of all provided reads assuming single-end or pre-filtered data.", ArgParseArgument::INTEGER));
    setMinValue(parser, "ur", "1");
    setMaxValue(parser, "ur", "2");

    addSection(parser, "Options for incorporating covariates");

    addOption(parser, ArgParseOption("is", "is", "Covariates file: position-wise values, e.g. smoothed reads start counts (KDEs) from input data. ", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "is", ".bed");
    addOption(parser, ArgParseOption("ibam", "ibam", "File containing mapped reads from control experiment, e.g. eCLIP input.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "ibam", ".bam");
    addOption(parser, ArgParseOption("ibai", "ibai", "File containing BAM index corresponding to mapped reads from control experiment", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "ibai", ".bai");

    addOption(parser, ArgParseOption("fis", "fis", "Fimo input motif score covariates file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "fis", ".bed");
    addOption(parser, ArgParseOption("nim", "nim", "Max. motif ID to use. Default: Only covariates with motif ID 1 are used.", ArgParseArgument::INTEGER));

    addSection(parser, "Advanced user options");

    addOption(parser, ArgParseOption("upe", "upe", "Use (n dependent) pseudo emission probabilities for crosslink state."));

    addOption(parser, ArgParseOption("m", "mibr", "Maximum number of iterations within BRENT algorithm.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mibr", "1");
    setMaxValue(parser, "mibr", "1000");
    addOption(parser, ArgParseOption("w", "mibw", "Maximum number of iterations within Baum-Welch algorithm.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mibw", "0");
    setMaxValue(parser, "mibw", "500");
    addOption(parser, ArgParseOption("g1kmin", "g1kmin", "Minimum shape k of 'non-enriched' gamma distribution (g1.k).", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("g1kmax", "g1kmax", "Maximum shape k of 'non-enriched' gamma distribution (g1.k).", ArgParseArgument::DOUBLE));
    setMinValue(parser, "g1kmin", "1.5");   
    addOption(parser, ArgParseOption("g2kmin", "g2kmin", "Minimum shape k of 'enriched' gamma distribution (g2.k).", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("g2kmax", "g2kmax", "Maximum shape k of 'enriched' gamma distribution (g2.k).", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("fk", "fk", "When incorporating input signal, do not constrain 'non-enriched' shape parameter k <= 'enriched' gamma parameter k."));

    addOption(parser, ArgParseOption("mkn", "mkn", "Max. k/N ratio (read start sites/N) used to learn truncation probabilities for 'non-crosslink' and 'crosslink' emission probabilities (high ratios might originate from mapping artifacts that can disturb parameter learning). Default: 1.0", ArgParseArgument::DOUBLE));
    setMinValue(parser, "mkn", "0.5");
    setMaxValue(parser, "mkn", "1.5");
    addOption(parser, ArgParseOption("b1p", "b1p", "Initial value for binomial probability parameter of 'non-crosslink' state. Default: 0.01.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("b2p", "b2p", "Initial value for binomial probability parameter of 'crosslink' state. Default: 0.15.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("mtp", "mtp", "Min. transition probability from state '2' to '3' (helpful for poor data, where no clear distinction between 'enriched' and 'non-enriched' is possible). Default: 0.0001.", ArgParseArgument::DOUBLE));

    addOption(parser, ArgParseOption("mk", "mkde", "Minimum KDE value used for fitting left-truncated gamma distributions. Default: corresponding to singleton read start.", ArgParseArgument::DOUBLE));

    addOption(parser, ArgParseOption("ntp", "ntp", "Only sites with n >= ntp are used to learn binomial probability parameters (bin1.p, bin2.p). Default: 10", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("ntp2", "ntp2", "Only sites with n >= ntp2 are used to learn probability of transition from state '2' to '2' or '3'. Useful for data with low truncation rate at crosslink sites or in general high fraction of non-coinciding read starts. Default: 0", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("antp", "antp", "Automatically choose n threshold (-ntp, -ntp2) to estimate parameters linked to crosslink states based on expected read start count at crosslink sites."));

    addOption(parser, ArgParseOption("pa", "pat", "Length threshold for internal poly-X stretches to get excluded.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("ea1", "epal", "Exclude intervals containing poly-A stretches from learning."));
    addOption(parser, ArgParseOption("ea2", "epaa", "Exclude intervals containing poly-A stretches from analysis."));
    addOption(parser, ArgParseOption("et1", "eptl", "Exclude intervals containing poly-U stretches from learning."));
    addOption(parser, ArgParseOption("et2", "epta", "Exclude intervals containing poly-U stretches from analysis."));
 
    addOption(parser, ArgParseOption("mrtf", "mrtf", "Fit gamma shape k only for positions with min. covariate value.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("mtc", "mtc", "Maximum number of read starts at one position used for learning. For sites with counts above threshold the whole covered regions will be ignored for learning! Default: 500.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mtc", "50");
    setMaxValue(parser, "mtc", "50000");
    addOption(parser, ArgParseOption("mtc2", "mtc2", "Maximum number of read starts at one position stored. For sites with counts above threshold the count will be truncated. Influences k and n. Default: 65000.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mtc2", "5000");
    setMaxValue(parser, "mtc2", "65000"); 

    addOption(parser, ArgParseOption("pet", "pet", "Prior enrichment threshold: a KDE threshold corresponding to 7 read start counts at one position will be used for initial classification of 'non-enriched' and 'enriched' site. Default: 7", ArgParseArgument::INTEGER));
    setMinValue(parser, "pet", "2");
    setMaxValue(parser, "pet", "50");

    addSection(parser, "General user options");
    addOption(parser, ArgParseOption("nt", "nt", "Number of threads used for learning.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("nta", "nta", "Number of threads used for applying learned parameters. Increases memory usage, if greater than number of chromosomes used for learning, since HMM will be build for multiple chromosomes in parallel. Default: min(nt, no. of chromosomes/transcripts used for learning).", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("oa", "oa", "Outputs all sites with at least one read start in extended output format."));
    addOption(parser, ArgParseOption("oe", "oe", "Outputs additionally all sites that are 'enriched' and contain at least one read start."));
    hideOption(parser, "oe");

    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));


    //addOption(parser, ArgParseOption("enk", "enk", "Estimate binomial N from KDEs values."));
    //addOption(parser, ArgParseOption("ulr", "ulr", "Use log RPKM values as input."));
    //addOption(parser, ArgParseOption("dis", "dis", "Discard intervals containing only one read start."));


    addTextSection(parser, "Parameter settings for proteins with different binding characteristics");
    addText(parser, "By default, the parameters are set to values optimized for proteins binding to short defined binding regions, e.g. proteins binding to short specific motifs such as PUM2 and RBFOX2. "
                    "With the \\fB-bc\\fP option this behaviour can be changed:");
    addListItem(parser, "\\fB0\\fP", "\\fBShort defined\\fP. Default. Equivalent to: \\fB-bdwn 50 -ntp 10 -ntp2 0 -b1p 0.01 -b2p 0.15\\fP.");
    addListItem(parser, "\\fB1\\fP", "\\fBLarger clusters\\fP. Proteins causing larger crosslink clusters with relatively lower read start counts, e.g. proteins binding to low complexity motifs. Equivalent to: \\fB-bdwn 100 -antp -b2p 0.01 -b2p 0.1\\fP. ");
    addText(parser, "");
    addText(parser, "In case of different binding characteristics adjust parameters \\fB-bdw\\fP, \\fB-bdwn\\fP, \\fB-b1p\\fP, \\fB-b2p\\fP, \\fB-antp\\fP or see \\fIhttp://pureclip.readthedocs.io/en/latest/PureCLIPTutorial/userOptions.html\\fP for more information. ");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBpureclip\\fP \\fB-i target.bam\\fP \\fB-bai target.bai\\fP \\fB-g ref.fasta\\fP \\fB-o called_crosslinksites.bed\\fP \\fB-nt 10\\fP  \\fB-iv '1;2;3;'\\fP", "Learn HMM parameters on chromosomes 1-3, use 10 threads for learning and otherwise default parameters.");
    addListItem(parser, "\\fBpureclip\\fP \\fB-i target.rep1.bam\\fP \\fB-bai target.rep1.bai\\fP \\fB-i target.rep2.bam\\fP \\fB-bai target.rep2.bai\\fP \\fB-g ref.fasta\\fP \\fB-o called_crosslinksites.bed\\fP \\fB-nt 10\\fP", "Include individual replicates (currently only supported for two), while learning parameters on whole datasets.");
    addListItem(parser, "\\fBpureclip\\fP \\fB-i target.bam\\fP \\fB-bai target.bai\\fP \\fB-g ref.fasta\\fP \\fB-o called_crosslinksites.bed\\fP \\fB-nt 10\\fP  \\fB-iv '1;2;3;'\\fP \\fB-bc 1\\fP ", "Use parameter settings for proteins causing larger crosslink clusters");
    addListItem(parser, "\\fBpureclip\\fP \\fB-i target.bam\\fP \\fB-bai target.bai\\fP \\fB-g ref.fasta\\fP \\fB-o called_crosslinksites.bed\\fP \\fB-nt 10\\fP  \\fB-iv '1;2;3;'\\fP \\fB-bc 1 -b2p 0.03\\fP ", "Use parameter settings for proteins causing larger crosslink clusters and decrease initial probability parameter for 'crosslink' state for data with high fraction of non-coinciding read starts.");
    addListItem(parser, "\\fBpureclip\\fP \\fB-i target.bam\\fP \\fB-bai target.bai\\fP \\fB-g ref.fasta\\fP \\fB-o called_crosslinksites.bed\\fP \\fB-nt 10\\fP  \\fB-iv '1;2;3;'\\fP \\fB-bdw 25\\fP ", "Use decreased bandwidth of 25 bp to access enrichment.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    unsigned repNo = getOptionValueCount(parser, "in");
    if (repNo != getOptionValueCount(parser, "bai"))
    {
        std::cout << "ERROR: currently only support for <= 2 replicates!" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    if (repNo != getOptionValueCount(parser, "bai"))
    {
        std::cout << "ERROR: number of BAI files must be the same as of BAM files!" << std::endl;
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
        std::cout << "ERROR: If using background signal as covariates, either -is or -ibam and -ibai must be given!" << std::endl;
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
    }   // can be overwritten afterwards ....

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
    if (isSet(parser, "epal"))
        options.excludePolyAFromLearning = true;
    if (isSet(parser, "epaa"))
        options.excludePolyA = true;
    if (isSet(parser, "eptl"))
        options.excludePolyTFromLearning = true;
    if (isSet(parser, "epta"))
        options.excludePolyT = true;
 
    getOptionValue(options.minRPKMtoFit, parser, "mrtf");
    if (isSet(parser, "mrtf"))
        options.mrtf_kdeSglt = false;

    getOptionValue(options.maxTruncCount, parser, "mtc");
    getOptionValue(options.maxTruncCount2, parser, "mtc2");

    getOptionValue(options.nInputMotifs, parser, "nim");

    getOptionValue(options.prior_enrichmentThreshold, parser, "pet");

    getOptionValue(options.numThreads, parser, "nt");
    getOptionValue(options.numThreadsA, parser, "nta");

    if (isSet(parser, "oa"))
        options.outputAll = true;
 
    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;


    //if (isSet(parser, "dis"))
    //    options.discardSingletonIntervals = true;
    getOptionValue(options.applyChr_str, parser, "chr");

    return ArgumentParser::PARSE_OK;
}


template <typename TOptions>
bool doIt(TOptions &options)
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
            // for each replicate
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
            // for each replicate
            String<ModelParams<GAMMA_REG, ZTBIN> > modelParams_reps;
            resize(modelParams_reps, repNo, modelParams);

            return doIt(modelParams_reps, options);  
        }
    }
    else
    {
        options.g1_kMax = 1.0;
        if (options.verbosity > 1) std::cout << "Note: set max. value of g1.k (shape parameter of 'non-enriched' gamma distribution) to 1.0." << std::endl;

        if (options.useFimoScore)
        {
            ModelParams<GAMMA, ZTBIN_REG> modelParams;
            modelParams.gamma1.tp = options.useKdeThreshold;
            modelParams.gamma2.tp = options.useKdeThreshold;             
            modelParams.bin1.b0 = log(options.p1/(1.0 - options.p1)); 
            modelParams.bin2.b0 = log(options.p2/(1.0 - options.p2));
            resize(modelParams.bin1.regCoeffs, options.nInputMotifs, 0.0, Exact());
            resize(modelParams.bin2.regCoeffs, options.nInputMotifs, 0.0, Exact());
            // for each replicate
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
            // for each replicate
            String<ModelParams<GAMMA, ZTBIN> > modelParams_reps;
            resize(modelParams_reps, repNo, modelParams);
            return doIt(modelParams_reps, options); 
        }
    }
}


int main(int argc, char const ** argv)
{
    AppOptions options;

    // Save defaults so we can detect which fields CLI explicitly set
    AppOptions defaults;

    // 1. Parse CLI first (flags not required — TOML or auto-derive handles missing)
    ArgumentParser parser;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // 2. Load TOML config, but only fill in fields that CLI didn't already set
    //    (detected by comparing against the default AppOptions values)
    std::string config_path;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if ((arg == "-c" || arg == "--config") && i + 1 < argc) {
            config_path = argv[i + 1];
            break;
        }
    }
    if (!config_path.empty()) {
        // Merge: load TOML, then restore any fields that CLI explicitly set
        AppOptions cli_opts = options;  // what CLI set
        options = defaults;              // reset to defaults
        if (!load_toml_config(config_path, options)) {
            return 1;
        }
        // Restore CLI-set fields (CLI overrides TOML)
        if (cli_opts.bandwidth != defaults.bandwidth) options.bandwidth = cli_opts.bandwidth;
        if (cli_opts.bandwidthN != defaults.bandwidthN) options.bandwidthN = cli_opts.bandwidthN;
        if (cli_opts.numThreads != defaults.numThreads) options.numThreads = cli_opts.numThreads;
        if (cli_opts.numThreadsA != defaults.numThreadsA) options.numThreadsA = cli_opts.numThreadsA;
        if (cli_opts.verbosity != defaults.verbosity) options.verbosity = cli_opts.verbosity;
        if (cli_opts.score_type != defaults.score_type) options.score_type = cli_opts.score_type;
        if (cli_opts.distMerge != defaults.distMerge) options.distMerge = cli_opts.distMerge;
        if (cli_opts.maxIter_bw != defaults.maxIter_bw) options.maxIter_bw = cli_opts.maxIter_bw;
        if (cli_opts.maxIter_simplex != defaults.maxIter_simplex) options.maxIter_simplex = cli_opts.maxIter_simplex;
        if (cli_opts.g1_kMin != defaults.g1_kMin) options.g1_kMin = cli_opts.g1_kMin;
        if (cli_opts.g1_kMax != defaults.g1_kMax) options.g1_kMax = cli_opts.g1_kMax;
        if (cli_opts.g2_kMin != defaults.g2_kMin) options.g2_kMin = cli_opts.g2_kMin;
        if (cli_opts.g2_kMax != defaults.g2_kMax) options.g2_kMax = cli_opts.g2_kMax;
        if (cli_opts.p1 != defaults.p1) options.p1 = cli_opts.p1;
        if (cli_opts.p2 != defaults.p2) options.p2 = cli_opts.p2;
        if (cli_opts.minTransProbCS != defaults.minTransProbCS) options.minTransProbCS = cli_opts.minTransProbCS;
        if (cli_opts.prior_enrichmentThreshold != defaults.prior_enrichmentThreshold) options.prior_enrichmentThreshold = cli_opts.prior_enrichmentThreshold;
        if (cli_opts.useHighPrecision != defaults.useHighPrecision) options.useHighPrecision = cli_opts.useHighPrecision;
        if (cli_opts.nThresholdForP != defaults.nThresholdForP) options.nThresholdForP = cli_opts.nThresholdForP;
        if (cli_opts.nThresholdForTransP != defaults.nThresholdForTransP) options.nThresholdForTransP = cli_opts.nThresholdForTransP;
        if (cli_opts.get_nThreshold != defaults.get_nThreshold) options.get_nThreshold = cli_opts.get_nThreshold;
        if (cli_opts.maxkNratio != defaults.maxkNratio) options.maxkNratio = cli_opts.maxkNratio;
        if (cli_opts.maxTruncCount != defaults.maxTruncCount) options.maxTruncCount = cli_opts.maxTruncCount;
        if (cli_opts.maxTruncCount2 != defaults.maxTruncCount2) options.maxTruncCount2 = cli_opts.maxTruncCount2;
        if (!empty(cli_opts.bamFileNames)) options.bamFileNames = cli_opts.bamFileNames;
        if (!empty(cli_opts.baiFileNames)) options.baiFileNames = cli_opts.baiFileNames;
        if (cli_opts.refFileName != "") options.refFileName = cli_opts.refFileName;
        if (cli_opts.outFileName != "") options.outFileName = cli_opts.outFileName;
        if (cli_opts.outRegionsFileName != "") options.outRegionsFileName = cli_opts.outRegionsFileName;
        if (cli_opts.parFileName != "") options.parFileName = cli_opts.parFileName;
        if (cli_opts.intervals_str != "") options.intervals_str = cli_opts.intervals_str;
        if (cli_opts.applyChr_str != "") options.applyChr_str = cli_opts.applyChr_str;
        if (cli_opts.rpkmFileName != "") options.rpkmFileName = cli_opts.rpkmFileName;
        if (cli_opts.inputBamFileName != "") options.inputBamFileName = cli_opts.inputBamFileName;
        if (cli_opts.inputBaiFileName != "") options.inputBaiFileName = cli_opts.inputBaiFileName;
        if (cli_opts.fimoFileName != "") options.fimoFileName = cli_opts.fimoFileName;
        if (cli_opts.useCov_RPKM != defaults.useCov_RPKM) options.useCov_RPKM = cli_opts.useCov_RPKM;
        if (cli_opts.useFimoScore != defaults.useFimoScore) options.useFimoScore = cli_opts.useFimoScore;
        if (cli_opts.nInputMotifs != defaults.nInputMotifs) options.nInputMotifs = cli_opts.nInputMotifs;
        if (cli_opts.outputAll != defaults.outputAll) options.outputAll = cli_opts.outputAll;
        if (cli_opts.crosslinkAtTruncSite != defaults.crosslinkAtTruncSite) options.crosslinkAtTruncSite = cli_opts.crosslinkAtTruncSite;
        if (cli_opts.selectRead != defaults.selectRead) options.selectRead = cli_opts.selectRead;
        if (cli_opts.use_pseudoEProb != defaults.use_pseudoEProb) options.use_pseudoEProb = cli_opts.use_pseudoEProb;
        if (cli_opts.g1_k_le_g2_k != defaults.g1_k_le_g2_k) options.g1_k_le_g2_k = cli_opts.g1_k_le_g2_k;
        if (cli_opts.maxIter_brent != defaults.maxIter_brent) options.maxIter_brent = cli_opts.maxIter_brent;
        // Re-apply binding characteristics if CLI set it
        if (cli_opts.bandwidthN == 100 && defaults.bandwidthN == 0 && cli_opts.get_nThreshold == true)
            ; // bc=1 was set by CLI, already captured above
    }

    // 3. Validate required fields (may come from TOML or CLI)
    if (empty(options.bamFileNames)) {
        std::cerr << "ERROR: no BAM files specified.  Use -i or set [bam] in config." << std::endl;
        return 1;
    }
    if (empty(options.baiFileNames)) {
        std::cerr << "ERROR: no BAI files specified.  Use -bai or set [bai] in config." << std::endl;
        return 1;
    }
    if (options.refFileName == "") {
        std::cerr << "ERROR: no genome file specified.  Use -g or set [genome] in config." << std::endl;
        return 1;
    }

    // 4. Auto-derive output name if not set
    if (options.outFileName == "") {
        std::string bam = toCString(options.bamFileNames[0]);
        size_t slash = bam.find_last_of("/\\");
        size_t dot = bam.find_last_of('.');
        std::string base = (slash != std::string::npos) ? bam.substr(slash + 1, dot - slash - 1)
                                                         : bam.substr(0, dot);
        options.outFileName = (base + "_sites.bed").c_str();
        if (options.outRegionsFileName == "")
            options.outRegionsFileName = (base + "_regions.bed").c_str();
        if (options.parFileName == "")
            options.parFileName = (base + "_params.txt").c_str();
        if (options.verbosity >= 1)
            std::cout << "Auto-derived output prefix from BAM: " << base << std::endl;
    }

    std::cout << "Protein-RNA crosslink site detection \n"
              << "===============\n\n";

    // Print the command line arguments back to the user.
    /*if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << "\n\n";
    }*/
    /////////////////////////////////////////
#if SEQAN_HAS_ZLIB
    if (options.verbosity > 1) std::cout << "SEQAN_HAS_ZLIB" << std::endl;
#else
    std::cout << "WARNING: zlib not available !" << std::endl;
#endif

    return doIt(options);
}






