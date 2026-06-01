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

#include <iostream>
#include <string>

#include "cli/cli.h"

int main(int argc, char const ** argv)
{
    AppOptions options;

    // Save defaults for CLI-override detection
    AppOptions defaults;

    // 1. Parse CLI arguments
    ArgumentParser parser;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // 2. Load TOML config (if -c/--config provided), CLI values take precedence
    std::string config_path;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if ((arg == "-c" || arg == "--config") && i + 1 < argc) {
            config_path = argv[i + 1];
            break;
        }
    }
    if (!config_path.empty()) {
        AppOptions cli_opts = options;  // save what CLI set
        options = defaults;              // reset to defaults
        if (!load_toml_config(config_path, options)) {
            return 1;
        }
        // Restore CLI-set fields (CLI > TOML > defaults)
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
    }

    // 3. Validate required fields
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

    // 4. Auto-derive output prefix from BAM name if not specified
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

#if SEQAN_HAS_ZLIB
    if (options.verbosity > 1) std::cout << "SEQAN_HAS_ZLIB" << std::endl;
#else
    std::cout << "WARNING: zlib not available !" << std::endl;
#endif

    return doIt(options);
}
