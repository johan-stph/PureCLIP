// ======================================================================
// PureCLIP — BED output: crosslink sites and binding regions
// ======================================================================

#ifndef PURECLIP_IO_BED_H_
#define PURECLIP_IO_BED_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <seqan/bed_io.h>
#include <seqan/store.h>

#include "types.h"
#include "util.h"

using namespace seqan;

// ── Scoring utility ─────────────────────────────────────────────────────────

inline double getCrosslinkSiteScore(double postProb0, double postProb1,
                                     double postProb2, double postProb3,
                                     unsigned score_type)
{
    if (score_type == 1)        // log(3/2) "crosslink focussed"
    {
        return log(postProb3 / std::max(postProb2, pureclip::float_min()));
    }
    else if (score_type == 2)   // log(3/1) "enrichment focussed"
    {
        return log(postProb3 / std::max(postProb1, pureclip::float_min()));
    }
    else if (score_type == 3)   // balanced: log(enriched/non-enriched) + log(crosslinked/non-crosslinked)
    {
        return (log((postProb2 + postProb3) / std::max(postProb0 + postProb1, pureclip::float_min()))
              + log((postProb1 + postProb3) / std::max(postProb0 + postProb2, pureclip::float_min())));
    }
    else   // default: log posterior probability ratio (best / second best)
    {
        double secondBest = postProb2;
        secondBest = std::max(secondBest, postProb1);
        secondBest = std::max(secondBest, postProb0);
        return log(postProb3 / std::max(secondBest, pureclip::float_min()));
    }
}

// ── Write crosslink sites BED ───────────────────────────────────────────────

inline void writeStates(String<BedRecord<Bed6> > &bedRecords_sites,
                        Data &data,
                        FragmentStore<> &store,
                        unsigned contigId,
                        AppOptions &options)
{
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            for (unsigned t = 0; t < data.setObs[s][i].length(); ++t)
            {
                if (options.outputAll && data.setObs[s][i].truncCounts[t] >= 1 && !data.setObs[s][i].discard)
                {
                    BedRecord<Bed6> record;
                    record.ref = store.contigNameStore[contigId];

                    if (s == 0)
                    {
                        record.beginPos = options.crosslinkAtTruncSite
                            ? t + data.setPos[s][i]
                            : t + data.setPos[s][i] - 1;
                    }
                    else
                    {
                        record.beginPos = options.crosslinkAtTruncSite
                            ? length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1
                            : length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                    }
                    record.endPos = record.beginPos + 1;

                    std::stringstream ss;
                    ss << (int)data.states[s][i][t];
                    record.name = ss.str();
                    ss.str(""); ss.clear();

                    ss << getCrosslinkSiteScore(data.statePosteriors[s][0][i][t],
                                                data.statePosteriors[s][1][i][t],
                                                data.statePosteriors[s][2][i][t],
                                                data.statePosteriors[s][3][i][t],
                                                options.score_type);
                    record.score = ss.str();
                    ss.str(""); ss.clear();

                    record.strand = (s == 0) ? '+' : '-';

                    ss << 0 << ";"
                       << (int)data.setObs[s][i].truncCounts[t] << ";"
                       << (int)data.setObs[s][i].nEstimates[t] << ";"
                       << (double)data.setObs[s][i].kdes[t] << ";"
                       << (double)data.statePosteriors[s][3][i][t] << ";"
                       << (options.useCov_RPKM ? (double)data.setObs[s][i].rpkms[t] : 0.0) << ";"
                       << log((data.statePosteriors[s][2][i][t] + data.statePosteriors[s][3][i][t])
                            / (data.statePosteriors[s][0][i][t] + data.statePosteriors[s][1][i][t])) << ";";
                    record.data = ss.str();

                    appendValue(bedRecords_sites, record);
                }
                else if (data.setObs[s][i].discard && options.outputAll && data.setObs[s][i].truncCounts[t] >= 1)
                {
                    BedRecord<Bed6> record;
                    record.ref = store.contigNameStore[contigId];

                    if (s == 0)
                    {
                        record.beginPos = options.crosslinkAtTruncSite
                            ? t + data.setPos[s][i]
                            : t + data.setPos[s][i] - 1;
                    }
                    else
                    {
                        record.beginPos = options.crosslinkAtTruncSite
                            ? length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1
                            : length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                    }
                    record.endPos = record.beginPos + 1;

                    std::stringstream ss;
                    ss << (int)0;
                    record.name = ss.str();
                    ss.str(""); ss.clear();
                    record.score = "NA";
                    record.strand = (s == 0) ? '+' : '-';

                    ss << 0 << ";"
                       << (int)data.setObs[s][i].truncCounts[t] << ";"
                       << (int)data.setObs[s][i].nEstimates[t] << ";"
                       << (double)data.setObs[s][i].kdes[t] << ";"
                       << "NA;"
                       << (options.useCov_RPKM ? (double)data.setObs[s][i].rpkms[t] : 0.0) << ";"
                       << "NA;";
                    record.data = ss.str();

                    appendValue(bedRecords_sites, record);
                }
                else if (!data.setObs[s][i].discard && data.states[s][i][t] == 3)
                {
                    BedRecord<Bed6> record;
                    record.ref = store.contigNameStore[contigId];

                    if (s == 0)
                    {
                        record.beginPos = options.crosslinkAtTruncSite
                            ? t + data.setPos[s][i]
                            : t + data.setPos[s][i] - 1;
                    }
                    else
                    {
                        record.beginPos = options.crosslinkAtTruncSite
                            ? length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1
                            : length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                    }
                    record.endPos = record.beginPos + 1;

                    std::stringstream ss;
                    ss << (int)data.states[s][i][t];
                    record.name = ss.str();
                    ss.str(""); ss.clear();

                    ss << getCrosslinkSiteScore(data.statePosteriors[s][0][i][t],
                                                data.statePosteriors[s][1][i][t],
                                                data.statePosteriors[s][2][i][t],
                                                data.statePosteriors[s][3][i][t],
                                                options.score_type);
                    record.score = ss.str();
                    ss.str(""); ss.clear();
                    record.strand = (s == 0) ? '+' : '-';

                    ss << "[score_CL=" << getCrosslinkSiteScore(data.statePosteriors[s][0][i][t],
                                                                 data.statePosteriors[s][1][i][t],
                                                                 data.statePosteriors[s][2][i][t],
                                                                 data.statePosteriors[s][3][i][t], 1) << ";"
                       << "score_E=" << getCrosslinkSiteScore(data.statePosteriors[s][0][i][t],
                                                               data.statePosteriors[s][1][i][t],
                                                               data.statePosteriors[s][2][i][t],
                                                               data.statePosteriors[s][3][i][t], 2) << ";"
                       << "score_B=" << getCrosslinkSiteScore(data.statePosteriors[s][0][i][t],
                                                               data.statePosteriors[s][1][i][t],
                                                               data.statePosteriors[s][2][i][t],
                                                               data.statePosteriors[s][3][i][t], 3) << ";"
                       << "score_UC=" << getCrosslinkSiteScore(data.statePosteriors[s][0][i][t],
                                                                data.statePosteriors[s][1][i][t],
                                                                data.statePosteriors[s][2][i][t],
                                                                data.statePosteriors[s][3][i][t], 0) << "]";
                    record.data = ss.str();

                    appendValue(bedRecords_sites, record);
                }
            }
        }
    }
}

// ── Write binding regions BED ───────────────────────────────────────────────

inline void writeRegions(String<BedRecord<Bed6> > &bedRecords_regions,
                         Data &data,
                         FragmentStore<> &store,
                         unsigned contigId,
                         AppOptions &options)
{
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.states[s]); ++i)
        {
            for (unsigned t = 0; t < length(data.states[s][i]); ++t)
            {
                if (!data.setObs[s][i].discard && data.states[s][i][t] == 3)
                {
                    BedRecord<Bed6> record;
                    record.ref = store.contigNameStore[contigId];

                    if (s == 0)
                    {
                        record.beginPos = options.crosslinkAtTruncSite
                            ? t + data.setPos[s][i]
                            : t + data.setPos[s][i] - 1;
                    }
                    else
                    {
                        record.beginPos = options.crosslinkAtTruncSite
                            ? length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1
                            : length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                    }
                    record.endPos = record.beginPos + 1;
                    record.strand = (s == 0) ? '+' : '-';

                    unsigned prev_cs = t;
                    double score = getCrosslinkSiteScore(data.statePosteriors[s][0][i][t],
                                                          data.statePosteriors[s][1][i][t],
                                                          data.statePosteriors[s][2][i][t],
                                                          data.statePosteriors[s][3][i][t],
                                                          options.score_type);
                    double scoresSum = score;
                    std::stringstream ss_indivScores;
                    ss_indivScores << score << ';';

                    while ((t + 1) < length(data.states[s][i]) && (t + 1 - prev_cs) <= options.distMerge)
                    {
                        ++t;
                        if (data.states[s][i][t] == 3)
                        {
                            if (s == 0)
                            {
                                record.endPos = options.crosslinkAtTruncSite
                                    ? t + data.setPos[s][i] + 1
                                    : t + data.setPos[s][i];
                            }
                            else
                            {
                                record.beginPos = options.crosslinkAtTruncSite
                                    ? length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1
                                    : length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                            }

                            score = getCrosslinkSiteScore(data.statePosteriors[s][0][i][t],
                                                           data.statePosteriors[s][1][i][t],
                                                           data.statePosteriors[s][2][i][t],
                                                           data.statePosteriors[s][3][i][t],
                                                           options.score_type);
                            scoresSum += score;
                            ss_indivScores << score << ';';
                            prev_cs = t;
                        }
                    }

                    std::stringstream ss;
                    ss << scoresSum;
                    record.score = ss.str();
                    record.name = ss_indivScores.str();

                    appendValue(bedRecords_regions, record);
                }
            }
        }
    }
}

#endif // PURECLIP_IO_BED_H_
