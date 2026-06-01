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
// =======================================================================
// Author: Sabrina Krakau <krakau@molgen.mpg.de>
// =======================================================================

#ifndef APPS_HMMS_HMM_1_H_
#define APPS_HMMS_HMM_1_H_


#include <iostream>
#include <fstream>

#include "types.h"
#include "result.h"

#include "model/gamma.h"
#include "model/gamma_reg.h"
#include "model/ztbin.h"
#include "model/ztbin_reg.h"
#include <math.h>

using namespace seqan;


template <typename TGAMMA, typename TBIN>
class HMM {

public:

    uint8_t                                 K;                  // no. of sates
    String<String<String<Float> > >   initProbs;          // initial probabilities

    String<String<Observations> >           & setObs;          // workaround for partial specialization
    String<String<unsigned> >               & setPos;
    unsigned                                contigLength;
    String<String<Float> >            transMatrix;

    HMM(int K_, String<String<Observations> > & setObs_, String<String<unsigned> > & setPos_, unsigned &contigLength_): K(K_), setObs(setObs_), setPos(setPos_), contigLength(contigLength_)
    {
        // initialize transition probabilities
        resize(transMatrix, K, Exact());
        for (unsigned i = 0; i < K; ++i)
            resize(transMatrix[i], K, 0.25, Exact());

        resize(initProbs, 2, Exact());
        resize(eProbs, 2, Exact());
        resize(statePosteriors, 2, Exact());
        for (unsigned s = 0; s < 2; ++s)
        {
            resize(initProbs[s], length(setObs[s]), Exact());
            resize(eProbs[s], length(setObs[s]), Exact());
            resize(statePosteriors[s], K, Exact());
            for (unsigned k = 0; k < K; ++k)
                resize(statePosteriors[s][k], length(setObs[s]), Exact());

            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                // set initial probabilities to uniform
                resize(initProbs[s][i], K, Exact());
                for (unsigned k = 0; k < K; ++k)
                    initProbs[s][i][k] = 1.0/K;

                unsigned T = setObs[s][i].length();
                resize(eProbs[s][i], T, Exact());
                for (unsigned k = 0; k < K; ++k)
                    resize(statePosteriors[s][k][i], T, Exact());

                for (unsigned t = 0; t < T; ++t)
                {
                    resize(eProbs[s][i][t], K, Exact());
                }
            }
        }
    }
    // Copy constructor
//     HMM<TGAMMA, TBIN>(const HMM<TGAMMA, TBIN> &hmm2) {x = p2.x; y = p2.y; }

    HMM();
    ~HMM();

    Result<void> computeEmissionProbs(ModelParams<TGAMMA, TBIN> &modelParams, bool learning, AppOptions &options);
    Result<void> iForward(String<String<Float> > &alphas_1, unsigned s, unsigned i, String<String<Float> > &logA, AppOptions &options);
    Result<void> iBackward(String<String<Float> > &betas_1, unsigned s, unsigned i, String<String<Float> > &logA, AppOptions &options);
    Result<void> computeStatePosteriorsFB(AppOptions &options);
    Result<void> computeStatePosteriorsFBupdateTrans(AppOptions &options);
    bool updateTransAndPostProbs(AppOptions &options);
    Result<void> updateDensityParams(TGAMMA &gamma1, TGAMMA &gamma2, unsigned &iter, unsigned &trial, AppOptions &options);
    Result<void> updateDensityParams(ModelParams<TGAMMA, TBIN> &modelParams, AppOptions &options);
    Result<void> baumWelch(ModelParams<TGAMMA, TBIN> &modelParams, CharString learnTag, AppOptions &options);
    Result<void> applyParameters(ModelParams<TGAMMA, TBIN> &modelParams, AppOptions &/*options*/);
    void posteriorDecoding(String<String<String<uint8_t> > > &states);
    void rmBoarderArtifacts(String<String<String<uint8_t> > > &states, String<Data> &data_replicates, String<ModelParams<TGAMMA, TBIN> > &modelParams);

    // for each F/R,interval,t, state ....
    String<String<String<String<Float> > > > eProbs;           // emission/observation probabilities  P(Y_t | S_t) -> precompute for each t given Y_t = (C_t, T_t) !!!
    String<String<String<String<Float> > > > statePosteriors;  // for each k: for each covered interval string of posteriors
};


template<typename TGAMMA, typename TBIN>
HMM<TGAMMA, TBIN>::~HMM()
{
    clear(this->eProbs);
    clear(this->statePosteriors);
    clear(this->initProbs);
    clear(this->transMatrix);
   // do not touch observations
}


/////////////////////////////////////////////////////////////////
// functionalities for computations in log space
/////////////////////////////////////////////////////////////////

inline Float myLog(Float x)
{
    if (x == 0) return pureclip::float_quiet_nan();
    return log(x);
}

inline Float myExp(Float x)
{
    if (std::isnan(x)) return 0.0;
    return exp(x);
}

// log-sum-exp trick
inline Float get_logSumExp(Float &f1, Float &f2, LogSumExp_lookupTable &lookUp)
{
    if (std::isnan(f1)) return f2;
    if (std::isnan(f2)) return f1;

    if (std::isinf(f1)) return f1;
    if (std::isinf(f2)) return f2;

    return lookUp.logSumExp_add(f1, f2);
}

// log-sum-exp trick
inline Float get_logSumExp_states(Float f1, Float f2, Float f3, Float f4, LogSumExp_lookupTable &lookUp)
{
    Float sum;
    sum = get_logSumExp(f1, f2, lookUp);
    sum = get_logSumExp(sum, f3, lookUp);
    sum = get_logSumExp(sum, f4, lookUp);
    return sum;
}

// log-sum-exp trick for string
inline Float get_logSumExp(String<Float> &fs, LogSumExp_lookupTable &lookUp)
{
    Float sum = pureclip::float_quiet_nan();

    for (unsigned i = 0; i < length(fs); ++i)
        sum = get_logSumExp(sum, fs[i], lookUp);

    return sum;
}

// log-sum-exp trick for String of String
inline Float get_logSumExp(String<String<Float> > &fs, LogSumExp_lookupTable &lookUp)
{
    Float sum = pureclip::float_quiet_nan();

    for (unsigned i = 0; i < length(fs); ++i)
        for (unsigned j = 0; j < length(fs[i]); ++j)
            sum = get_logSumExp(sum, fs[i][j], lookUp);

    return sum;
}


/////////////////////////////////////////////////////////////////
// emission probabilities
/////////////////////////////////////////////////////////////////


// Unified computeEProb — replaces 4 overloads using if constexpr (C++23).
// Dispatches on whether TGAMMA/TBIN are GAMMA/GAMMA_REG and ZTBIN/ZTBIN_REG.
template<typename TEProbs, typename TSetObs, typename TGamma1, typename TGamma2, typename TBin1, typename TBin2>
bool computeEProb(TEProbs &eProbs, TSetObs &setObs, TGamma1 &gamma1, TGamma2 &gamma2, TBin1 &bin1, TBin2 &bin2, unsigned t, AppOptions &options)
{
    // ── Gamma component ──
    Float gamma1_d = 1.0;
    Float gamma2_d = 0.0;
    if (setObs.kdes[t] >= gamma1.tp)
    {
        if constexpr (std::is_same_v<TGamma1, GAMMA>)
        {
            gamma1_d = gamma1.getDensity(setObs.kdes[t]);
            gamma2_d = gamma2.getDensity(setObs.kdes[t]);
        }
        else
        {
            Float x = std::max(setObs.rpkms[t], options.minRPKMtoFit);
            Float gamma1_pred = exp(gamma1.b0 + gamma1.b1 * x);
            Float gamma2_pred = exp(gamma2.b0 + gamma2.b1 * x);
            gamma1_d = gamma1.getDensity(setObs.kdes[t], gamma1_pred, options);
            gamma2_d = gamma2.getDensity(setObs.kdes[t], gamma2_pred, options);
        }
    }

    // ── Binomial component ──
    Float bin1_d = 1.0;
    Float bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        if constexpr (std::is_same_v<TBin1, ZTBIN>)
        {
            bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);
            bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);
        }
        else
        {
            unsigned mId = setObs.motifIds[t];
            Float bin1_pred = 1.0/(1.0+exp(-bin1.b0 - bin1.regCoeffs[mId]*setObs.fimoScores[t]));
            Float bin2_pred = 1.0/(1.0+exp(-bin2.b0 - bin2.regCoeffs[mId]*setObs.fimoScores[t]));
            bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin1_pred, options);
            bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin2_pred, options);
        }
    }

    // log-space
    eProbs[0] = myLog(gamma1_d) + myLog(bin1_d);
    eProbs[1] = myLog(gamma1_d) + myLog(bin2_d);
    eProbs[2] = myLog(gamma2_d) + myLog(bin1_d);
    eProbs[3] = myLog(gamma2_d) + myLog(bin2_d);

    // check if valid
    if ((gamma1_d + gamma2_d == 0.0) || (bin1_d + bin2_d == 0.0) ||
       (std::isnan(eProbs[0]) && std::isnan(eProbs[1]) && std::isnan(eProbs[2]) && std::isnan(eProbs[3])) )
    {
        if (options.verbosity >= 2)
        {
            SEQAN_OMP_PRAGMA(critical)
                std::cout << "WARNING: emission probabilities 0.0!" << std::endl;
            SEQAN_OMP_PRAGMA(critical)
                std::cout << "       fragment coverage (kde): " << setObs.kdes[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical)
                std::cout << "       read start count: " << (int)setObs.truncCounts[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical)
                std::cout << "       estimated n: " << setObs.nEstimates[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical)
                std::cout << "       emission probability 'non-enriched' gamma: " << gamma1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical)
                std::cout << "       emission probability 'enriched' gamma: " << gamma2_d << std::endl;
            SEQAN_OMP_PRAGMA(critical)
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical)
                std::cout << "       emission probability 'crosslink' binomial: " << bin2_d << std::endl;
        }
        eProbs[0] = 0.0;
        eProbs[1] = pureclip::float_quiet_nan();
        eProbs[2] = pureclip::float_quiet_nan();
        eProbs[3] = pureclip::float_quiet_nan();
        return false;
    }
    return true;
}


template<typename TGAMMA, typename TBIN>
Result<void> HMM<TGAMMA, TBIN>::computeEmissionProbs(ModelParams<TGAMMA, TBIN> &modelParams, bool learning, AppOptions &options)
{
    // Invalidate gamma caches — parameters may have changed since last call
    modelParams.gamma1.invalidateCache();
    modelParams.gamma2.invalidateCache();
    bool stop = false;
    for (unsigned s = 0; s < 2; ++s)
    {
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1))
#endif
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            bool discardInterval = false;
            for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
            {
                if (this->setObs[s][i].kdes[t] == 0.0)
                {
                    std::cerr << "ERROR: KDE is 0.0 at i " << i << " t: " << t << std::endl;
                    SEQAN_OMP_PRAGMA(critical)
                    stop = true;
                }

                if (!computeEProb(this->eProbs[s][i][t], this->setObs[s][i], modelParams.gamma1, modelParams.gamma2, modelParams.bin1, modelParams.bin2, t, options))
                {
                    SEQAN_OMP_PRAGMA(critical)
                    discardInterval = true;
                }
            }
            if (learning && discardInterval)
            {
                SEQAN_OMP_PRAGMA(critical)
                std::cout << "ERROR: Emission probability became 0.0! This might be due to artifacts or outliers." << std::endl;
                SEQAN_OMP_PRAGMA(critical)
                if (options.verbosity >= 2)
                {
                    if (s == 0)
                        std::cout << " Interval: [" << (this->setPos[s][i]) << ", " << (this->setPos[s][i] + this->setObs[s][i].length()) << ") on forward strand." << std::endl;
                    else
                        std::cout << " Interval: [" << (this->contigLength - this->setPos[s][i] - 1) << ", " << (this->contigLength - this->setPos[s][i] - 1 + this->setObs[s][i].length()) << ") on reverse strand." << std::endl;
                }
                stop = true;
#ifdef PURECLIP_HIGH_PRECISION
                // already using long double
#else
                SEQAN_OMP_PRAGMA(critical)
                std::cout << "NOTE: Try running PureCLIP in high floating-point precision mode (long double), build with -DPURECLIP_HIGH_PRECISION." << std::endl;
#endif
            }
            else if (!learning && discardInterval)
            {
                this->setObs[s][i].discard = true;
                SEQAN_OMP_PRAGMA(critical)
                std::cout << "Warning: discarding interval on forward strand due to emission probabilities of 0.0 (set to state 'non-enriched + non-crosslink')." << std::endl;
                if (options.verbosity >= 2)
                {
                    SEQAN_OMP_PRAGMA(critical)
                    if (s == 0)
                        std::cout << " Interval [" << (this->setPos[s][i]) << ", " << (this->setPos[s][i] + this->setObs[s][i].length()) << ") on forward strand. " << std::endl;
                    else
                        std::cout << " Interval [" << (this->contigLength - this->setPos[s][i] - 1) << ", " << (this->contigLength - this->setPos[s][i] - 1 + this->setObs[s][i].length()) << ") on reverse strand." << std::endl;
                }
#ifdef PURECLIP_HIGH_PRECISION
                // already using long double
#else
                SEQAN_OMP_PRAGMA(critical)
                std::cout << "NOTE: If this happens frequently, rerun PureCLIP in high floating-point precision mode, build with -DPURECLIP_HIGH_PRECISION." << std::endl;
#endif
            }
        }
    }
    if (stop) return Result<void>("Could not compute emission probabilities");
    return Result<void>();
}



/////////////////////////////////////////////////////////////////
// forward-backward algorithm parts
/////////////////////////////////////////////////////////////////

// for one interval only
// Forward-algorithm: log-space
template<typename TGAMMA, typename TBIN>
Result<void> HMM<TGAMMA, TBIN>::iForward(String<String<Float> > &alphas_1, unsigned s, unsigned i, String<String<Float> > &logA, AppOptions &options)
{
    // NOTE
    // in log-space: alphas_1, eProbs
    // trans. probs, init porbs, state post. probs. not in log-space

    // for t = 1
    for (unsigned k = 0; k < this->K; ++k)
    {
        alphas_1[0][k] = myLog(this->initProbs[s][i][k]) + this->eProbs[s][i][0][k];    // log ? initProbs should not become 0.0!
    }

    // for t = 2:T
    for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
    {
        for (unsigned k = 0; k < this->K; ++k)
        {
            Float f1 = alphas_1[t-1][0] + logA[0][k] + this->eProbs[s][i][t][k];
            Float f2 = alphas_1[t-1][1] + logA[1][k] + this->eProbs[s][i][t][k];
            Float f3 = alphas_1[t-1][2] + logA[2][k] + this->eProbs[s][i][t][k];
            Float f4 = alphas_1[t-1][3] + logA[3][k] + this->eProbs[s][i][t][k];

            alphas_1[t][k] = get_logSumExp_states(f1, f2, f3, f4, options.lookUp);
#ifndef NDEBUG
            if (std::isinf(alphas_1[t][k]))
            {
                std::cout << "ERROR: alphas_1[" << t << "][" << k << "] is " << alphas_1[t][k] << std::endl;
                std::cout << "       f1 " << f1 << " f2 " << f2 << " f3 " << f3 << " f4 " << f4 << std::endl;
                std::cout << "       alphas_1[t-1][0] " << alphas_1[t-1][0] << " logA[0][k] " << logA[0][k] << " this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                std::cout << "       alphas_1[t-1][1] " << alphas_1[t-1][1] << " logA[1][k] " << logA[1][k] << " this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                std::cout << "       alphas_1[t-1][2] " << alphas_1[t-1][2] << " logA[2][k] " << logA[2][k] << " this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                std::cout << "       alphas_1[t-1][3] " << alphas_1[t-1][3] << " logA[3][k] " << logA[3][k] << " this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                return Result<void>("Forward: alphas became infinite");
            }
#endif
        }
    }
    return Result<void>();
}


// Backward-algorithm: log-space
template<typename TGAMMA, typename TBIN>
Result<void> HMM<TGAMMA, TBIN>::iBackward(String<String<Float> > &betas_1, unsigned s, unsigned i, String<String<Float> > &logA, AppOptions &options)
{
    unsigned T = this->setObs[s][i].length();
    // for t = T
    for (unsigned k = 0; k < this->K; ++k)
       betas_1[T - 1][k] = log(1.0);

    // for t = 2:T
    for (int t = this->setObs[s][i].length() - 2; t >= 0; --t)
    {
        for (unsigned k = 0; k < this->K; ++k)
        {
            // sum over following states
            Float f1 = betas_1[t+1][0] + logA[k][0] + this->eProbs[s][i][t+1][0];
            Float f2 = betas_1[t+1][1] + logA[k][1] + this->eProbs[s][i][t+1][1];
            Float f3 = betas_1[t+1][2] + logA[k][2] + this->eProbs[s][i][t+1][2];
            Float f4 = betas_1[t+1][3] + logA[k][3] + this->eProbs[s][i][t+1][3];

            betas_1[t][k] = get_logSumExp_states(f1, f2, f3, f4, options.lookUp);
#ifndef NDEBUG
            if (std::isinf(betas_1[t][k]))
            {
                std::cout << "ERROR: betas_1[" << t << "][" << k << "] is " << betas_1[t][k] << std::endl;
                return Result<void>("Backward: betas became infinite");
            }
#endif
        }
    }
    return Result<void>();
}


// for log-space
// interval-wise to avoid storing alpha_1 and beta_1 values for whole genome
// TODO learn 2-> 2/3 only above threshold !?
template<typename TGAMMA, typename TBIN>
Result<void> HMM<TGAMMA, TBIN>::computeStatePosteriorsFBupdateTrans(AppOptions &options)
{
    String<String<Float> > logA = this->transMatrix;
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
            logA[k_1][k_2] = log(this->transMatrix[k_1][k_2]);

    String<String<Float> > p;
    resize(p, this->K, Exact());
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
    {
        SEQAN_OMP_PRAGMA(critical)
        resize(p[k_1], this->K, Exact());
        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
            p[k_1][k_2] = 0.0;
    }
    Float p_2_2 = 0.0;     // for separate learning of trans. prob from '2' -> '2'
    Float p_2_3 = 0.0;     // for separate learning of trans. prob from '2' -> '3'

    for (unsigned s = 0; s < 2; ++s)
    {
        bool stop = false;
        // Pre-scan: find largest interval for buffer pre-allocation
        unsigned maxT = 0;
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
            if (setObs[s][i].length() > maxT)
                maxT = setObs[s][i].length();

#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel)
        {
            // ── Thread-local pre-allocated buffers ──
            String<String<Float> > alphas_1;
            resize(alphas_1, maxT, Exact());
            for (unsigned t = 0; t < maxT; ++t)
                resize(alphas_1[t], this->K, Exact());

            String<String<Float> > betas_1;
            resize(betas_1, maxT, Exact());
            for (unsigned t = 0; t < maxT; ++t)
                resize(betas_1[t], this->K, Exact());

            // Pre-allocate transition accumulators
            String<String<Float> > xis;
            resize(xis, this->K, Exact());
            String<String<Float> > p_i;
            resize(p_i, this->K, Exact());
            for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
            {
                resize(xis[k_1], this->K, Exact());
                resize(p_i[k_1], this->K, Exact());
            }

            SEQAN_OMP_PRAGMA(for schedule(dynamic, 1))
            for (unsigned i = 0; i < length(this->setObs[s]); ++i)
#else
            String<String<Float> > alphas_1;
            resize(alphas_1, maxT, Exact());
            for (unsigned t = 0; t < maxT; ++t)
                resize(alphas_1[t], this->K, Exact());
            String<String<Float> > betas_1;
            resize(betas_1, maxT, Exact());
            for (unsigned t = 0; t < maxT; ++t)
                resize(betas_1[t], this->K, Exact());
            String<String<Float> > xis;
            resize(xis, this->K, Exact());
            String<String<Float> > p_i;
            resize(p_i, this->K, Exact());
            for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
            {
                resize(xis[k_1], this->K, Exact());
                resize(p_i[k_1], this->K, Exact());
            }
            for (unsigned i = 0; i < length(this->setObs[s]); ++i)
#endif
        {
            unsigned T = setObs[s][i].length();
            // forward probabilities (uses pre-allocated alphas_1, writes rows 0..T-1)
            if (!iForward(alphas_1, s, i, logA, options))
            {
                stop = true;
                continue;
            }

            // backward probabilities (uses pre-allocated betas_1, writes rows 0..T-1)
            if (!iBackward(betas_1, s, i, logA, options))
            {
                stop = true;
                continue;
            }

            // compute state posterior probabilities
            for (unsigned t = 0; t < T; ++t)
            {
                Float f1 = alphas_1[t][0] + betas_1[t][0];
                Float f2 = alphas_1[t][1] + betas_1[t][1];
                Float f3 = alphas_1[t][2] + betas_1[t][2];
                Float f4 = alphas_1[t][3] + betas_1[t][3];

                Float norm = get_logSumExp_states(f1, f2, f3, f4, options.lookUp);

                for (unsigned k = 0; k < this->K; ++k)
                {
                    this->statePosteriors[s][k][i][t] = myExp(alphas_1[t][k] + betas_1[t][k] - norm);     // store not in log-space!
#ifndef NDEBUG
                    if (std::isnan(this->statePosteriors[s][k][i][t]) || std::isinf(this->statePosteriors[s][k][i][t]) ||
                        this->statePosteriors[s][k][i][t] < 0.0 || this->statePosteriors[s][k][i][t] > 1.0)
                    {
                        std::cout << "ERROR: state posterior probability is " << this->statePosteriors[s][k][i][t] << "." << std::endl;
                        std::cout << "       s: " << s << " i: " << i << " t: " << t << " k:" << k << std::endl;
                        std::cout << "       alphas_1[t][k]: " << alphas_1[t][k] << " betas_1[t][k]: " << betas_1[t][k] << " norm: " << norm << std::endl;
                        stop = true;
                        continue;
                    }
#endif
                }
            }

            // update initial probabilities
            for (unsigned k = 0; k < this->K; ++k)
                this->initProbs[s][i][k] = this->statePosteriors[s][k][i][0];

            // compute xi values for interval in preparation for new trans. probs
            // Zero-out pre-allocated accumulators
            for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
                for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                {
                    xis[k_1][k_2] = 0.0;
                    p_i[k_1][k_2] = 0.0;
                }
            Float p_2_2_i = 0.0;
            Float p_2_3_i = 0.0;
            //
            for (unsigned t = 1; t < T; ++t)
            {
                Float norm = pureclip::float_quiet_nan();
                for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
                {
                    for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                    {
                        xis[k_1][k_2] = alphas_1[t-1][k_1] + logA[k_1][k_2] + this->eProbs[s][i][t][k_2] + betas_1[t][k_2];
                        norm = get_logSumExp(norm, xis[k_1][k_2], options.lookUp);
                    }
                }
                for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
                    for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                        p_i[k_1][k_2] += myExp(xis[k_1][k_2] - norm);

                // learn p[2->2/3] for region over nThresholdForP
                if (options.nThresholdForTransP > 0 && setObs[s][i].nEstimates[t] >= options.nThresholdForTransP)
                {
                    p_2_2_i += myExp(xis[2][2] - norm);
                    p_2_3_i += myExp(xis[2][3] - norm);
                }
            }
            // add to global sum
            for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
                for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                    SEQAN_OMP_PRAGMA(critical)
                        p[k_1][k_2] += p_i[k_1][k_2];

            SEQAN_OMP_PRAGMA(critical)
                p_2_2 += p_2_2_i;
            SEQAN_OMP_PRAGMA(critical)
                p_2_3 += p_2_3_i;
        }
#if HMM_PARALLEL
        }  // end omp parallel
#endif
        if (stop) return Result<void>("Forward-backward update failed");
    }

    // update transition matrix
    String<String<Float> > A = this->transMatrix;
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
    {
        Float denumerator = 0.0;
        for (unsigned k_3 = 0; k_3 < this->K; ++k_3)
            denumerator += p[k_1][k_3];

        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
        {
            A[k_1][k_2] = p[k_1][k_2] / denumerator;
            if (A[k_1][k_2] <= 0.0) A[k_1][k_2] = pureclip::float_min();          // make sure not getting zero
        }
    }
    // Fix p[2->2/3] using only trans. probs. for region over nThresholdForP, while keeping sum of p[2->2] and p[2->3] constant
    if (options.nThresholdForTransP > 0)
    {
        Float sum_2_23 = A[2][2] + A[2][3];
        A[2][2] = sum_2_23 * p_2_2/(p_2_2 + p_2_3);
        A[2][3] = sum_2_23 * p_2_3/(p_2_2 + p_2_3);
    }
    // keep transProb of '2' -> '3' on min. value
    if (A[2][3] < options.minTransProbCS)
    {
        A[2][3] = options.minTransProbCS;

        if (A[3][3] < options.minTransProbCS) A[3][3] = options.minTransProbCS;
        std::cout << "NOTE: Prevented transition probability '2' -> '3' from dropping below min. value of " << options.minTransProbCS << ". Set for transitions '2' -> '3' (and if necessary also for '3'->'3') to " << options.minTransProbCS << "." << std::endl;
    }
    this->transMatrix = A;
    return Result<void>();
}


// without updating transition probabilities: log space
template<typename TGAMMA, typename TBIN>
Result<void> HMM<TGAMMA, TBIN>::computeStatePosteriorsFB(AppOptions &options)
{
    String<String<Float> > logA = this->transMatrix;
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
            logA[k_1][k_2] = log(this->transMatrix[k_1][k_2]);

    for (unsigned s = 0; s < 2; ++s)
    {
        bool stop = false;
        // Pre-scan: find largest interval for buffer pre-allocation
        unsigned maxT = 0;
        for (unsigned ii = 0; ii < length(this->setObs[s]); ++ii)
            if (setObs[s][ii].length() > maxT)
                maxT = setObs[s][ii].length();

#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel)
        {
            String<String<Float> > alphas_1;
            resize(alphas_1, maxT, Exact());
            for (unsigned t = 0; t < maxT; ++t)
                resize(alphas_1[t], this->K, Exact());
            String<String<Float> > betas_1;
            resize(betas_1, maxT, Exact());
            for (unsigned t = 0; t < maxT; ++t)
                resize(betas_1[t], this->K, Exact());

            SEQAN_OMP_PRAGMA(for schedule(dynamic, 1))
            for (unsigned i = 0; i < length(this->setObs[s]); ++i)
#else
        {
            String<String<Float> > alphas_1;
            resize(alphas_1, maxT, Exact());
            for (unsigned t = 0; t < maxT; ++t)
                resize(alphas_1[t], this->K, Exact());
            String<String<Float> > betas_1;
            resize(betas_1, maxT, Exact());
            for (unsigned t = 0; t < maxT; ++t)
                resize(betas_1[t], this->K, Exact());
            for (unsigned i = 0; i < length(this->setObs[s]); ++i)
#endif
        {
            unsigned T = setObs[s][i].length();
            if (!iForward(alphas_1, s, i, logA, options))
            {
                stop = true;
                continue;
            }
            if (!iBackward(betas_1, s, i, logA, options))
            {
                stop = true;
                continue;
            }

            for (unsigned t = 0; t < T; ++t)
            {
                Float f1 = alphas_1[t][0] + betas_1[t][0];
                Float f2 = alphas_1[t][1] + betas_1[t][1];
                Float f3 = alphas_1[t][2] + betas_1[t][2];
                Float f4 = alphas_1[t][3] + betas_1[t][3];

                Float norm = get_logSumExp_states(f1, f2, f3, f4, options.lookUp);

                for (unsigned k = 0; k < this->K; ++k)
                {
                    this->statePosteriors[s][k][i][t] = myExp(alphas_1[t][k] + betas_1[t][k] - norm);
#ifndef NDEBUG
                    if (std::isnan(this->statePosteriors[s][k][i][t])) std::cout << "ERROR: statePosterior is nan! " << std::endl;
#endif
                }
            }

            for (unsigned k = 0; k < this->K; ++k)
                this->initProbs[s][i][k] = this->statePosteriors[s][k][i][0];
        }
#if HMM_PARALLEL
        }  // end omp parallel
#endif
        if (stop) return Result<void>("Forward-backward failed");
    }
    return Result<void>();
}



bool updateDensityParams2(String<String<String<double> > > &statePosteriors1, String<String<String<double> > > &statePosteriors2, String<String<Observations> > &setObs, 
                          GAMMA &gamma1, GAMMA &gamma2, 
                          unsigned & /*iter*/, unsigned & /*trial*/,
                          AppOptions &options)
{
    if (!gamma1.updateThetaAndK(statePosteriors1, setObs, options.g1_kMin, options.g1_kMax, options))
        return false;

    if (!gamma2.updateThetaAndK(statePosteriors2, setObs, options.g2_kMin, options.g2_kMax, options))         // make sure g1k <= g2k
        return false;

    // make sure gamma1.mu < gamma2.mu
    checkOrderG1G2(gamma1, gamma2, options);
    return true;
}

bool updateDensityParams2(String<String<String<double> > > &statePosteriors1, String<String<String<double> > > &statePosteriors2, String<String<Observations> > &setObs,
                          GAMMA_REG &gamma1, GAMMA_REG &gamma2,
                          unsigned &iter, unsigned &trial,
                          AppOptions &options)
{
    if (!gamma1.updateRegCoeffsAndK(statePosteriors1, setObs, options.g1_kMin, options.g1_kMax, options))
        return false;

    double g2_kMin = options.g2_kMin;
    if (options.g1_k_le_g2_k)
        g2_kMin = std::max(gamma1.k, options.g2_kMin);

    if (!gamma2.updateRegCoeffsAndK(statePosteriors2, setObs, g2_kMin, options.g2_kMax, options))
        return false;

    // make sure gamma1.mu < gamma2.mu
    checkOrderG1G2(gamma1, gamma2, iter, trial, options);
    return true;
}


template<typename TGAMMA, typename TBIN>
Result<void> HMM<TGAMMA, TBIN>::updateDensityParams(TGAMMA &gamma1, TGAMMA &gamma2, unsigned &iter, unsigned &trial, AppOptions &options)
{
    String<String<String<double> > > statePosteriors1;
    String<String<String<double> > > statePosteriors2;
    resize(statePosteriors1, 2, Exact());
    resize(statePosteriors2, 2, Exact());
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(statePosteriors1[s], length(this->statePosteriors[s][0]), Exact());
        resize(statePosteriors2[s], length(this->statePosteriors[s][0]), Exact());
        for (unsigned i = 0; i < length(this->statePosteriors[s][0]); ++i)
        {
            resize(statePosteriors1[s][i], length(this->statePosteriors[s][0][i]), Exact());
            resize(statePosteriors2[s][i], length(this->statePosteriors[s][0][i]), Exact());
            for (unsigned t = 0; t < length(this->statePosteriors[s][0][i]); ++t)
            {
                statePosteriors1[s][i][t] = this->statePosteriors[s][0][i][t] + this->statePosteriors[s][1][i][t];
                statePosteriors2[s][i][t] = this->statePosteriors[s][2][i][t] + this->statePosteriors[s][3][i][t];
            }
        }
    }

    updateDensityParams2(statePosteriors1, statePosteriors2, this->setObs, gamma1, gamma2, iter, trial, options);

    return Result<void>();
}


template<typename TGAMMA, typename TBIN>
Result<void> HMM<TGAMMA, TBIN>::updateDensityParams(ModelParams<TGAMMA, TBIN> &modelParams, AppOptions &options)
{
    String<String<String<double> > > statePosteriors1;
    String<String<String<double> > > statePosteriors2;
    resize(statePosteriors1, 2, Exact());
    resize(statePosteriors2, 2, Exact());
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(statePosteriors1[s], length(this->statePosteriors[s][0]), Exact());
        resize(statePosteriors2[s], length(this->statePosteriors[s][0]), Exact());
        for (unsigned i = 0; i < length(this->statePosteriors[s][0]); ++i)
        {
            resize(statePosteriors1[s][i], length(this->statePosteriors[s][0][i]), Exact());
            resize(statePosteriors2[s][i], length(this->statePosteriors[s][0][i]), Exact());
            for (unsigned t = 0; t < length(this->statePosteriors[s][0][i]); ++t)
            {
                statePosteriors1[s][i][t] = this->statePosteriors[s][2][i][t];
                statePosteriors2[s][i][t] = this->statePosteriors[s][3][i][t];
            }
        }
    }

    // truncation counts
    modelParams.bin1.updateP(statePosteriors1, this->setObs, options);
    modelParams.bin2.updateP(statePosteriors2, this->setObs, options);

    // make sure bin1.p < bin2.p
    checkOrderBin1Bin2(modelParams.bin1, modelParams.bin2);

    return Result<void>();
}



// Baum-Welch
// in log-space (using log-sum-exp trick)
template<typename TGAMMA, typename TBIN>
Result<void> HMM<TGAMMA, TBIN>::baumWelch(ModelParams<TGAMMA, TBIN> &modelParams, CharString learnTag, AppOptions &options)
{
    TGAMMA prev_gamma1 = modelParams.gamma1;
    TGAMMA prev_gamma2 = modelParams.gamma2;
    TBIN prev_bin1 = modelParams.bin1;
    TBIN prev_bin2 = modelParams.bin2;
    unsigned trial = 0;
    for (unsigned iter = 0; iter < options.maxIter_bw; ++iter)
    {
        std::cout << ".. " << iter << "th iteration " << std::endl;
        std::cout << "                        computeEmissionProbs() " << std::endl;
        auto eprobsRes = computeEmissionProbs(modelParams, true, options);
        if (!eprobsRes)
        {
            std::cerr << "ERROR: Could not compute emission probabilities! " << eprobsRes.error() << std::endl;
            return Result<void>("Emission probability computation failed");
        }
        std::cout << "                        computeStatePosteriorsFB() " << std::endl;
        auto fbRes = computeStatePosteriorsFBupdateTrans(options);
        if (!fbRes)
        {
            std::cerr << "ERROR: Could not compute forward-backward algorithm! " << fbRes.error() << std::endl;
            return Result<void>("Forward-backward failed");
        }

        std::cout << "                        updateDensityParams() " << std::endl;

        if (learnTag == "LEARN_BINOMIAL")
        {
            auto upRes = updateDensityParams(modelParams, options);
            if (!upRes)
            {
                std::cerr << "ERROR: Could not update parameters! " << upRes.error() << std::endl;
                return Result<void>("Density param update failed");
            }
        }
        else
        {
            auto upResG = updateDensityParams(modelParams.gamma1, modelParams.gamma2, iter, trial, options);
            if (!upResG)
            {
                std::cerr << "ERROR: Could not update parameters! " << upResG.error() << std::endl;
                return Result<void>("Density param update failed");
            }
            if (trial > 10)
            {
                std::cerr << "ERROR: Could not learn gamma parameters, exceeded max. number of reseedings! " << std::endl;
                return Result<void>("Exceeded max reseedings");
            }

        }

        if (learnTag == "LEARN_GAMMA" && checkConvergence(modelParams.gamma1, prev_gamma1, options) && checkConvergence(modelParams.gamma2, prev_gamma2, options) )
        {
            std::cout << " **** Convergence ! **** " << std::endl;
            break;
        }
        else if (learnTag != "LEARN_GAMMA" && checkConvergence(modelParams.bin1, prev_bin1, options) && checkConvergence(modelParams.bin2, prev_bin2, options) )
        {
            std::cout << " **** Convergence ! **** " << std::endl;
            break;
        }
        prev_gamma1 = modelParams.gamma1;
        prev_gamma2 = modelParams.gamma2;
        prev_bin1 = modelParams.bin1;
        prev_bin2 = modelParams.bin2;

        myPrint(modelParams.gamma1);
        myPrint(modelParams.gamma2);

        std::cout << "*** Transition probabilitites ***" << std::endl;
        for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
        {
            std::cout << "    " << k_1 << ": ";
            for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                std::cout << this->transMatrix[k_1][k_2] << "  ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
        if (learnTag != "LEARN_GAMMA")
         {
            myPrint(modelParams.bin1);
            myPrint(modelParams.bin2);
        }
    }
    return Result<void>();
}


template<typename TGAMMA, typename TBIN>
Result<void> HMM<TGAMMA, TBIN>::applyParameters(ModelParams<TGAMMA, TBIN> &modelParams, AppOptions &options)
{
    auto eprobsRes = computeEmissionProbs(modelParams, false, options);
    if (!eprobsRes)
    {
        std::cerr << "ERROR: Could not compute emission probabilities! " << eprobsRes.error() << std::endl;
        return Result<void>("Emission probability computation failed");
    }
    auto fbRes = computeStatePosteriorsFB(options);
    if (!fbRes)
    {
        std::cerr << "ERROR: Could not compute forward-backward algorithm! " << fbRes.error() << std::endl;
        return Result<void>("Forward-backward failed");
    }
    return Result<void>();
}



template<typename TGAMMA, typename TBIN>
void HMM<TGAMMA, TBIN>::posteriorDecoding(String<String<String<uint8_t> > > &states)
{
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(states[s], length(this->setObs[s]), Exact());
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            if (!this->setObs[s][i].discard)
            {
                resize(states[s][i], this->setObs[s][i].length(), Exact());
                for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
                {
                    Float max_p = 0.0;
                    unsigned max_k = 0;
                    for (unsigned k = 0; k < this->K; ++k)
                    {
                        if (this->statePosteriors[s][k][i][t] > max_p)
                        {
                            max_p = this->statePosteriors[s][k][i][t];
                            max_k = k;
                        }
                    }
                    states[s][i][t] = max_k;
                }
            }
        }
    }
}


void rmBoarderArtifacts2(String<String<String<uint8_t> > > &states,
                         String<String<Observations> > & /*setObs*/,
                         String<Data> &data_replicates,
                         String<ModelParams<GAMMA_REG, ZTBIN> > &modelParams)
{
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data_replicates[0].setObs[s]); ++i)
        {
            if (!data_replicates[0].setObs[s][i].discard)
            {
                for (unsigned t = 0; t < data_replicates[0].setObs[s][i].length(); ++t)
                {
                    if (states[s][i][t] >= 2)
                    {
                        double x1 = data_replicates[0].setObs[s][i].rpkms[t];   // same for all replicates
                        bool rm = true;
                        // rm if not > gamma1 expected value in all replicates
                        for (unsigned rep = 0; rep < length(data_replicates); ++rep)
                        {
                            double gamma1_pred = exp(modelParams[rep].gamma1.b0 + modelParams[rep].gamma1.b1 * x1);
                            if (data_replicates[rep].setObs[s][i].kdes[t] > gamma1_pred) rm = false;
                        }
                        if (rm) states[s][i][t] -= 2;
                    }
                }
            }
        }
    }
}

void rmBoarderArtifacts2(String<String<String<uint8_t> > > &states,
                         String<String<Observations> > & /*setObs*/,
                         String<Data> &data_replicates,
                         String<ModelParams<GAMMA_REG, ZTBIN_REG> > &modelParams)
{
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data_replicates[0].setObs[s]); ++i)
        {
            if (!data_replicates[0].setObs[s][i].discard)
            {
                for (unsigned t = 0; t < data_replicates[0].setObs[s][i].length(); ++t)
                {
                    if (states[s][i][t] >= 2)
                    {
                        double x1 = data_replicates[0].setObs[s][i].rpkms[t];   // same for all replicates
                        bool rm = true;
                        // rm if not > gamma1 expected value in all replicates
                        for (unsigned rep = 0; rep < length(data_replicates); ++rep)
                        {
                            double gamma1_pred = exp(modelParams[rep].gamma1.b0 + modelParams[rep].gamma1.b1 * x1);
                            if (data_replicates[rep].setObs[s][i].kdes[t] > gamma1_pred) rm = false;
                        }
                        if (rm) states[s][i][t] -= 2;
                    }
                }
            }
        }
    }
}

// do nothing
void rmBoarderArtifacts2(String<String<String<uint8_t> > > & /*states*/,
                         String<String<Observations> > & /*setObs*/,
                         String<Data> & /*data_replicates*/,
                         String<ModelParams<GAMMA, ZTBIN> > & /*modelParams*/){}

void rmBoarderArtifacts2(String<String<String<uint8_t> > > & /*states*/,
                         String<String<Observations> > & /*setObs*/,
                         String<Data> & /*data_replicates*/,
                         String<ModelParams<GAMMA, ZTBIN_REG> > & /*modelParams*/){}


// for GLM with input signal:
// when using free gamma shapes, i.e. gamma1.k can be > gamma2.k
// make sure sites with fragment coverage (KDE) below gamma1.mean are classified as 'non-enriched'
template<typename TGAMMA, typename TBIN>
void HMM<TGAMMA, TBIN>::rmBoarderArtifacts(String<String<String<uint8_t> > > &states,
                                           String<Data> &data_replicates,
                                           String<ModelParams<TGAMMA, TBIN> > &modelParams)
{
    rmBoarderArtifacts2(states, this->setObs, data_replicates, modelParams);
}


// BED-output functions moved to io/bed.h
#include "io/bed.h"


template<typename TGAMMA, typename TBIN>
void myPrint(HMM<TGAMMA, TBIN> &hmm)
{
    std::cout << "*** Transition probabilitites ***" << std::endl;
    for (unsigned k_1 = 0; k_1 < hmm.K; ++k_1)
    {
        std::cout << "    " << k_1 << ": ";
        for (unsigned k_2 = 0; k_2 < hmm.K; ++k_2)
            std::cout << hmm.transMatrix[k_1][k_2] << "  ";
        std::cout << std::endl;
    }
}


template<typename TOut>
void printParams(TOut &out, String<String<Float> > &transMatrix)
{
    out << "Transition probabilities:" << std::endl;
    for (unsigned k_1 = 0; k_1 < 4; ++k_1)
    {
        for (unsigned k_2 = 0; k_2 < 4; ++k_2)
            out << transMatrix[k_1][k_2] << "\t";
        out << std::endl;
    }
}

#endif
