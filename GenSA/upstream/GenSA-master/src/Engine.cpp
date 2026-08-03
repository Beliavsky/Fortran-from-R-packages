// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################

#include "Utils.h"
#include "Engine.h"
#include <stdexcept>

using namespace std;

Engine::Engine()
{
}

Engine::~Engine()
{
}

int Engine::initialize()
{
    // For Tracing
    std::vector < std::string > keys;
    keys.push_back("currentEnergy");
    keys.push_back("minEnergy");
    keys.push_back("nSteps");
    keys.push_back("temperature");
    tracer_.clear();
    tracer_.setKeyList(keys);

    if (isVerbose())
    {
        Rprintf("Initialization...\n");
    }
    // Init x related vectors
    try {
        xRange_.resize(x_.size());
        xBackup_.resize(x_.size());
        xMini_.resize(x_.size());
        xBuffer_.resize(x_.size());
        g_.resize(x_.size());
    }
    catch (std::length_error& le) {
        Rprintf("Engine: Length error: %s\n",le.what());
    }

    itSoftMax_ = x_.size() * 6;
    factr_ = 1000;
    pgTol_ = 1.e-6;
    reps_ = 1.e-6;

    nbFctCall_ = 0;
    idum_ = -100377;

    // Check markov chain length
    if (0 != markovLength_ % x_.size())
    {
        Rprintf(
                "LMarkov should be size of 'x' (recommended) or 2*n or 3*n ... since component.change is 1\n");
        return -1;
    }
    //	if (lsEnd_)
    //	{
    //		markovLength_ = 200 * x_.size();
    //		if (markovLength_ < 1000)
    //		{
    //			markovLength_ = 1000;
    //		}
    //		else if (markovLength_ > 10000)
    //		{
    //			markovLength_ = 10000;
    //		}
    //	}

    if (isVerbose())
    {
        Rprintf("LMarkov= %i\n", markovLength_);
    }

    for (unsigned int i = 0; i < x_.size(); ++i)
    {
        xRange_[i] = upper_[i] - lower_[i];
    }

    if (isVerbose())
    {
        Rprintf("xrange: ");
        printVect(xRange_);
    }

    // Check if starting point is in constraint
    bool inConstraint = true;
    bool initError = true;
    while (initError)
    {
        if (inConstraint)
        {
            if (hasConstraint_)
            {
                inConstraint = judgeConstraint();
                while (!inConstraint)
                {
                    coordin(idum_, x_);
                    inConstraint = judgeConstraint();
                }
            }
        }
        if (isVerbose())
        {
            Rprintf("The random intial x coordinates:\n");
            printVect(x_);
        }
        energy(x_);
        if (isVerbose())
        {
            Rprintf("The energy of initial x = %.10g\n", etot_);
        }

        if (etot_ >= BIG_VALUE)
        {
            if (isVerbose())
            {
                Rprintf("x: ");
                printVect(x_);
                Rprintf(
                        " give NaN, NA, or inf, generating new starting point\n");
            }
            double rd = 0;
            for (unsigned int i=0; i < x_.size(); ++i)
            {
                // lower + runif(Rf_length(lower))*(upper-lower)
                rd = Utils::ran2(&idum_);
                x_[i] = lower_[i] + rd * (upper_[i] - lower_[i]);
            }
        }
        else
        {
            initError = false;
        }
    }

    return 0;
}

int Engine::startSearch()
{
    if (isVerbose())
    {
        Rprintf("Starting...\n");
    }
    double temQa;
    double visit, a, b;
    int itNew = 0, itDev;
    indTrace_ = 0.;
    double s1, s, t1, t2, r, pqa, pqa1;
    bool inConstraint = false;
    bool eMini_NotChanged = true;
    double eMiniMarkov = 0;
    int indexNoEminiUpdate = 0;
    int indexTolEminiUpdate = 1000;
    dVec xMiniMarkov(x_.size());

    if (getIsSimpleFunction())
    {
        indexTolEminiUpdate = x_.size();
    }

    //	if (x_.size() <= 2)
    //	{
    //		indexTolEminiUpdate = 3;
    //	}
    //	else if (x_.size() > 2 && x_.size() <= 4)
    //	{
    //		indexTolEminiUpdate = 4 * x_.size();
    //	}
    //	else if (x_.size() > 4 && x_.size() <= 10)
    //	{
    //		indexTolEminiUpdate = 4 * x_.size();
    //	}

    startTime_ = clock();
    eMini_ = etot_;
    xMini_ = x_;
    etot0_ = etot_;

    // Initialize etot0 and temp
    if (!lsEnd_)
    {
        etot_ = lsEnergy(x_);
        if (etot_ < eMini_)
        {
            eMini_ = etot_;
            xMini_ = x_;
            if (useTraceMat())
            {
                tracer_.addValue("currentEnergy", etot0_);
                tracer_.addValue("minEnergy", eMini_);
                tracer_.addValue("nSteps", itNew);
                tracer_.addValue("temperature", temSta_);
            }
        }

    }
    if (etot_ < eMini_)
    {
        eMini_ = etot_;
        xMini_ = x_;
    }
    etot0_ = etot_;
    tem_ = temSta_;

    if (isVerbose())
    {
        Rprintf("The transformed xinitial x: \n");
        printVect(x_);
        Rprintf("The energy of  transformed initial x = %.10g\n", etot_);
    }

    if (isVerbose())
    {
        Rprintf("At the beginning, etot0= %.10g\n", etot0_);
        Rprintf("Emini= %.10g\n", eMini_);
        Rprintf("Current x: ");
        printVect(x_);
        Rprintf("Current xmini: ");
        printVect(xMini_);
    }
    if (checkStoping())
    {
        stopSearch();
        return 0;
    }

    if (isVerbose())
    {
        Rprintf("Number of function call: %i\n", nbFctCall_);
    }
    int stepRecord = 0;
L2435:
    // Main loop
    for (int i = 0; i < maxStep_; ++i)
    {
        itNew = i + 1;
        indTrace_++;
        s1 = (double) itNew;
        s = s1 + 1.;
        t1 = exp((qv_ - 1.) * log(2.)) - 1.;
        t2 = exp((qv_ - 1.) * log(s)) - 1.;
        tem_ = temSta_ * t1 / t2;
        stepRecord += 1;
        if (stepRecord == maxStep_)
        {
            break;
        }
        if (tem_ < temRestart_)
        {
            //printf("*\n");
            goto L2435;
        }
        temQa = tem_ / (double) itNew;

        indexNoEminiUpdate++;

        // Markov loop
        for (unsigned int j = 0; j < (unsigned) markovLength_; ++j)
        {
            if (j == 0)
            {
                eMini_NotChanged = true;
            }
            if (j == 0 && i == 0)
            {
                eMini_NotChanged = false;
            }

            xBackup_ = x_;
            inConstraint = false;
            while (!inConstraint)
            {
                // Loop on coordinates
                if (j < x_.size())
                {
                    for (unsigned int k = 0; k < x_.size(); ++k)
                    {
                        if (isVerbose())
                        {
                            Rprintf("IDUM before visit: %ld\n", idum_);
                        }
                        visit = visita(&qv_, &tem_, &idum_);
                        if (visit > 1.e8)
                        {
                            visit = 1.e8 * Utils::ran2(&idum_);
                        }
                        else if (visit < -1.e8)
                        {
                            visit = -1.e8 * Utils::ran2(&idum_);
                        }

                        x_[k] = visit + xBackup_[k];
                        a = x_[k] - lower_[k];
                        b = Utils::dMod(&a, &xRange_[k]) + xRange_[k];
                        x_[k] = Utils::dMod(&b, &xRange_[k]) + lower_[k];
                        if (fabs(x_[k] - lower_[k]) < 1.e-10)
                        {
                            x_[k] += 1.e-10;
                        }
                        if (isVerbose())
                        {
                            Rprintf(
                                    "visit: %.10g a: %.10g b: %.10g idum: %ld x: %.10g\n",
                                    visit, a, b, idum_, x_[k]);
                        }
                    } // end coordinates loop
                }
                else
                {
                    // Now change only one component at a time
                    visit = visita(&qv_, &tem_, &idum_);
                    if (visit > 1.e8)
                    {
                        visit = 1.e8 * Utils::ran2(&idum_);
                    }
                    else if (visit < -1.e8)
                    {
                        visit = -1.e8 * Utils::ran2(&idum_);
                    }
                    int index = j - x_.size();
                    x_[index] = visit + xBackup_[index];
                    a = x_[index] - lower_[index];
                    b = Utils::dMod(&a, &xRange_[index]) + xRange_[index];
                    x_[index] = Utils::dMod(&b, &xRange_[index])
                        + lower_[index];
                    if (fabs(x_[index] - lower_[index]) < 1.e-10)
                    {
                        x_[index] += 1.e-10;
                    }
                }
                if (isVerbose())
                {
                    Rprintf("\ntem: %.10g temqa: %.10g itnew: %d markov j=%d\n",
                            tem_, temQa, itNew, j);
                    Rprintf("fx are: ");
                    printVect(xBackup_);
                    Rprintf("x are: ");
                    printVect(x_);
                }

                if (hasConstraint_)
                {
                    inConstraint = judgeConstraint();
                }
                else
                {
                    inConstraint = true;
                }
                if (inConstraint)
                {
                    if (lsEnd_)
                    {
                        if (isVerbose())
                        {
                            Rprintf("Calling energy\n");
                        }
                        energy(x_);
                    }
                    else
                    {
                        if (isVerbose())
                        {
                            Rprintf("Calling lsEnergy\n");
                        }
                        etot_ = lsEnergy(x_);
                    }
                    if (isVerbose())
                    {
                        Rprintf("Before judge, etot0= %.10g etot= %.10g\n",
                                etot0_, etot_);
                    }
                    if (etot_ < etot0_)
                    {
                        etot0_ = etot_;
                        if (isVerbose())
                        {
                            Rprintf("etot is smaller than etot0\n");
                        }
                        if (etot_ < eMini_)
                        {
                            eMini_ = etot_;
                            xMini_ = x_;
                            eMini_NotChanged = false;
                            indexNoEminiUpdate = 0;
                        }
                    }
                    else
                    {
                        r = Utils::ran2(&idum_);
                        pqa1 = (qa_ - 1.) * (etot_ - etot0_) / temQa + 1.;
                        /* write(*,*)' etot0=',etot0,', etot=',etot,', pqa1=',pqa1 ! r */
                        if (pqa1 < 0.)
                        {
                            pqa = 0.;
                        }
                        else
                        {
                            pqa = exp(log(pqa1) / (1. - qa_));
                        }
                        if (isVerbose())
                        {
                            Rprintf("pqa= %.10g r= %.10g \n", pqa, r);
                        }
                        if (r > pqa)
                        {
                            x_ = xBackup_;
                        }
                        else
                        {
                            etot0_ = etot_;
                        }
                        // if (useTraceMat())
                        // {
                        //     tracer_.addValue("currentEnergy", etot0_);
                        //     tracer_.addValue("minEnergy", eMini_);
                        //     tracer_.addValue("nSteps", indTrace_);
                        //     tracer_.addValue("temperature", tem_);
                        // }
                    }

                    if (checkStoping())
                    {
                        stopSearch();
                        return 0;
                    }
                } // end else hasConstraint
            } // end while !inconstraint
            if (indexNoEminiUpdate >= indexTolEminiUpdate - 1)
            {
                if (j == 0)
                {
                    eMiniMarkov = etot0_;
                    std::copy(x_.begin(), x_.end(), xMiniMarkov.begin());
                }
                else
                {
                    if (etot0_ < eMiniMarkov)
                    {
                        eMiniMarkov = etot0_;
                        std::copy(x_.begin(), x_.end(), xMiniMarkov.begin());
                    }
                }
            }
        } // end markov chain loop
        if (lsEnd_)
        {
            if (!eMini_NotChanged)
            {
                dVec temp(x_.size());
                std::copy(xMini_.begin(), xMini_.end(), temp.begin());
                //Rprintf("Xmini:\n");
                //printVect(xMini_);
                //				if (isUserVerbose())
                //				{
                //					Rprintf("Before lsEnergy, itNew: %d eTemp: %.15g x: %.15g %15g\n", itNew, eMini_, temp[0], temp[1]);
                //				}
                double eTemp = lsEnergy(temp);
                //				if (isUserVerbose())
                //				{
                //					Rprintf("lsEnergy called, itNew: %d eTemp: %.15g x: %.15g %.15g\n", itNew, eTemp, temp[0], temp[1]);
                //				}

                if (eTemp < eMini_)
                {
                    if (isUserVerbose())
                    {
                        Rprintf("It: %d, obj value (lsEnd): %.10g indTrace: %d\n", itNew, eTemp, indTrace_);
                    }
                    std::copy(temp.begin(), temp.end(), xMini_.begin());
                    eMini_ = eTemp;
                    indexNoEminiUpdate = 0;
                    if (useTraceMat())
                    {
                        tracer_.addValue("currentEnergy", etot0_);
                        tracer_.addValue("minEnergy", eMini_);
                        tracer_.addValue("nSteps", indTrace_);
                        tracer_.addValue("temperature", tem_);
                    }
                }
            }
            if (indexNoEminiUpdate >= indexTolEminiUpdate)
            {
                //Rprintf("Before lsEnergy, itNew: %d x: %.15g %.15g\n", itNew, xMiniMarkov[0], xMiniMarkov[1]);
                eMiniMarkov = lsEnergy(xMiniMarkov);
                //Rprintf("lsEnergy called, itNew: %d eMiniMarkov: %.15g x: %.15g %.15g\n", itNew, eMiniMarkov, xMiniMarkov[0], xMiniMarkov[1]);
                indexNoEminiUpdate = 0;
                indexTolEminiUpdate = x_.size();
                if (eMiniMarkov < eMini_)
                {
                    std::copy(xMiniMarkov.begin(), xMiniMarkov.end(),
                            xMini_.begin());
                    eMini_ = eMiniMarkov;
                    if (useTraceMat())
                    {
                        tracer_.addValue("currentEnergy", etot0_);
                        tracer_.addValue("minEnergy", eMini_);
                        tracer_.addValue("nSteps", indTrace_);
                        tracer_.addValue("temperature", tem_);
                    }
                    if (isUserVerbose())
                    {
                        Rprintf("It: %d, obj value: %.10g indTrace: %d\n", itNew, eMini_, indTrace_);
                    }
                    if (checkStoping())
                    {
                        stopSearch();
                        return 0;
                    }
                }
            }
        }

        itDev = itNew % interval_;
        if (0 == itDev)
        {
            if (isVerbose())
            {
                Rprintf(">After one search %d %.10g %.10g %.10g <-- Emini\n",
                        itNew, tem_, etot0_, eMini_);
                Rprintf("Current x: ");
                printVect(x_);
                Rprintf("\nCurrent xmini: ");
                printVect(xMini_);
                Rprintf(
                        "\n__________________________________________________\n");
            }
        }
    } // end step loop
    stopSearch();
    if (isUserVerbose())
    {
        Rprintf("Algorithm reached max number of iterations.\n");
    }
    return 0;
}

bool Engine::checkStoping()
{
    bool canStop = false;
    if (knowRealEnergy_)
    {
        canStop = eMini_ <= realEnergyThreshold_;
        if (canStop)
        {
            if (isUserVerbose())
            {
                Rprintf(
                        "Have got accurate energy %.10g <= %.10g in smooth search\n",
                        eMini_, realEnergyThreshold_);
            }
            return true;
        }
    }
    endTime_ = clock();
    timeSpan_ = (double) (endTime_ - startTime_) / CLOCKS_PER_SEC;
    if (timeSpan_ > maxTime_)
    {
        if (isUserVerbose())
        {
            Rprintf("timeSpan = %.10g maxTime = %.10g\n", timeSpan_, maxTime_);
        }
        return true;
    }
    if (nbFctCall_ >= maxFctCall_)
    {
        if (isUserVerbose())
        {
            Rprintf("Stop. Nb function call=%d max function call=%d.\n",
                    nbFctCall_, maxFctCall_);
        }
        return true;
    }

    if (useTraceMat() && noImprovementStop_ != -1)
    {
        if (tracer_.getTracerLength() == 0)
        {
            return false;
        }
        double last_iter_best_value = tracer_.getLastValue("nSteps");
        if ((indTrace_ - last_iter_best_value) >= noImprovementStop_ &&
            noImprovementStop_ != -1)
        {
            if (isUserVerbose())
            {
                Rprintf("No improvement in %i iterations\n", noImprovementStop_);
            }
            return true;
        }
    }
    return false;
}

void Engine::stopSearch()
{
    if (isUserVerbose())
    {
        Rprintf("Emini is: %.10g\n", eMini_);
        Rprintf("xmini are:\n");
        printVect(xMini_);
    }
    endTime_ = clock();
    timeSpan_ = (double) (endTime_ - startTime_) / CLOCKS_PER_SEC;

    if (isUserVerbose())
    {
        Rprintf("Totally it used %.10g secs\n", timeSpan_);
        Rprintf("No. of function call is: %d\n", nbFctCall_);
    }
}

void Engine::coordin(long int idum, dVec& x)
{
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        x[i] = Utils::ran2(&idum) * xRange_[i] + lower_[i];
    }
}

int Engine::energy(const dVec& x)
{
    double penalty = 0.;
    double deltaEnergy = 0.;
    bool inConstraint = false;
    double kSpring = 1.e8;

    penalty = 0.;

    if (hasConstraint_)
    {
        inConstraint = judgeConstraint();
        if (!inConstraint)
        {
            etot_ = BIG_VALUE;
            return 0;
        }
    }
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        if (x[i] >= lower_[i] && x[i] <= upper_[i])
        {
            deltaEnergy = 0.;
        }
        else
        {
            if (x[i] < lower_[i])
            {
                deltaEnergy = fabs(x[i] - lower_[i]) * kSpring;
            }
            if (x[i] > upper_[i])
            {
                deltaEnergy = fabs(x[i] - upper_[i]) * kSpring;
            }
        }
        penalty += deltaEnergy;
    }
    etot_ = fn(x);
    nbFctCall_++;
    etot_ = etot_ + penalty;

    if (ISNA(etot_) || !R_FINITE(etot_))
    {
        etot_ = BIG_VALUE;
    }
    return 0;
}

double Engine::lsEnergy(dVec& x)
{
    // Saving given x to the buffer
    //	Rprintf("Xvector is:\n");
    //	printVect(x);
    std::copy(x.begin(), x.end(), xBuffer_.begin());
    //	Rprintf("Xbuffer is:\n");
    //	printVect(xBuffer_);

    // Filling XBigBuffer with ones
    //	for (unsigned int i = 0; i < x_.size(); ++i)
    //	{
    //		xBigBuffer_[i * x_.size() + i] = 1.;
    //	}

    // Switch between methods
    if (method_ == Engine::SMOOTH)
    {
        smoothSearch();
    }
    else
    {
        hardSearch();
    }
    std::copy(xBuffer_.begin(), xBuffer_.end(), x.begin());
    return fValue_;
}

double Engine::visita(double* q, double* temp, long int *idum)
{
    static double x, y, pi, den;
    static double fator1, fator2, fator3, fator4, fator5, fator6, sigmax;
    double ret_val, d__1;
    pi = asin(1.) * 2.;
    fator1 = exp(log(*temp) / (*q - 1.));
    fator2 = exp((4. - *q) * log(*q - 1.));
    fator3 = exp((2. - *q) * log(2.) / (*q - 1.));
    fator4 = sqrt(pi) * fator1 * fator2 / (fator3 * (3. - *q));
    fator5 = 1. / (*q - 1.) - .5;
    // calculates the gamma function using the reflection formula for
    // 0<arg<1
    d__1 = 2. - fator5;
    fator6 = pi * (1. - fator5) / sin(pi * (1. - fator5))
        / exp(Utils::yyGaml(&d__1));
    sigmax = exp(-(*q - 1.) * log(fator6 / fator4) / (3. - *q));
    x = sigmax * Utils::yyGas(idum);
    y = Utils::yyGas(idum);
    den = exp((*q - 1.) * log((fabs(y))) / (3. - *q));
    ret_val = x / den;
    return ret_val;
}

double Engine::fObjective(const dVec& x)
{
    std::copy(x.begin(), x.end(), x_.begin());
    energy(x_);
    return etot_;
}

double Engine::fn(const dVec& x)
{
    SEXP x4R, val;
    double res = 0;
    if (isVerbose())
    {
        Rprintf(".");
    }

    // Allocate vector for R which is size of the vector in the R context.
    PROTECT(x4R = Rf_allocVector(REALSXP, x.size()));
    if (rEnv_->xNames)
        Rf_setAttrib(x4R, R_NamesSymbol, rEnv_->xNames);

    for (unsigned int i = 0; i < x.size(); i++)
    {
        if (!R_FINITE(x[i]))
        {
            if (isVerbose())
            {
                Rprintf("x[%i] is NAN: %.10g\n", i, x[i]);
            }
            REAL(x4R)[i] = 0.;
        }
        else
        {
            REAL(x4R)[i] = x[i];
        }
    }

    SETCADR(rEnv_->R_fn, x4R);
    val = Rf_eval(rEnv_->R_fn, rEnv_->R_env);
    res = REAL(val)[0];
    UNPROTECT(1);
    return res;
}

bool Engine::judgeConstraint()
{
    SEXP x4R, val;
    int res;

    // Allocate vector for R which is size of the vector in the R context.
    PROTECT(x4R = Rf_allocVector(REALSXP, x_.size()));
    if (rEnv_->xNames)
        Rf_setAttrib(x4R, R_NamesSymbol, rEnv_->xNames);

    for (unsigned int i = 0; i < x_.size(); i++)
    {
        if (!R_FINITE(x_[i]))
        {
            if (isVerbose())
            {
                Rprintf("x[%i] is NAN: %.10g\n", i, x_[i]);
            }
            REAL(x4R)[i] = 0;
        }
        else
        {
            REAL(x4R)[i] = x_[i];
        }
    }

    SETCADR(rEnv_->R_jc, x4R);
    val = Rf_eval(rEnv_->R_jc, rEnv_->R_env);
    res = LOGICAL(val)[0];
    UNPROTECT(1);

    return res;
}

int Engine::hardSearch()
{
    SEXP s, t, val;
    SEXP uiMatrix, ciVector;
    SEXP thetaVector;
    SEXP xlow;
    SEXP xhigh;

    //int lsConvergence = 0;
    int counts = 0;
    double mu;
    int xSize = x_.size();

    PROTECT(uiMatrix = Rf_allocMatrix(REALSXP, xSize * 2, xSize));
    //protect 1
    PROTECT(ciVector = Rf_allocVector(REALSXP, xSize * 2));
    //protect 2
    PROTECT(thetaVector = Rf_allocVector(REALSXP, xSize));
    //protect 3
    PROTECT(xlow = Rf_allocVector(REALSXP, xSize));
    // protect 4
    PROTECT(xhigh = Rf_allocVector(REALSXP, xSize));
    // protect 5
    mu = 1.e-4;
    // Initialize ui with zeros
    for (int i = 0; i < xSize * 2; ++i)
    {
        for (int j = 0; j < xSize; ++j)
        {
            REAL(uiMatrix)[i * xSize + j] = 0.;
        }
    }

    for (int i = 0; i < xSize; ++i)
    {
        // Initialize theta
        REAL(thetaVector)[i] = xBuffer_[i];
        // Initialize ci
        REAL(ciVector)[i * 2] = lower_[i];
        REAL(ciVector)[i * 2 + 1] = -upper_[i];
        REAL(uiMatrix)[i * 2 * xSize + i * 2] = 1.0;
        REAL(uiMatrix)[i * 2 * xSize + i * 2 + 1] = -1.0;
        REAL(xlow)[i] = lower_[i];
        REAL(xhigh)[i] = upper_[i];
    }

    PROTECT(t = s = Rf_allocLang(8));
    SETCAR(t, Rf_install("LSE"));
    t = CDR(t);
    SETCAR(t, thetaVector);
    SET_TAG(t, Rf_install("theta"));
    t = CDR(t);
    SETCAR(t, uiMatrix);
    SET_TAG(t, Rf_install("ui"));
    t = CDR(t);
    SETCAR(t, ciVector);
    SET_TAG(t, Rf_install("ci"));
    t = CDR(t);
    SETCAR(t, Rf_ScalarReal(mu));
    SET_TAG(t, Rf_install("mu"));
    t = CDR(t);
    SETCAR(t, xlow);
    SET_TAG(t, Rf_install("xlow"));
    t = CDR(t);
    SETCAR(t, xhigh);
    SET_TAG(t, Rf_install("xhigh"));
    t = CDR(t);
    SETCAR(t, Rf_ScalarInteger(nbFctCall_));
    SET_TAG(t, Rf_install("count"));

    for (unsigned int i = 0; i < xBuffer_.size(); ++i)
    {
        if (xBuffer_[i] < lower_[i] || xBuffer_[i] > upper_[i])
        {
            Rprintf("PROBLEM WITH x(%d):\n", i);
            printVect(xBuffer_);
        }
    }

    val = Rf_eval(s, rEnv_->R_env);
    fValue_ = REAL(VECTOR_ELT(val, 0))[0];

    //lsConvergence = INTEGER(VECTOR_ELT(val, 1))[0];

    for (unsigned int i = 0; i < xBuffer_.size(); ++i)
    {
        xBuffer_[i] = REAL(VECTOR_ELT(val, 2))[i];
    }
    counts = INTEGER(VECTOR_ELT(val, 3))[0];
    nbFctCall_ = counts;
    UNPROTECT(6);
    return 0;
}

int Engine::smoothSearch()
{
    char task[60];
    std::vector<int> nbd(xBuffer_.size());
    double f, dsave[29], *wa;
    int tr = -1, iter = 0, *iwa, isave[44], lsave[4];
    int xSize = xBuffer_.size();
    int m = 5;
    bool canstop = 0;

    wa = (double*) R_alloc(
            ((2 * m + 4) * xSize + 11 * m * m + 8 * m), sizeof(double));
    iwa = (int *) R_alloc(3 * xSize, sizeof(int));

    if (itSoftMax_ < 100)
        itSoftMax_ = 100;
    else if (itSoftMax_ > 1000)
        itSoftMax_ = 1000;

    strcpy(task, "START");

    for (int i = 0; i < xSize; ++i)
    {
        nbd[i] = 2;
    }
L111: if (iter >= itSoftMax_)
      {
          fValue_ = f;
          return 0;
      }
      if (isVerbose())
      {
          Rprintf("iter: %d itSoftMax: %d \n", iter, itSoftMax_);
      }

      Utils::setulb(xSize, m, &xBuffer_[0], &lower_[0], &upper_[0], &nbd[0], &f, &g_[0],
              factr_, &pgTol_, wa, iwa, task, tr, lsave, isave, dsave);
      iter++;
      if (strncmp(task, "FG", 2) == 0)
      {
          f = fObjective(xBuffer_);
          // Check if we have reach the threshold
          if (knowRealEnergy_)
          {
              canstop = f <= realEnergyThreshold_;
              if (canstop)
              {
                  if (isVerbose())
                  {
                      Rprintf(
                              "Have got accurate energy %.10g <= %.10g in smooth search\n",
                              f, realEnergyThreshold_);
                  }
                  fValue_ = f;
                  return 0;
              }
          }
          gradient();
          //Rprintf("LS TRACE: fValue: %.10g gradient: %.10g\n", f, g_[0]);
          goto L111;
      }
      if (strncmp(task, "NEW_X", 5) == 0)
      {
          goto L111;
      }
      // We should stop here
      fValue_ = f;
      if (fValue_ >= BIG_VALUE)
      {
          return -1;
      }
      return 0;
}

int Engine::gradient()
{
    double repsL, repsR;
    std::vector<double> x1(xBuffer_.size());
    std::vector<double> x2(xBuffer_.size());

    for (unsigned int i = 0; i < xBuffer_.size(); ++i)
    {
        std::copy(xBuffer_.begin(), xBuffer_.end(), x1.begin());
        std::copy(xBuffer_.begin(), xBuffer_.end(), x2.begin());
        repsL = reps_;
        repsR = reps_;

        x1[i] = xBuffer_[i] + repsR;
        if (x1[i] > upper_[i])
        {
            x1[i] = upper_[i];
            repsR = x1[i] - xBuffer_[i];
        }

        x2[i] = xBuffer_[i] - repsL;
        if (x2[i] < lower_[i])
        {
            x2[i] = lower_[i];
            repsL = xBuffer_[i] - x2[i];
        }
        g_[i] = (fObjective(x1) - fObjective(x2)) / (repsL + repsR);

        if (ISNA(g_[i]) || !R_FINITE(g_[i]))
        {
            g_[i] = 101.;
        }
    }
    return 0;
}

void Engine::printVect(const dVec& v)
{
    for (unsigned int i = 0; i < v.size(); ++i)
    {
        Rprintf("%.10g ", v[i]);
    }
    Rprintf("\n");
}
