// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#############################################################################

#ifndef ENGINE_H_
#define ENGINE_H_

#include <vector>
#include <math.h>
#include "Tracer.h"
#include "Utils.h"

#define BIG_VALUE	1.e13

class Engine
{
    private:
        double etot_;
        double etot0_;
        double eMini_;
        double fValue_;
        double temSta_;
        double tem_ ;
        double qv_;
        double qa_;
        double reps_;
        double pgTol_;
        double factr_;
        double realEnergyThreshold_;
        double maxTime_;
        double timeSpan_;
        double temRestart_;
        clock_t startTime_;
        clock_t endTime_;
        bool hasConstraint_;
        long int idum_;
        int markovLength_;
        int nbFctCall_;
        int maxFctCall_;
        int method_;
        int itSoftMax_;
        int indTrace_;
        int maxStep_;
        int interval_;
        int noImprovementStop_;
        bool knowRealEnergy_;
        bool lsEnd_ ;
        bool isSimpleFuction_;
        bool useTraceMat_;
        dVec x_;
        dVec xMini_;
        dVec xBuffer_;
        dVec xBackup_;
        dVec lower_;
        dVec upper_;
        dVec xRange_;
        dVec g_;
        Tracer tracer_;

    public:
        Engine();
        ~Engine();
        OptStruct rEnv_;

        enum
        {
            SMOOTH, HARD
        };

        bool isVerbose()
        {
            return false ;
        }


        bool isUserVerbose()
        {
            OptStruct OS = (OptStruct) rEnv_;
            return (bool)(OS->verbose);
        }

        void setIsSimpleFunction(bool b)
        {
            isSimpleFuction_ = b;
        }

        bool getIsSimpleFunction()
        {
            return isSimpleFuction_;
        }


        void setTemSta(double v)
        {
            temSta_ = v;
        }
        double getTemSta()
        {
            return temSta_;
        }

        void setMarkovLength(int v)
        {
            markovLength_ = v;
        }
        int getMarkovLenght()
        {
            return markovLength_;
        }

        void setMaxStep(int v)
        {
            maxStep_ = v;
        }
        int getMaxStep()
        {
            return maxStep_;
        }

        void setInterval(int v)
        {
            interval_ = v;
        }
        int getInterval()
        {
            return interval_;
        }

        void setKnowRealEnergy(bool v)
        {
            knowRealEnergy_ = v;
        }
        bool knowRealEnergy()
        {
            return knowRealEnergy_;
        }

        void setRealEnergyThreshold(double v)
        {
            realEnergyThreshold_ = v;
        }
        double getRealEnergyThreshold()
        {
            return realEnergyThreshold_;
        }

        void setMaxTime(double v)
        {
            maxTime_ = v;
        }
        double getMaxTime()
        {
            return maxTime_;
        }

        void setNoImprovementStop(int v)
        {
            noImprovementStop_ = v;
            if (!useTraceMat_ && v != -1) {
                Rf_warning("Warning, trace mat has to be enabled to make no improvement stop criteria working");
            }
        }
        int getNoImprovementStop()
        {
            return noImprovementStop_;
        }

        void setMethod(int v)
        {
            method_ = v;
        }
        int getMethod()
        {
            return method_;
        }

        void setMaxFctCall(int v)
        {
            maxFctCall_ = v;
        }
        int getMaxFctCall()
        {
            return maxFctCall_;
        }

        void setHasConstraint(bool v)
        {
            hasConstraint_ = v;
        }
        bool hasConstraint()
        {
            return hasConstraint_;
        }

        void setLower(int size, double* data)
        {
            lower_.assign(data, data+size) ;
        }
        const dVec& getLower()
        {
            return lower_;
        }

        void setUpper(int size, double* data)
        {
            upper_.assign(data, data+size) ;
        }
        const dVec& getUpper()
        {
            return upper_;
        }

        void setX(int size, double* data)
        {
            x_.assign(data, data+size) ;
        }
        const dVec& getX()
        {
            return x_;
        }

        void setInitialTemp(double v)
        {
            temSta_ = v;
        }
        double getInitialTemp()
        {
            return temSta_;
        }

        void setVisitingParam(double v)
        {
            qv_ = v;
        }
        double getVisitingParam()
        {
            return qv_;
        }

        void setAcceptanceParam(double v)
        {
            qa_ = v;
        }
        double getAcceptanceParam()
        {
            return qa_;
        }

        double getEmini()
        {
            return eMini_;
        }

        const dVec& getXMini()
        {
            return xMini_;
        }

        long int getNbFctCall()
        {
            return nbFctCall_;
        }

        void setREnv(const OptStruct env)
        {
            rEnv_ = env;
        }

        const Tracer& getTracer()
        {
            return tracer_;
        }

        void setLSEnd(bool b)
        {
            lsEnd_ = b ;
        }

        void setTemRestart(double v)
        {
            temRestart_ = v;
        }

        void setUseTraceMat(bool b)
        {
            useTraceMat_ = b;
            if (!b && noImprovementStop_ > 0)
            {
                Rf_warning("Trace mat has to be activated to have no improvment stop criteria to work.");
            }
        }

        bool useTraceMat()
        {
            return useTraceMat_;
        }
        
        void setSeed(long int seed)
        {
            idum_ = seed;
        }

        int initialize();
        int startSearch();
        bool checkStoping();
        void stopSearch();
        void coordin(long int, dVec&);
        int energy(const dVec&);
        double lsEnergy(dVec&);
        double visita(double*, double*, long int*);
        double fObjective(const dVec& x);
        double fn(const dVec& x);
        bool judgeConstraint();
        int smoothSearch();
        int hardSearch();
        int gradient();
        void printVect(const dVec&);

};

#endif
