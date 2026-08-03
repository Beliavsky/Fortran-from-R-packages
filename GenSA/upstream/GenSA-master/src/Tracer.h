// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################

#ifndef TRACER_H_
#define TRACER_H_

#include "Utils.h"
#define INIT_TRACE_SIZE 1000

class Tracer
{
    public:
        Tracer() {};
        virtual ~Tracer() {};
        void setKeyList(const strVec& keylist) ;
        const double* getVectorPtr(const std::string& key) ;
        double getLastValue(const std::string& key) ;
        void updateLastValue(const std::string& key, double value);
        void addValue(const std::string& key, double value) ;
        unsigned int getTracerLength() ;
        void clear() ;
    private:
        Map traceMap_ ;
};

#endif
