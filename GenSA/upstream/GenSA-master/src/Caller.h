// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#############################################################################

#ifndef CALLER_H
#define CALLER_H

#include "Engine.h"
#include "Tracer.h"
#include "Utils.h"

#define DEFAULT_SEED    -100377


class Caller
{
    public:
        Caller()
        {
        };
        ~Caller()
        {
        };
        void execute(SEXP x_R, SEXP lb_R, SEXP ub_R, SEXP fn_R, SEXP jc_R,
                SEXP controls_R, SEXP genSAEnvironment);
        SEXP getEnergy();
        SEXP getXMiniVector();
        SEXP getTraceMat(const char*);
        SEXP getTraceMatSize() ;
        SEXP getNbFuncCall();
        void release()
        {
        };

    protected:
        static SEXP getListElement(SEXP list, char* elementName);

    private:
        Engine engine_;
};

extern "C"
{
    SEXP createInstance(void);
    SEXP releaseInstance(SEXP R_instancePtr);
    SEXP execute(SEXP x_R, SEXP lb_R, SEXP ub_R, SEXP fn_R, SEXP jc_R,
            SEXP controls_R, SEXP genSAEnvironment, SEXP R_instancePtr);
    SEXP getREnergy(SEXP R_instancePtr);
    SEXP getRXMiniVector(SEXP R_instancePtr);
    SEXP getRTraceMatSize(SEXP R_instancePtr);
    SEXP getRTraceMat(SEXP R_instancePtr, SEXP R_str);
    SEXP getRNbFuncCall(SEXP R_instancePtr);

    SEXP errR(int err)
    {
        SEXP errorCode; // Error code
        int* errorCodePtr = 0;
        PROTECT(errorCode = NEW_INTEGER(1));
        errorCodePtr = INTEGER_POINTER(errorCode);
        errorCodePtr[0] = err;
        UNPROTECT(1);
        return (errorCode);
    }
} // end extern "C"

#endif //CALLER_H
