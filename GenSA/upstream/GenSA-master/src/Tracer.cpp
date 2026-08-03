// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################

#include "Tracer.h"
#include <limits>

void Tracer::setKeyList(const strVec& keyList)
{
    strVecIt it;

    for (it = keyList.begin(); it != keyList.end(); it++)
    {
        traceMap_.insert(std::make_pair(*it, dVec()));
    }
}

const double* Tracer::getVectorPtr(const std::string& key)
{
    MapIt it = traceMap_.find(key);

    if (it != traceMap_.end())
    {
        return (const double*) (&(it->second[0]));
    }
    else
    {
        return 0;
    }

}

double Tracer::getLastValue(const std::string& key)
{
    MapIt it = traceMap_.find(key);
    if (it != traceMap_.end())
    {
        const dVec& vec = it->second;
        return vec[vec.size() - 1];
    }
    else
    {
        return DBL_MIN;
    }
}

void Tracer::updateLastValue(const std::string& key, double value)
{
    MapIt it = traceMap_.find(key);
    if (it != traceMap_.end())
    {
        dVec& vec = it->second;
        vec[vec.size() - 1] = value;
    }
}

void Tracer::addValue(const std::string& key, double value)
{
    MapIt it = traceMap_.find(key);
    if (it != traceMap_.end())
    {
        it->second.push_back(value);
    }
}

void Tracer::clear()
{
    MapIt it;
    for (it = traceMap_.begin(); it != traceMap_.end(); it++)
    {
        it->second.clear();
    }
}

unsigned int Tracer::getTracerLength()
{
    unsigned int size = 0;
    MapIt it = traceMap_.begin();
    size = it->second.size();
    it++;
    for (; it != traceMap_.end(); it++)
    {
        if (size != it->second.size())
        {
            return 0;
        }
        else
        {
            size = it->second.size();
        }
    }
    return size;
}
