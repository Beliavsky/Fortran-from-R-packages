/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <streambuf>
#include <ostream>
#include <cstring>
#include "soplex/spxout.h"
#include "soplex/exceptions.h"
#include "soplex/spxalloc.h"

extern "C" {
extern void Rprintf(const char *, ...);
extern void REprintf(const char *, ...);
}

namespace
{

/// A streambuf that forwards output to Rprintf (stdout) or REprintf (stderr).
class RStreambuf : public std::streambuf
{
public:
   /// If \p useStderr is true, output goes to REprintf; otherwise Rprintf.
   explicit RStreambuf(bool useStderr) : m_useStderr(useStderr) {}

protected:
   int overflow(int c) override
   {
      if(c != EOF)
      {
         char ch = static_cast<char>(c);
         if(m_useStderr)
            REprintf("%c", ch);
         else
            Rprintf("%c", ch);
      }
      return c;
   }

   std::streamsize xsputn(const char* s, std::streamsize n) override
   {
      // Write in chunks to avoid issues with non-null-terminated data
      char buf[1024];
      std::streamsize remaining = n;
      while(remaining > 0)
      {
         std::streamsize chunk = (remaining < 1023) ? remaining : 1023;
         std::memcpy(buf, s, static_cast<size_t>(chunk));
         buf[chunk] = '\0';
         if(m_useStderr)
            REprintf("%s", buf);
         else
            Rprintf("%s", buf);
         s += chunk;
         remaining -= chunk;
      }
      return n;
   }

private:
   bool m_useStderr;
};

static RStreambuf s_errBuf(true);
static RStreambuf s_outBuf(false);
static std::ostream s_errStream(&s_errBuf);
static std::ostream s_outStream(&s_outBuf);

} // anonymous namespace

namespace soplex
{
/// constructor
SPxOut::SPxOut()
   : m_verbosity(VERB_ERROR)
   , m_streams(nullptr)
{
   spx_alloc(m_streams, VERB_INFO3 + 1);
   m_streams = new(m_streams) std::ostream*[VERB_INFO3 + 1];
   m_streams[ VERB_ERROR ] = m_streams[ VERB_WARNING ] = &s_errStream;

   for(int i = VERB_DEBUG; i <= VERB_INFO3; ++i)
      m_streams[ i ] = &s_outStream;
}

//---------------------------------------------------

// destructor
SPxOut::~SPxOut()
{
   spx_free(m_streams);
}

SPxOut& SPxOut::operator=(const SPxOut& base)
{
   if(this != &base)
      m_verbosity = base.m_verbosity;

   for(int i = VERB_DEBUG; i <= VERB_INFO3; ++i)
      m_streams[ i ] = base.m_streams[ i ];

   return *this;
}

SPxOut::SPxOut(const SPxOut& rhs)
{
   m_verbosity = rhs.m_verbosity;
   m_streams = nullptr;
   spx_alloc(m_streams, VERB_INFO3 + 1);
   m_streams = new(m_streams) std::ostream*[VERB_INFO3 + 1];
   m_streams[ VERB_ERROR ] = m_streams[ VERB_WARNING ] = rhs.m_streams[VERB_ERROR];

   for(int i = VERB_DEBUG; i <= VERB_INFO3; ++i)
      m_streams[ i ] = rhs.m_streams[ i ];
}

} // namespace soplex
