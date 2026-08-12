/*
 * Custom std::ostream replacements for std::cout and std::cerr.
 * These write through R's Rprintf/REprintf so that R CMD check
 * does not flag use of std::cout/std::cerr.
 *
 * Usage: replace std::cerr with r_cerr() and std::cout with r_cout().
 * The returned ostream& supports all << operators identically.
 */

#include <ostream>
#include <streambuf>
#include <cstring>
#include <algorithm>

extern "C" {
extern void Rprintf(const char *, ...);
extern void REprintf(const char *, ...);
}

class RStreamBuf : public std::streambuf {
    void (*m_fn)(const char *, ...);
protected:
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        /* Write in chunks to avoid stack overflow on large output */
        const std::streamsize max_chunk = 511;
        std::streamsize written = 0;
        char tmp[512];
        while (written < n) {
            std::streamsize chunk = std::min(n - written, max_chunk);
            std::memcpy(tmp, s + written, static_cast<size_t>(chunk));
            tmp[chunk] = '\0';
            m_fn("%s", tmp);
            written += chunk;
        }
        return n;
    }

    int overflow(int c) override {
        if (c != EOF) {
            char ch = static_cast<char>(c);
            m_fn("%c", ch);
        }
        return c;
    }

public:
    explicit RStreamBuf(void (*fn)(const char *, ...)) : m_fn(fn) {}
};

static RStreamBuf s_cout_buf(Rprintf);
static RStreamBuf s_cerr_buf(REprintf);
static std::ostream s_r_cout(&s_cout_buf);
static std::ostream s_r_cerr(&s_cerr_buf);

std::ostream& r_cout() { return s_r_cout; }
std::ostream& r_cerr() { return s_r_cerr; }
