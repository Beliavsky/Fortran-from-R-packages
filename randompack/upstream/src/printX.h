#ifndef PRINTX_H
#define PRINTX_H

#ifdef DEBUG_PRINT
#include <stdint.h>

void printOff(void);
void printOn(void);
int  printIsOff(void);
void printSetStatus(int noprint);
void printSetNdec(int n);
void printSetWidth(int w);
void printSetFmtChar(char f);
void printSetFmt(int w, int n, char f);
int printGetNdec(void);

void printNewl(void);                                  // print new line
void printMsg(const char *s);                          // print message
void printMsgUpper(const char *s);                     // print uppercased message
void printMsgUnderl(const char *s);                    // print underlined message

void printI(const char *name, int i);                  // print named int
void printP(const char *name, const void *p);          // print named pointer
void print64(const char *name, uint64_t x);            // print named uint64
void print32(const char *name, uint32_t x);            // print named uint32
void printD(const char *name, double d);               // print named double
void printV(const char *name, const double x[], int n);// print named vector
void printIV(const char *name, const int iv[], int n); // print named int vector
void printIVS(const char *name, const int iv[], int n); // summarized int vector
void printS(const char *name, const char *s);          // print named string

void printM(const char *name, const double A[], int nr, int nc);  // print matrix
void printMP(const char *name, const double *ap, int nr, int nc,  // print submat
             const double A[], int ldA);                          // (ap points into A)
void printMT(const char *name, const double A[], int nr, int nc); // print transpose

void print2I(const char *s, int i1, int i2);          // print 2 integers
void print3I(const char *s, int i1, int i2, int i3);  // print 3 integers
void print4I(const char *s, int i1, int i2, int i3, int i4); // print 4 integers
void print5I(const char *s, int i1, int i2, int i3, int i4, int i5); // 5 ints
void print6I(const char *s, int i1, int i2, int i3, int i4, int i5,
                           int i6); // 6 ints

#else

#define printOff()          ((void)0)
#define printOn()           ((void)0)
#define printIsOff()        (1)
#define printSetStatus(x)   ((void)0)
#define printSetNdec(x)     ((void)0)
#define printSetWidth(x)    ((void)0)
#define printSetFmtChar(x)  ((void)0)
#define printSetFmt(w,n,f)  ((void)0)
#define printGetNdec()      (0)

#define printNewl()                             ((void)0)
#define printMsg(s)                             ((void)0)
#define printMsgUpper(s)                        ((void)0)
#define printMsgUnderl(s)                       ((void)0)

#define printI(name, i)                         ((void)0)
#define printP(name, p)                         ((void)0)
#define print64(name, x)                        ((void)0)
#define print32(name, x)                        ((void)0)
#define printD(name, d)                         ((void)0)
#define printV(name, x, n)                      ((void)0)
#define printIV(name, iv, n)                    ((void)0)
#define printIVS(name, iv, n)                   ((void)0)
#define printS(name, s)                         ((void)0)

#define printM(name, A, nr, nc)                 ((void)0)
#define printMP(name, ap, nr, nc, A, ldA)       ((void)0)
#define printMT(name, A, nr, nc)                ((void)0)

#define print2I(s, i1, i2)                      ((void)0)
#define print3I(s, i1, i2, i3)                  ((void)0)
#define print4I(s, i1, i2, i3, i4)              ((void)0)
#define print5I(s, i1, i2, i3, i4, i5)          ((void)0)
#define print6I(s, i1, i2, i3, i4, i5, i6)      ((void)0)

#endif /* DEBUG_PRINT */

#endif /* PRINTX_H */
