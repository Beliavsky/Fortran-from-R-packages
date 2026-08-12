// [[Rcpp::plugins(cpp11)]]
# include <Rcpp.h>
# include <bitset>
using namespace Rcpp;
# define vec std::vector


// #ifndef uint64_t
// #define uint64_t std::uint64_t
// #endif


constexpr uint64_t ENDBIT = ((uint64_t)1) << (sizeof(uint64_t) * 8 - 1);
constexpr uint64_t _ENDBIT = ~ENDBIT;
constexpr uint64_t Nfull = sizeof(uint64_t) * 8;
constexpr uint64_t Nfull_1 = Nfull - 1;
constexpr uint64_t Nhalf = Nfull / 2;
constexpr uint64_t lowHalf = ((uint64_t)0 - 1) >> Nhalf;
constexpr uint64_t uppHalf = ((uint64_t)0 - 1) << Nhalf;


// z = x + y, x or/and y can be z.
void addHint(uint64_t *z, uint64_t *x, uint64_t *y, int size)
{
  z[0] = x[0] + y[0];
  for(int i = 1; i < size; ++i)
  {
    uint64_t carry = z[i - 1] >> Nfull_1;
    z[i - 1] &= _ENDBIT;
    z[i] = x[i] + y[i] + carry;
  }
}



// [[Rcpp::export]]
List testAddHintV(IntegerVector xv, IntegerVector yv)
{
  vec<int> x(xv.begin(), xv.end()), y(yv.begin(), yv.end());
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());
  int size = x.size() / 64;
  vec<uint64_t> a(size, 0), b(size, 0);
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      a[i] += ((uint64_t)x[k]) << j;
  }
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      b[i] += ((uint64_t)y[k]) << j;
  }


  if(false)
  {
    std::string astr = "";
    for(int i = a.size() - 1; i >= 0; --i)
      astr += std::bitset<64>(a[i]).to_string();
    std::string bstr = "";
    for(int i = b.size() - 1; i >= 0; --i)
      bstr += std::bitset<64>(b[i]).to_string();
    std::cout << "astr = " << astr << std::endl;
    std::cout << "bstr = " << bstr << std::endl;
  }


  vec<uint64_t> ab(a.size(), 0);
  addHint(&ab[0], &a[0], &b[0], a.size());
  // mulHint(&ab[0], ab.size(), &a[0], a.size(), &b[0], b.size());
  std::string abstr = "";
  for(int i = ab.size() - 1; i >= 0; --i)
    abstr += std::bitset<64>(ab[i]).to_string();


  vec<int> xnew(x.size(), 0), ynew(y.size(), 0);
  // Skip every 64th bit.
  for(int i = 0, iend = x.size() / 64, *xnewk = &xnew[0]; i < iend; ++i)
  {
    int *lb = &x[0] + i * 64;
    int *ub = lb + 63;
    std::copy(lb, ub, xnewk);
    xnewk += 63;
  }
  for(int i = 0, iend = y.size() / 64, *ynewk = &ynew[0]; i < iend; ++i)
  {
    int *lb = &y[0] + i * 64;
    int *ub = lb + 63;
    std::copy(lb, ub, ynewk);
    ynewk += 63;
  }
  x.swap(xnew); y.swap(ynew);


  vec<int> added(x.size(), 0);
  for(int i = 0, iend = x.size(); i < iend; ++i) added[i] = x[i] + y[i];


  for(int i = 0, iend = added.size() - 1; i < iend; ++i)
  {
    added[i + 1] += added[i] / 2;
    added[i] %= 2;
  }
  vec<int> addednew(added.size() * 2, 0);
  for(int i = 0, iend = (added.size() + 62) / 63; i < iend; ++i)
  {
    int *lb = &added[0] + i * 63;
    int *ub = std::min(lb + 63, &*added.end());
    std::copy(lb, ub, addednew.begin() + i * 64);
  }
  addednew.resize(added.size());
  added.swap(addednew);
  added.back() = added[added.size()];


  std::reverse(added.begin(), added.end());


  std::string correctRst(added.size(), '0');
  for(int i = 0, iend = added.size(); i < iend; ++i)
    correctRst[i] = added[i] == 1 ? '1' : '0';


  return List::create(Named("rst") = abstr, Named("correctRst") = correctRst);
}


// z = x - y, x or/and y can be z.
void subHint(uint64_t *z, uint64_t *x, uint64_t *y, int size)
{
  z[0] = x[0] - y[0];
  for(int i = 1; i < size; ++i)
  {
    uint64_t carry = z[i - 1] >> Nfull_1;
    z[i - 1] &= _ENDBIT;
    z[i] = x[i] - carry - y[i];
  }
}


// [[Rcpp::export]]
List testSubHintV(IntegerVector xv, IntegerVector yv)
{
  vec<int> x(xv.begin(), xv.end()), y(yv.begin(), yv.end());
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());
  int size = x.size() / 64;
  vec<uint64_t> a(size, 0), b(size, 0);
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      a[i] += ((uint64_t)x[k]) << j;
  }
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      b[i] += ((uint64_t)y[k]) << j;
  }


  if(false)
  {
    std::string astr = "";
    for(int i = a.size() - 1; i >= 0; --i)
      astr += std::bitset<64>(a[i]).to_string();
    std::string bstr = "";
    for(int i = b.size() - 1; i >= 0; --i)
      bstr += std::bitset<64>(b[i]).to_string();
    std::cout << "astr = " << astr << std::endl;
    std::cout << "bstr = " << bstr << std::endl;
  }


  vec<uint64_t> ab(a.size(), 0);
  subHint(&ab[0], &a[0], &b[0], a.size());
  // mulHint(&ab[0], ab.size(), &a[0], a.size(), &b[0], b.size());
  std::string abstr = "";
  for(int i = ab.size() - 1; i >= 0; --i)
    abstr += std::bitset<64>(ab[i]).to_string();


  vec<int> xnew(x.size(), 0), ynew(y.size(), 0);
  // Skip every 64th bit.
  for(int i = 0, iend = x.size() / 64, *xnewk = &xnew[0]; i < iend; ++i)
  {
    int *lb = &x[0] + i * 64;
    int *ub = lb + 63;
    std::copy(lb, ub, xnewk);
    xnewk += 63;
  }
  for(int i = 0, iend = y.size() / 64, *ynewk = &ynew[0]; i < iend; ++i)
  {
    int *lb = &y[0] + i * 64;
    int *ub = lb + 63;
    std::copy(lb, ub, ynewk);
    ynewk += 63;
  }
  x.swap(xnew); y.swap(ynew);


  vec<int> subed(x.size(), 0);
  for(int i = 0, iend = x.size(); i < iend; ++i) subed[i] = x[i] - y[i];


  for(int i = 0, iend = subed.size() - 1; i < iend; ++i)
  {
    // subed[i + 1] += subed[i] / 2;
    subed[i + 1] += (subed[i] - 1) / 2;
    subed[i] = std::abs(subed[i] % 2);
  }
  vec<int> subednew(subed.size() * 2, 0);
  for(int i = 0, iend = (subed.size() + 62) / 63; i < iend; ++i)
  {
    int *lb = &subed[0] + i * 63;
    int *ub = std::min(lb + 63, &*subed.end());
    std::copy(lb, ub, subednew.begin() + i * 64);
  }
  subednew.resize(subed.size());
  subed.swap(subednew);
  subed.back() = subed[subed.size()];


  std::reverse(subed.begin(), subed.end());


  std::string correctRst(subed.size(), '0');
  for(int i = 0, iend = subed.size(); i < iend; ++i)
    correctRst[i] = subed[i] == 1 ? '1' : '0';


  return List::create(Named("rst") = abstr, Named("correctRst") = correctRst);
}


// [[Rcpp::export]]
void seeWhatInsideNegativeMul(int x, int y)
{
  int z = x * y;
  std::cout << "x = " << std::bitset<32>(x).to_string() << std::endl;
  std::cout << "y = " << std::bitset<32>(y).to_string() << std::endl;
  std::cout << "z = " << std::bitset<32>(z).to_string() << std::endl;
}


// [[Rcpp::export]]
void seeWhatInsideMinus(int x, int y)
{
  int z = x - y;
  std::cout << "x = " << std::bitset<32>(x).to_string() << std::endl;
  std::cout << "y = " << std::bitset<32>(y).to_string() << std::endl;
  std::cout << "z = " << std::bitset<32>(z).to_string() << std::endl;
}


// Store x * y in 2 containers L (lower half) and U (upper half).
void mulHint(const uint64_t &x, const uint64_t &y, uint64_t &L, uint64_t &U)
{
  uint64_t xlow = x & lowHalf;
  uint64_t xhih = x >> Nhalf;
  uint64_t ylow = y & lowHalf;
  uint64_t yhih = y >> Nhalf;
  uint64_t xlow_ylow = xlow * ylow;
  uint64_t xhih_yhih = xhih * yhih;
  uint64_t xlow_yhih = xlow * yhih;
  uint64_t xhih_ylow = xhih * ylow;
  uint64_t mid = (xlow_yhih & lowHalf) + (xhih_ylow & lowHalf) +
    (xlow_ylow >> Nhalf);
  uint64_t carry = mid >> Nhalf;
  L = (mid << Nhalf) + (xlow_ylow & lowHalf);
  U = (xlow_yhih >> Nhalf) + (xhih_ylow >> Nhalf) + xhih_yhih + carry;
  U = (U << 1) | (L >> Nfull_1);
  L &= _ENDBIT;
}


// x and y are boolean values.
// [[Rcpp::export]]
List testMulHint(IntegerVector x, IntegerVector y)
{
  uint64_t a = 0;
  for(int k = 0, i = std::min<int> (64, x.size()) - 1; i >= 0; --i, ++k)
    a += (uint64_t)x[i] << k;
  uint64_t b = 0;
  for(int k = 0, i = std::min<int> (64, x.size()) - 1; i >= 0; --i, ++k)
    b += (uint64_t)y[i] << k;
  // std::cout << std::bitset<sizeof(uint64_t) * 8> (a) << ", " <<
  //   std::bitset<sizeof(uint64_t) * 8> (b) << std::endl;
  uint64_t L, U;
  mulHint(a, b, L, U);
  std::bitset<sizeof(uint64_t) * 8> Lb(L), Ub(U);
  std::string rst = Ub.to_string() + Lb.to_string();
  // std::reverse(rst.begin(), rst.end());


  vec<unsigned char> xtmp(x.begin(), x.end()), ytmp(y.begin(), y.end());
  std::reverse(xtmp.begin(), xtmp.end());
  std::reverse(ytmp.begin(), ytmp.end());
  vec<unsigned char> z(x.size() + y.size(), 0);
  for(int i = 0, iend = x.size(); i < iend; ++i)
  {
    for(int j = 0, jend = y.size(); j < jend; ++j)
    {
      z[i + j] += xtmp[i] * ytmp[j];
    }
  }
  std::reverse(z.begin(), z.end());
  for(int i = z.size() - 1; i >= 1; --i)
  {
    z[i - 1] += z[i] / 2;
    z[i] %= 2;
  }
  std::copy(z.begin() + 1, z.begin() + x.size(), z.begin());
  z[x.size() - 1] = z[x.size()];
  z[x.size()] = 0;


  std::string correctrst(z.size(), '0');
  // Rcout << "correctrst.size() = " << correctrst.size() << "\n";
  for(int i = 0, iend = z.size(); i < iend; ++i)
    correctrst[i] = z[i] == 1 ? '1' : '0';


  return List::create(Named("rst") = rst, Named("correctRst") = correctrst);
}


// z = x * y, unsigned.
// zsize is at least ysize + xsize.
// z does not have to be filled with 0.
// z cannot equal x or y.
void mulHint(uint64_t *z, int zsize,
             uint64_t *x, int xsize,
             uint64_t *y, int ysize)
{
  int i, j, k;
  uint64_t L, U, carry;


  for(i = xsize - 1; i >= 0 and x[i] == 0; --i);
  xsize = i + 1;
  for(j = ysize - 1; j >= 0 and y[j] == 0; --j);
  ysize = j + 1;
  std::fill(z, z + zsize, 0);


  for(i = 0; i < xsize; ++i)
  {
    for(j = 0; j < ysize; ++j)
    {
      mulHint(x[i], y[j], L, U);
      k = i + j;
      z[k] += L;
      carry = (z[k] >> Nfull_1) + U;
      z[k] &= _ENDBIT;
      for(int u = k + 1; carry != 0; ++u)
      {
        z[u] += carry;
        carry = z[u] >> Nfull_1;
      }
    }
  }
}




// // [[Rcpp::export]]
// List testMulHintV32(IntegerVector xv, IntegerVector yv)
// {
//   vec<uint32_t> x(xv.begin(), xv.end()), y(yv.begin(), yv.end());
//   std::reverse(x.begin(), x.end());
//   std::reverse(y.begin(), y.end());
//   int size = x.size() / 32;
//   vec<uint32_t> a(size, 0), b(size, 0);
//   for(int i = 0, k = 0; i < size; ++i)
//   {
//     for(int j = 0; j < 32; ++k, ++j)
//       a[i] += ((uint32_t)x[k]) << j;
//   }
//   for(int i = 0, k = 0; i < size; ++i)
//   {
//     for(int j = 0; j < 32; ++k, ++j)
//       b[i] += ((uint32_t)y[k]) << j;
//   }
//
//
//   if(false)
//   {
//     std::string astr = "";
//     for(int i = a.size() - 1; i >= 0; --i)
//       astr += std::bitset<32>(a[i]).to_string();
//     std::string bstr = "";
//     for(int i = b.size() - 1; i >= 0; --i)
//       bstr += std::bitset<32>(b[i]).to_string();
//     std::cout << "astr = " << astr << std::endl;
//     std::cout << "bstr = " << bstr << std::endl;
//   }
//
//
//   vec<uint32_t> ab(a.size() + b.size(), 0);
//   mulHint(&ab[0], ab.size(), &a[0], a.size(), &b[0], b.size());
//   std::string abstr = "";
//   for(int i = ab.size() - 1; i >= 0; --i)
//     abstr += std::bitset<32>(ab[i]).to_string();
//
//
//   vec<int> conv(x.size() + y.size(), 0);
//   for(int i = 0, iend = x.size(); i < iend; ++i)
//   {
//     for(int j = 0, jend = y.size(); j < jend; ++j)
//       conv[i + j] += x[i] * y[j];
//   }
//   for(int i = 0, iend = conv.size() - 1; i < iend; ++i)
//   {
//     conv[i + 1] += conv[i] / 2;
//     conv[i] %= 2;
//   }
//   for(int i = 0, iend = conv.size() / 64; i < iend; ++i)
//   {
//     int k = (i + 1) * 64 - 1;
//     if(k - conv.size() < 0)
//     {
//       conv.insert(conv.begin() + k, 0);
//       conv.resize(conv.size() - 1);
//     }
//   }
//
//
//   std::reverse(conv.begin(), conv.end());
//   for(int i = 0, iend = conv.size(); i < iend; i += 64)
//   {
//     std::copy(conv.begin() + i + 1, conv.begin() + i + 64, conv.begin() + i);
//     if(i + 64 < iend) conv[i + 63] = conv[i + 64];
//   }
//
//
//   std::string correctRst(conv.size(), '0');
//   for(int i = 0, iend = conv.size(); i < iend; ++i)
//     correctRst[i] = conv[i] == 1 ? '1' : '0';
//
//
//   return List::create(Named("rst") = abstr, Named("correctRst") = correctRst);
// }


/*
// [[Rcpp::export]]
List testMulHintV(IntegerVector xv, IntegerVector yv)
{
  vec<int> x(xv.begin(), xv.end()), y(yv.begin(), yv.end());
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());
  int size = x.size() / 64;
  vec<uint64_t> a(size, 0), b(size, 0);
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      a[i] += ((uint64_t)x[k]) << j;
  }
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      b[i] += ((uint64_t)y[k]) << j;
  }


  if(false)
  {
    std::string astr = "";
    for(int i = a.size() - 1; i >= 0; --i)
      astr += std::bitset<64>(a[i]).to_string();
    std::string bstr = "";
    for(int i = b.size() - 1; i >= 0; --i)
      bstr += std::bitset<64>(b[i]).to_string();
    std::cout << "astr = " << astr << std::endl;
    std::cout << "bstr = " << bstr << std::endl;
  }


  vec<uint64_t> ab(a.size() + b.size(), 0);
  mulHint(&ab[0], ab.size(), &a[0], a.size(), &b[0], b.size());
  std::string abstr = "";
  for(int i = ab.size() - 1; i >= 0; --i)
    abstr += std::bitset<64>(ab[i]).to_string();


  if(false)
  {
    vec<int> xnew(x.size(), 0), ynew(y.size(), 0);
    // Skip every 64th bit.
    for(int i = 0, iend = x.size() / 64, *xnewk = &xnew[0]; i < iend; ++i)
    {
      int *lb = &x[0] + i * 64;
      int *ub = lb + 63;
      std::copy(lb, ub, xnewk);
      xnewk += 63;
    }
    for(int i = 0, iend = y.size() / 64, *ynewk = &ynew[0]; i < iend; ++i)
    {
      int *lb = &y[0] + i * 64;
      int *ub = lb + 63;
      std::copy(lb, ub, ynewk);
      ynewk += 63;
    }
    x.swap(xnew); y.swap(ynew);
  }


  if(false)
  {
    for(int i = 0, iend = x.size(); i < iend; ++i)
      Rcout << x[i] << ",";
    Rcout << "\n";
    for(int i = 0, iend = y.size(); i < iend; ++i)
      Rcout << y[i] << ",";
  }


  vec<int> conv(x.size() + y.size(), 0);
  for(int i = 0, iend = x.size(); i < iend; ++i)
  {
    for(int j = 0, jend = y.size(); j < jend; ++j)
    {
      conv[i + j] += x[i] * y[j];
    }
  }
  for(int i = 0, iend = conv.size() - 1; i < iend; ++i)
  {
    if((i + 1) % 64 != 0)
    {
      conv[i + 1] += conv[i] / 2;
      conv[i] %= 2;
    }
    else
    {
      conv[i + 1] += conv[i];
      conv[i] = 0;
    }
  }


  if(false)
  {
    vec<int> convnew(conv.size(), 0);
    if(false)
    {
      for(int i = 0, iend = (conv.size() + 62) / 63; i < iend; ++i)
      {
        int *lb = &conv[0] + i * 63;
        int *ub = std::min(lb + 63, &*conv.end());
        if(convnew.begin() + i * 64 >= convnew.end()) break;
        std::copy(lb, ub, convnew.begin() + i * 64);
      }
      convnew.resize(conv.size());
      conv.swap(convnew);
    }
  }


  std::reverse(conv.begin(), conv.end());


  std::string correctRst(conv.size(), '0');
  for(int i = 0, iend = conv.size(); i < iend; ++i)
    correctRst[i] = conv[i] == 1 ? '1' : '0';


  return List::create(Named("rst") = abstr, Named("correctRst") = correctRst);
}
*/


void intString2hugeInt(vec<uint64_t> &rst, char *s, int ssize)
{
  if(ssize <= 0) return;
  int Nbit = int(ssize * 3.33) + 1;
  int nbitPeruint64_t = sizeof(uint64_t) * 8 - 1;
  int rstSize = (Nbit + nbitPeruint64_t - 1) / nbitPeruint64_t;
  rst.resize(rstSize, 0);
  vec<uint64_t> buffer(rstSize * 3, 0);
  uint64_t *power10 = &buffer[0], *power10new = power10 + rstSize;
  uint64_t *adder = power10new + rstSize;
  power10[0] = 1;
  uint64_t ten = 10;
  for(int i = ssize - 1; ; --i)
  {
    uint64_t a = s[i] - 48;
    mulHint(adder, rstSize, &a, 1, power10, rstSize);
    addHint(&rst[0], &rst[0], adder, rstSize);
    if(i <= 0) break;
    mulHint(power10new, rstSize, power10, rstSize, &ten, 1);
    std::swap(power10new, power10);
  }
}


// [[Rcpp::export]]
IntegerVector testIntString2hugeInt(std::string &ss)
{
  vec<uint64_t> rst;
  char *s = &ss[0];
  intString2hugeInt(rst, s, ss.size());


  if(true)
  {
    int i = rst.size() - 1;
    for(; i >= 0 and rst[i] == 0; --i);
    rst.resize(i + 1);
  }
  IntegerVector result(rst.size() * sizeof(uint64_t) * 8);
  int Nbits = sizeof(uint64_t) * 8;
  uint64_t masks[Nbits];
  uint64_t one = 1;
  for(int i = 0; i < Nbits; ++i) masks[i] = one << i;


  for(int i = 0, iend = rst.size(), k = result.size() - 1; i < iend; ++i)
  {
    for(int u = 0, uend = Nbits - 1; u < uend; ++u)
    {
      result[k] = (rst[i] & masks[u]) != 0;
      --k;
    }
  }


  return result;
}


// [[Rcpp::export]]
IntegerVector testMulHintV(std::string &xv, std::string &yv)
{
  vec<uint64_t> x, y;
  intString2hugeInt(x, &xv[0], xv.size());
  intString2hugeInt(y, &yv[0], yv.size());
  if(true)
  {
    int i = x.size() - 1;
    for(; i >= 0 and x[i] == 0; --i);
    x.resize(i + 1);
    i = y.size() - 1;
    for(; i >= 0 and y[i] == 0; --i);
    y.resize(i + 1);
  }


  vec<uint64_t> rst(x.size() + y.size());
  mulHint(&rst[0], rst.size(), &x[0], x.size(), &y[0], y.size());


  IntegerVector result(rst.size() * sizeof(uint64_t) * 8);
  int Nbits = sizeof(uint64_t) * 8;
  uint64_t masks[Nbits];
  uint64_t one = 1;
  for(int i = 0; i < Nbits; ++i) masks[i] = one << i;


  for(int i = 0, iend = rst.size(), k = result.size() - 1; i < iend; ++i)
  {
    for(int u = 0, uend = Nbits - 1; u < uend; ++u)
    {
      result[k] = (rst[i] & masks[u]) != 0;
      --k;
    }
  }


  return result;
}























































