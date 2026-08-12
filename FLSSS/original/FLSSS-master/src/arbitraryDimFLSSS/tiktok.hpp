// [[Rcpp::plugins(cpp17)]]
#include <chrono>


// timetype =
// std::chrono::hours,
// std::chrono::minutes,
// std::chrono::seconds,
// std::chrono::milliseconds,
// std::chrono::microseconds,
// std::chrono::nanoseconds
template<typename timetype>
struct tiktok
{
  std::chrono::time_point<std::chrono::steady_clock> start;
  // Return time passed since tik.
  size_t tik() { start = std::chrono::steady_clock::now(); return 0; }
  // Return time passed since tok.
  size_t tok()
  {
    return std::chrono::duration_cast<timetype> (
        std::chrono::steady_clock::now() - start).count();
  }
};










