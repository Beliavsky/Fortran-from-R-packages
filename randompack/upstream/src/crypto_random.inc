// -*- mode: c -*-
// Cross-platform cryptographically secure system random number generation
//
// Platform support:
// - Windows: RtlGenRandom (uses AES-CTR)
// - macOS/iOS: arc4random_buf (uses ChaCha20)
// - OpenBSD: arc4random_buf (uses ChaCha20)
// - Linux/FreeBSD/DragonFly: getrandom() syscall (Linux uses ChaCha20)
// - Fallback: /dev/urandom for other Unix-like systems
// - Final fallback: An empty function that returns false (should not happen).

#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE
#endif

#include <stdbool.h>
#include <stddef.h>

#if (defined(__linux__) || defined(__FreeBSD__) || defined(__DragonFly__)) \
    && defined(SYS_getrandom)
  #define USE_GETRANDOM 1
#endif

#if defined(_WIN32) // Windows
  #include <windows.h>
  #define RtlGenRandom SystemFunction036
  BOOLEAN NTAPI RtlGenRandom(PVOID RandomBuffer, ULONG RandomBufferLength);
  #if defined(_MSC_VER)
  #pragma comment(lib, "advapi32.lib")
  #endif
#elif defined(__APPLE__) || defined(__OpenBSD__) // macOS, iOS, OpenBSD
  #include <stdlib.h>
#elif defined(USE_GETRANDOM)  // Linux, FreeBSD, DragonFly with getrandom support
  #include <fcntl.h>
  #include <unistd.h>
  #include <errno.h>
  #include <sys/syscall.h>
#elif defined(__unix__)  // Other Unixes
  #include <fcntl.h>
  #include <unistd.h>
  #include <errno.h>
#else
  #error "Unsupported platform"
#endif

// Fill buffer with cryptographically secure random bits. True = success, false = failure
static bool crypto_random(void *buffer, size_t size) {
  if (!buffer || size == 0) return false;

#if defined(_WIN32)
  return RtlGenRandom(buffer, (ULONG)size) != 0;

#elif defined(__APPLE__) || defined(__OpenBSD__)
  arc4random_buf(buffer, size);  // never fails
  return true;

#elif defined(USE_GETRANDOM)
  ssize_t result = syscall(SYS_getrandom, buffer, size, 0);
  if (result == size) return true;  // fallback to /dev/urandom
  
#endif

#ifdef __unix__  // Fallback or other unixes
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) return false;

  unsigned char *buf = buffer;
  size_t remaining = size;
  
  while (remaining > 0) {
    ssize_t n = read(fd, buf, remaining);
    if (n <= 0) {
      if (errno == EINTR) continue;
      close(fd);
      return false;
    }
    buf += n;
    remaining -= n;
  }
  
  close(fd);
  return true;
#else
  return false;
#endif
}
