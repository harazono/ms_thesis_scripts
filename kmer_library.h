#ifndef _KMER_LIBRARY_HEADER
#define _KMER_LIBRARY_HEADER

#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "stackdump.h"
#include "cpas_debug.h"
#include "cpas_tsv.h"

typedef unsigned long ulong;
typedef unsigned int uint;
typedef unsigned char Base;
typedef std::vector<Base> BString; ///< string of bases

const uint NUM_CHARS_FOR_BASES = 5;

constexpr ulong ipow(ulong i, uint N)
{
  return N <= 1 ? i : i * ipow(i, N - 1);
}

inline Base char2Base(char inch)
{
  switch(inch) {
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    case '-':
      return 4;
    default:
      MYASSERT_NEVERREACH_WD(DUMP(inch));
  }
}

inline char base2Char(Base b)
{
  MYASSERT_WMD("Base must be 0 to 4", 0 <= b && b < 5, DUMP(b));
  return "ACGT-"[b];
}

/**
 * k-mer integer class
 */
template<uint K>
class KInt {
  uint64_t kint;
  static const size_t KMERINT_COUNT = ipow(NUM_CHARS_FOR_BASES, K);

public:
  /// Put a new base at LSB, and lift the other stuff toward MSB.
  /// The overflown digit is thrown away.
  inline void ShiftIn(Base b) {
    kint = (kint * NUM_CHARS_FOR_BASES + b) % KMERINT_COUNT;
  }

  inline KInt() : kint(0) {}
  inline KInt(uint64_t v) : kint(v) {}
  /// Construct by string
  /// e.g.) KInt k("ACG-");
  inline KInt(const char* s) : kint(0) {
    MYASSERT_WMD("KInt must be constructed with a string of size K", std::strlen(s) == K, DUMP(s));
    for(ulong i = 0; i < K; i++) {
      ShiftIn(char2Base(s[i]));
    }
  }

  /// Output string for debug
  inline std::string str() const {
    std::string retval;
    uint64_t v = kint;
    for(ulong i = 0; i < K; i++) {
      retval += char2Base(v % NUM_CHARS_FOR_BASES);
      v /= NUM_CHARS_FOR_BASES;
    }
    for(size_t i = 0; i < K / 2; i++) {
      std::swap(retval[i], retval[K - 1 - i]);
    }
    return retval;
  }
  inline operator uint64_t () const {
    return kint;
  }
};

#endif // #ifndef _KMER_LIBRARY_HEADER
