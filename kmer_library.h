#ifndef _KMER_LIBRARY_HEADER
#define _KMER_LIBRARY_HEADER

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cctype>
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

const Base GAP_BASE = 4;
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
      return GAP_BASE;
    default:
      MYASSERT_NEVERREACH_WD(DUMP(int(inch)));
  }
}

inline char base2Char(Base b)
{
  MYASSERT_WMD("Base must be 0 to 4", 0 <= b && b < 5, DUMP(b));
  return "ACGT-"[b];
}

inline std::string BString2String(const BString& b)
{
  std::string retval;
  retval.resize(b.size());
  for(uint i = 0; i < b.size(); i++) {
    retval[i] = base2Char(b[i]);
  }
  return retval;
}

inline BString String2BString(const std::string& s)
{
  BString retval;
  retval.resize(s.size());
  for(uint i = 0; i < s.size(); i++) {
    retval[i] = char2Base(s[i]);
  }
  return retval;
}



std::string splitBy1stSpace(const std::string s) {
  std::string retval;
  retval.resize(s.size());
  for (int i = 0; i < s.size(); i++) {
    if (s[i] != ' ') {
        retval[i] = s[i];
      }else{
        retval[i] = '\0';
        return retval;
        break;
      }
  }
  return retval;
}




typedef std::string SequenceName;
typedef std::map<SequenceName, BString> MultiFASTA;

inline MultiFASTA loadFromFASTA(const std::string& inputFASTAFileName)
{
  MultiFASTA retval;
  std::ifstream ifs(inputFASTAFileName.c_str());
  if(ifs.fail()) {
    std::cerr << "ERROR: Cannot open file '" << inputFASTAFileName << "'\n";
    exit(2);
  }
  std::string tmp;
  std::string currentSequenceName;
  BString currentBString;
  while(std::getline(ifs, tmp)) {
    if(tmp.empty()) continue;
    if(tmp[0] == '>') {
      if(!currentSequenceName.empty()) {
        retval[currentSequenceName] = currentBString;
      }
      currentSequenceName = splitBy1stSpace(tmp.substr(1));
      currentBString.resize(0);
    } else {
      for(uint i = 0; i < tmp.size(); i++) {
        currentBString.push_back(char2Base(tmp[i]));
      }
    }
  }
  if(!currentSequenceName.empty()) {
    retval[currentSequenceName] = currentBString;
  }
  return retval;
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
      retval += base2Char(v % NUM_CHARS_FOR_BASES);
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

struct CIGAROp {
  char op;
  int len;
  inline CIGAROp() {}
  inline CIGAROp(char op, int len) : op(op), len(len) {}
};

typedef std::vector<CIGAROp> CIGAROPS;

inline std::string CIGAROps2String(const CIGAROPS& cops) {
  std::string retval;
  char buf[8];
  for(auto& x: cops) {
    std::sprintf(buf, "%d%c", x.len, x.op);
    retval += buf;
  }
  return retval;
}

inline CIGAROPS parseCIGARString(const std::string& cigarString) {
  CIGAROPS retval;
  const std::string& s = cigarString;
  for(size_t i = 0; i < s.size(); ++i) {
    // Parse digits
    int num = 0;
    while(std::isdigit(s[i])) {
      num = num * 10 + s[i] - '0';
      i++;
    }
    if(num <= 0) {
      std::cerr << "CIGAR op size must be positive (num = " << num << ", cigar = '" << s << "'" << std::endl;
      std::exit(2);
    }
    char op = s[i];
    switch(op) {
    case 'I':
    case 'D':
    case 'M':
    case 'X':
    case '=':
    case 'H':
    case 'S':
    case 'P':
    case 'N':
      break;
    default:
      std::cerr << "CIGAR type '" << op << "' is not valid." << std::endl;
      exit(2);
    }
    retval.push_back(CIGAROp(op, num));
  }
  return retval;
}

inline void generateAlignmentSequencesFromCIGARAndSeqs(
  const BString& refbs,
  const BString& querybs,
  const CIGAROPS& ops,
  int refStartPos,
  int queryStartPos,
  BString& ras,
  BString& qas
) {
  size_t alignmentSize = 0;
  for(auto& x: ops) {
    switch(x.op) {
    case 'M':
    case 'I':
    case 'D':
    case 'X':
    case '=':
    case 'P':
      alignmentSize += x.len;
    }
  }
  ras.resize(alignmentSize);
  qas.resize(alignmentSize);
  size_t ai = 0;
  size_t rp = refStartPos;
  size_t qp = queryStartPos;
  for(auto& x: ops) {
    switch(x.op) {
    case 'M':
    case '=':
    case 'X':
      for(int i = 0; i < x.len; ++i) {
        ras[ai] = refbs[rp++];
        qas[ai] = querybs[qp++];
        ai++;
      }
      break;
    case 'I':
      for(int i = 0; i < x.len; ++i) {
        ras[ai] = GAP_BASE;
        qas[ai] = querybs[qp++];
        ai++;
      }
      break;
    case 'D':
      for(int i = 0; i < x.len; ++i) {
        ras[ai] = refbs[rp++];
        qas[ai] = GAP_BASE;
        ai++;
      }
      break;
    case 'S':
      qp += x.len;
      break;
    case 'H':
      break;
    case 'P':
    case 'N':
      // not implemented
      std::cout << "ERROR: P/N in CIGAR are not supported." << std::endl;
      std::exit(2);
    default:
      std::cout << "ERROR: Logic error." << std::endl;
      std::exit(2);
    }
  }
}


struct SAMRecord {
  std::string qname;
  int flag;
  std::string rname;
  int pos;   ///< 0-origin position
  int mapQ;
  std::string cigar;
  std::string rnext;
  int pnext;
  int tlen;
  std::string seq;
  std::string qual;

  /// fill SAMRecord from the current line
  /// @return true when the parsing succeeds
  bool fill(FastTSVParse& ftp) {
    // See: https://samtools.github.io/hts-specs/SAMv1.pdf for the detailed specification for SAM.
    // Col Field Type
    // 1 QNAME String
    // 2 FLAG Int
    // 3 RNAME String
    // 4 POS Int
    // 5 MAPQ Int
    // 6 CIGAR String
    // 7 RNEXT String
    // 8 PNEXT Int
    // 9 TLEN Int
    // 10 SEQ String
    // 11 QUAL String
    if(ftp.size() < 11) return false;
    qname = ftp.getString(0);
    flag  = ftp.getInteger(1);
    rname = ftp.getString(2);
    pos   = ftp.getInteger(3) - 1;
    mapQ  = ftp.getInteger(4);
    cigar = ftp.getString(5);
    rnext = ftp.getString(6);
    pnext = ftp.getInteger(7);
    tlen  = ftp.getInteger(8);
    seq   = ftp.getString(9);
    qual  = ftp.getString(10);
    return true;
  }
};


inline Base Base2CompBase(Base inb)
{
  switch(inb) {
    case 0:     //when got 'A'
      return 3; //return   'T'

    case 1:     //when got 'C'
      return 2; //retern   'G'

    case 2:     //when got 'G'
      return 1; //return   'C'

    case 3:     //when got 'T'
      return 0; //return   'A'
    //case '-':
      //return GAP_BASE;
    default:
      MYASSERT_NEVERREACH_WD(DUMP(int(inb)));
  }
}


inline BString revCompBString(BString b)
{
  BString retval;
  size_t stringLength = b.size();
  retval.resize(stringLength);
  for(size_t idx = 0; idx < stringLength; idx++){
    retval[idx] = Base2CompBase(b[idx]);
  }
  for(size_t i = 0; i < stringLength / 2; i++) {
    std::swap(retval[i], retval[stringLength - 1 - i]);
  }
  return retval;
}


#endif // #ifndef _KMER_LIBRARY_HEADER
