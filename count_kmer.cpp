#include <stdio.h>
#include <iostream>
#include <vector>
#include "stackdump.h"
#include "cpas_debug.h"
#include "cpas_tsv.h"

using namespace std;

typedef unsigned char Base;


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
      ASSERT_NEVERREACH_WD(DUMP(inch + 2));
  }
}

int main(int agrc, char *argv[]){
  GDB_On_SEGV g(argv[0]);
  const Base f = char2Base('&');
  return 0;
}
