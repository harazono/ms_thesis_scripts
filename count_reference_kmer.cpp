#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <getopt.h>
#include <math.h>
#include "kmer_library.h"
#include "cpas_debug.h"

using namespace std;

template<uint KMERSIZE>
struct FrequencyTable {
  static const size_t tablesize = ipow(5, KMERSIZE);
  typedef int Frequency;
  typedef int Score;
  vector<Frequency> kmer_table; ///< Reference side

  void printKTable(){
    for(int i = 0; i < tablesize; i++){
      fprintf(stdout, "%8d ", kmer_table[i]);
      if(i != tablesize - 1) fprintf(stdout, ", ");
    }
    fprintf(stdout, "\n\n");
  }


  public:
  FrequencyTable() : kmer_table(tablesize, 0) {}
  void countKmerFrequencies (
      const char* FASTAFileName,
      uint KmerSize,
      )
  {
    fprintf(stderr, "Loading from the FASTA file ...\r");
    const MultiFASTA multiFASTA = loadFromFASTA(FASTAFileName);
    fprintf(stderr, "Done                           \n");

    for(int i = 0; i < multiFASTA.size(); i++){
      string seqname = multiFASTA[i]->first.c_str();
      fprintf(stdout, "%s\n", seqname);
      //EXPECT_STREQ(BString2String(mf.begin()->second).c_str(), "CGACTATTCC");
    }
  }
}






int main(int argc, char *argv[]){
  GDB_On_SEGV g(argv[0]);

  struct option longopts[] = {
    // { "delete" , optional_argument , NULL , 'd' } ,
    { "kmer"      , required_argument , NULL , 'k' } ,
    { 0           , 0                 , 0    , 0  }  ,
  };

  /// PARAMETERS ///
  int kmer_size = 1;
  //////////////////
  int opt;
  int longindex;
  while ((opt = getopt_long(argc, argv, "k:", longopts, &longindex)) != -1) {
    switch (opt) {
    case 'k':
      kmer_size = atoi(optarg);
      if(kmer_size < 1 || 6 < kmer_size) {
        error("kmer size must be 1-6");
      }
      break;
    default:
      MYASSERT_NEVERREACH();
    }
  }
  const int NUM_REQUIRED_ARGUMENTS = 1;
  if(optind + NUM_REQUIRED_ARGUMENTS != argc) {
    printUsageAndExit();
  }
  const char* fasta_file_name = argv[optind + 0];
  const char* sam_file_name   = argv[optind + 1];

  #define FT() ft.countKmerFrequencies(fasta_file_name, kmer_size)
  switch(kmer_size) {
    case 1: { FrequencyTable<1> ft; FT(); } break;
    case 2: { FrequencyTable<2> ft; FT(); } break;
    case 3: { FrequencyTable<3> ft; FT(); } break;
    case 4: { FrequencyTable<4> ft; FT(); } break;
    case 5: { FrequencyTable<5> ft; FT(); } break;
    case 6: { FrequencyTable<6> ft; FT(); } break;
    default: MYASSERT_NEVERREACH();
  }
  #undef FT
  return 0;
}

