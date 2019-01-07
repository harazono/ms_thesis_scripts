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

void error(const char* msg) {
  fprintf(stderr, "ERROR: %s\n", msg);
  exit(2);
}


template<uint KMERSIZE>
struct FrequencyTable {
  static const size_t tablesize = ipow(4, KMERSIZE);
  typedef int Frequency;
  typedef int Score;
  vector<Frequency> kmer_table; ///< Reference side

  void printKTable(){
    for(int i = 0; i < tablesize; i++){
      fprintf(stdout, "%d\n", kmer_table[i]);
      if(i != tablesize - 1) fprintf(stdout, ", ");
    }
    fprintf(stdout, "\n\n");
  }


  int kmer2index(string str){
    int retval = 0;
    int k = str.size();
    for(int i = 0; i < k; i++){
    retval = (retval * 5 + char2Base(str[i])) % tablesize;
    }
    return retval;
  }


  public:
  FrequencyTable() : kmer_table(tablesize, 0) {}
  void countKmerFrequencies (
      const char* FASTAFileName,
      uint KmerSize
      )
  {
    fprintf(stderr, "Loading from the FASTA file ...\r");
    const MultiFASTA multiFASTA = loadFromFASTA(FASTAFileName);
    fprintf(stderr, "Done                           \n");

    const int kmer_size = KmerSize;
    for(auto itr = multiFASTA.begin(); itr != multiFASTA.end(); ++itr) {
      cerr << "begin to proceeding " << itr->first << "\"" << endl;
      cerr << "size of " << itr->first << ": " << itr->second.size() << endl;

      int chrlen = itr->second.size();
      string refseq = BString2String(itr->second).c_str();
      for(int i = 0; i < chrlen - KmerSize; i++){
        string kmerstr = refseq.substr(i, kmer_size);
        int idx = kmer2index(kmerstr);
        kmer_table[idx] += 1;
      }
    }// end of for(auto itr = multiFASTA.begin(); itr != multiFASTA.end(); ++itr)
  printKTable();
  }
};






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
    fprintf(stderr, "few args\n");
  }
  const char* fasta_file_name = argv[optind + 0];

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

