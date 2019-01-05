#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <getopt.h>
#include "kmer_library.h"

using namespace std;

void error(const char* msg) {
  fprintf(stderr, "ERROR: %s\n", msg);
  exit(2);
}

void printUsageAndExit() {
  fprintf(stderr, "Usage: count_kmer [options] <FASTA file name> <SAM file name>\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t--kmer,-k <SIZE>\tSpecifies the k-mer size\n");
  fprintf(stderr, "\n");
  exit(2);
}

void countKmerFrequencies (
  const char* FASTAFileName,
  const char* SAMFileName,
  uint KmerSize
)
{
  fprintf(stderr, "===Parameters===\n");
  fprintf(stderr, "Reference FASTA: %s\n", FASTAFileName);
  fprintf(stderr, "Input SAM file : %s\n", SAMFileName);
  fprintf(stderr, "K-mer size     : %d\n", KmerSize);
  fprintf(stderr, "================\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "Loading from the FASTA file ...\r");
  const MultiFASTA multiFASTA = loadFromFASTA(FASTAFileName);
  fprintf(stderr, "Done                           \n");
  FastTSVParse ftp(SAMFileName);
  if(!ftp) {
    fprintf(stderr, "ERROR: Cannot open SAM file '%s'\n", SAMFileName);
    exit(2);
  }


  int tablesize = 1;
  for(int i = 0; i < KmerSize; i++){
    tablesize *= 5;
  }

  int *kmer_kmer_table;
  kmer_kmer_table = new int[tablesize * tablesize];
  for(int i = 0; i < tablesize * tablesize; i++){
    kmer_kmer_table[i] = 0;
  }
  int *kmer_table;
  kmer_table = new int[tablesize];
  for(int i = 0; i < tablesize;  i++){
    kmer_table[i] = 0;
  }

  SAMRecord record;
  size_t recordCount = 0;
  while(ftp.readNextLine()){
    const char* line = ftp.c_str();
    const bool isEmptyLine = line[0] == '\0';
    if(isEmptyLine) continue;
    const bool isCommentLine = line[0] == '@';
    if(isCommentLine) continue;
    if(!record.fill(ftp)) {
      cerr << "SAM Parsing ERROR at line " << ftp.getLineNumber();
      exit(2);
    }
    if(record.rname == "*") continue;
    if(record.seq == "*") continue;
    // use only primary alignment and supplementary alignment and reverse compliment of them.
    // sample has many supplimentary alignment.
    if((record.flag & 2064) != record.flag) continue;


    if(!multiFASTA.count(record.rname.c_str())) {
      cerr << "SAM record says RNAME = '" << record.rname << "', but the reference genome does not have '" << record.rname << "'" << endl;
      exit(2);
    }


    // get aligned sequence at here

    const CIGAROPS cops     = parseCIGARString(record.cigar);

    BString refBS           = multiFASTA.at(record.rname);
    BString queryBS         = String2BString(record.seq);
    const int refStartPos   = record.pos;
    const int queryStartPos = 0; // generateAlignmentSequencesFromCIGARAndSeqs() will manege first Softclip / Hardclip
    BString ras, qas;
    generateAlignmentSequencesFromCIGARAndSeqs(refBS, queryBS, cops, refStartPos, queryStartPos, ras, qas);

    // if record flag has revcomp flag, modify aligned reference sequence and read sequence.
    if((record.flag & 16) != 16){
      revCompBString(ras);
      revCompBString(qas);
    }
    // count up # of kmer-kmer occurence.
    size_t aligned_len = BString2String(ras).size();
    for(size_t idx = 0; idx < aligned_len - KmerSize; idx++){
      //ref_kmer = ras[idx]~ras[idx + KmerSize];
      //KInt<KmerSize> refidx("AAA");

    }
    // divide # of

    cerr << ++recordCount << " processed\r" << flush;
  }
  cerr << endl << "Done." << endl;
}


int main(int argc, char *argv[]){
  GDB_On_SEGV g(argv[0]);

  struct option longopts[] = {
    // { "add",    no_argument,       NULL, 'k' },
    // { "delete", optional_argument, NULL, 'd' },
    { "kmer",  required_argument, NULL,     'k' },
    { 0,        0,                   0,      0  },
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
  const int NUM_REQUIRED_ARGUMENTS = 2;
  if(optind + NUM_REQUIRED_ARGUMENTS != argc) {
    printUsageAndExit();
  }
  const char* fasta_file_name = argv[optind + 0];
  const char* sam_file_name   = argv[optind + 1];

  countKmerFrequencies(fasta_file_name, sam_file_name, kmer_size);
  return 0;
}

