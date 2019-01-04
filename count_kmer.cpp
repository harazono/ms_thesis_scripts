#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <getopt.h>
#include "kmer_library.h"

using namespace std;

void error(const char* msg) {
  fprintf(stderr, "ERROR: %s\n", msg);
  exit(2);
}

void print_usage_and_exit() {
  fprintf(stderr, "Usage: count_kmer [options] <FASTA file name> <SAM file name>\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t--kmer,-k <SIZE>\tSpecifies the k-mer size\n");
  fprintf(stderr, "\n");
  exit(2);
}

void count_kmer_frequencies(
  const char* fasta_file_name,
  const char* sam_file_name,
  uint kmer_size
)
{
  fprintf(stderr, "===Parameters===\n");
  fprintf(stderr, "Reference FASTA: %s\n", fasta_file_name);
  fprintf(stderr, "Input SAM file : %s\n", sam_file_name);
  fprintf(stderr, "K-mer size     : %d\n", kmer_size);
  fprintf(stderr, "================\n");
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
    print_usage_and_exit();
  }
  const char* fasta_file_name = argv[optind + 0];
  const char* sam_file_name   = argv[optind + 1];

  count_kmer_frequencies(fasta_file_name, sam_file_name, kmer_size);
  return 0;
}
