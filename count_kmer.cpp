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

void printUsageAndExit() {
  fprintf(stderr, "Usage: count_kmer [options] <FASTA file name> <SAM file name>\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t--kmer,-k <SIZE>\tSpecifies the k-mer size\n");
  fprintf(stderr, "\n");
  exit(2);
}

template<uint KMERSIZE>
struct FrequencyTable {
  static const size_t tablesize = ipow(5, KMERSIZE);
  typedef int Frequency;
  typedef int Score;
  vector<Frequency> kmer_table; ///< Reference side
  vector<Frequency> kmer_kmer_table; ///< The first dimension is reference, the second is query.
  vector<double>    kmer_kmer_prob_table;
  vector<Score>     score_table;
  //vector<Score>     score_table(tablesize * tablesize, 0); ///< divided by kmer_table
  inline Frequency& kk(size_t r, size_t q) {
    MYASSERT_WMD("Out of range (r)", r < tablesize, DUMP(r));
    MYASSERT_WMD("Out of range (q)", q < tablesize, DUMP(q));
    return kmer_kmer_table[r * tablesize + q];
  }
  inline double&   kkp(size_t r, size_t q) {
    MYASSERT_WMD("Out of range (r)", r < tablesize, DUMP(r));
    MYASSERT_WMD("Out of range (q)", q < tablesize, DUMP(q));
    return kmer_kmer_prob_table[r * tablesize + q];
  }

  inline Score&    scr(size_t r, size_t q) {
    MYASSERT_WMD("Out of range (r)", r < tablesize, DUMP(r));
    MYASSERT_WMD("Out of range (q)", q < tablesize, DUMP(q));
    return score_table[r * tablesize + q];
  }

  void outputAsBinaryTable(const string& outputFileName)
  {
    const char *fn = outputFileName.c_str();
    FILE *fle = fopen(fn, "wb");
    if(fle == NULL) {
      fprintf(stderr, "ERROR: Cannot open binary table file '%s'\n", fn);
      exit(2);
    }
    //// FILE FORMAT:
    ////  offset  0: 'BINFREQT' (8 bytes)
    ////  offset  8: k-mer size (8 bytes)
    ////  offset 16: sizeof(Frequency) (8 bytes)
    ////  offset 24: table_size (8 bytes)
    ////  offset 32: table (8 x table_size bytes)
    static const char MAGIC_STRING[] = "BINFREQT";
    MYASSERT_WMD("MAGIC_STRING must be 8 bytes", sizeof(MAGIC_STRING) - 1 == 8, DUMP(MAGIC_STRING));
    fwrite(MAGIC_STRING, sizeof(MAGIC_STRING) - 1, 1, fle);
    size_t kmer_size = KMERSIZE;
    MYASSERT_WMD("sizeof(size_t) must be 8", sizeof(size_t) == 8, DUMP(sizeof(size_t)));
    fwrite(&kmer_size, sizeof(kmer_size), 1, fle);
    size_t size_of_frequency = sizeof(Frequency);
    fwrite(&size_of_frequency, sizeof(size_of_frequency), 1, fle);
    size_t var_table_size = tablesize * tablesize;
    fwrite(&var_table_size, sizeof(var_table_size), 1, fle);
    //MYASSERT_WMD("var_table_size == kmer_kmer_table.size()", var_table_size == kmer_kmer_table.size(), DUMP(var_table_size, tablesize, kmer_kmer_table.size()));
    fwrite(&*kmer_kmer_table.begin(), sizeof(Frequency) * kmer_kmer_table.size(), 1, fle);
    fclose(fle);
  }
  void countAlignment(
    const BString& ras,
    const BString& qas
  )
  {
    MYASSERT_WMD("ras.size() must be qas.size()", ras.size() == qas.size(), DUMP(ras.size(), qas.size()));
    const size_t aligned_len = ras.size();
    if(aligned_len < KMERSIZE) return;
    KInt<KMERSIZE> rki, qki;
    for(size_t i = 0; i < KMERSIZE - 1; ++i) {
      rki.ShiftIn(ras[i]);
      qki.ShiftIn(qas[i]);
    }
    for(size_t i = KMERSIZE - 1; i <= aligned_len - KMERSIZE; i++){
      rki.ShiftIn(ras[i]);
      qki.ShiftIn(qas[i]);
      kmer_table[rki]++;
      kk(rki, qki)++;
    }
  }
  void printTable(){
    vector<double>    kmer_kmer_prob_table(tablesize * tablesize, 0); ///< divided by kmer_table
    auto kkp = [&kmer_kmer_prob_table](size_t r, size_t q) -> double& {
      MYASSERT_WMD("Out of range (r)", r < tablesize, DUMP(r));
      MYASSERT_WMD("Out of range (q)", q < tablesize, DUMP(q));
      return kmer_kmer_prob_table[r * tablesize + q];
    };
    for(int i = 0; i < tablesize; i++){
      for(int j = 0; j < tablesize; j++){
        kkp(i, j) = static_cast<double>(kk(i, j)) / static_cast<double>(kmer_table[i]);
      }
    }
    for(int i = 0; i < tablesize; i++){
      for(int j = 0; j < tablesize; j++){
        fprintf(stdout, "%lf", log10(kkp(j, i)));
        if(j != tablesize) fprintf(stdout, ", ");
      }
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
  }


  void printKKTable(){
    for(int i = 0; i < tablesize; i++){
      for(int j = 0; j < tablesize; j++){
        fprintf(stdout, "%8d", kk(j, i));
        if(j != tablesize - 1) fprintf(stdout, ", ");
      }
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
  }

  
  void printKTable(){
    for(int i = 0; i < tablesize; i++){
      fprintf(stdout, "%8d, ", kmer_table[i]);
    }
    fprintf(stdout, "\n\n");
  }



  void scorerize(){
    vector<Frequency> kmer_ins_table(tablesize, 0); //<used for normalize gap containing reference line.
    vector<double>    kmer_kmer_prob_table(tablesize * tablesize, 0); ///< divided by kmer_table
    auto kkp = [&kmer_kmer_prob_table](size_t r, size_t q) -> double& {
      MYASSERT_WMD("Out of range (r)", r < tablesize, DUMP(r));
      MYASSERT_WMD("Out of range (q)", q < tablesize, DUMP(q));
      return kmer_kmer_prob_table[r * tablesize + q];
    };
    for(int qi = 0; qi < tablesize; qi++){
      for(int ri = 0; ri < tablesize; ri++){
        kmer_ins_table[qi] += kk(ri, qi);
      }
    }
    // almost all cells are normalized by sum of reference k-mer frequency.
    for(int i = 0; i < tablesize; i++){
      for(int j = 0; j < tablesize; j++){
        kkp(i, j) = static_cast<double>(kk(i, j)) / static_cast<double>(kmer_table[i]);//divide by reference k-mer frequency.
      }
    }
    // cells they contains '-' in reference k-mer should be normalized sum of frequency of query k-mer.
    for(int qi = 0; qi < tablesize; qi++){
      for(int ri = 4; ri < tablesize; ri+=5){
        kkp(ri, qi) = static_cast<double>(kk(ri, qi)) / static_cast<double>(kmer_ins_table[qi]); // divede by query kmer frequency.
      }
    }

    for(int i = 0; i < tablesize; i++){
      for(int j = 0; j < tablesize; j++){
        if(kkp(j, i) != 0){
          scr(j, i) = static_cast<int>(round(100 * log10(kkp(j, i))));
          fprintf(stdout, "%d",static_cast<int>(round(100 * log10(kkp(j, i)))));
        }else{
          scr(j. i) = -1 * ipow(2, 10);
          fprintf(stdout, "%d", -1 * ipow(2, 10));
        }
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "%d", scr(1,1));
  }

public:
  FrequencyTable() : kmer_table(tablesize, 0), kmer_kmer_table(tablesize * tablesize, 0), score_table(tablesize * tablesize, 0) {}
  void countKmerFrequencies (
    const char* FASTAFileName,
    const char* SAMFileName,
    uint KmerSize,
    const bool outputInCSV,
    const string binaryOutputFileName ///< empty() if --binary is not given
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
      countAlignment(ras, qas);
      // divide # of

      ++recordCount;
      if(recordCount % 100 == 0) {
        cerr << recordCount << " processed\r" << flush;
      }
    }
    cerr << recordCount << " processed\n";
    cerr << "Done." << endl;
    scorerize();
    if(outputInCSV) {
      printTable();
    }
    if(!binaryOutputFileName.empty()) {
      outputAsBinaryTable(binaryOutputFileName);
    }
  }


};


int main(int argc, char *argv[]){
  GDB_On_SEGV g(argv[0]);

  struct option longopts[] = {
    { "csv"       , no_argument       , NULL , 'c' } ,
    // { "delete" , optional_argument , NULL , 'd' } ,
    { "kmer"      , required_argument , NULL , 'k' } ,
    { "binary"    , required_argument , NULL , 'b' } ,
    { 0           , 0                 , 0    , 0  }  ,
  };

  /// PARAMETERS ///
  int kmer_size = 1;
  bool output_in_csv = false;
  string binary_output_file_name;
  //////////////////
  int opt;
  int longindex;
  while ((opt = getopt_long(argc, argv, "k:", longopts, &longindex)) != -1) {
    switch (opt) {
    case 'b':
      binary_output_file_name = optarg;
      break;
    case 'c':
      output_in_csv = true;
      break;
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

  #define FT() ft.countKmerFrequencies(fasta_file_name, sam_file_name, kmer_size, output_in_csv, binary_output_file_name)
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

