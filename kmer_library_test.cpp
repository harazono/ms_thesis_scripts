#include "gtest/gtest.h"
#include <string>
#include "kmer_library.h"

using namespace std;

class MyLibraryTest : public ::testing::Test {
  protected:
        // データメンバーの初期化
  virtual void SetUp() {
  }
};


TEST_F(MyLibraryTest, KInt_Simple1_Test) {
  KInt<1> k1;
  EXPECT_EQ(k1, 0);
  EXPECT_STREQ(k1.str().c_str(), "A");
  KInt<1> k2("A");
  EXPECT_EQ(k2, 0);
  KInt<1> k3("C");
  EXPECT_EQ(k3, 1);
  KInt<1> k4("G");
  EXPECT_EQ(k4, 2);
  KInt<1> k5("T");
  EXPECT_EQ(k5, 3);
  KInt<1> k6("-");
  EXPECT_EQ(k6, 4);
}

TEST_F(MyLibraryTest, KInt_Simple2_Test) {
  KInt<2> k1;
  EXPECT_EQ(k1, 0);
  KInt<2> k2("AA");
  EXPECT_EQ(k2, 0);
  KInt<2> k3("AC");
  EXPECT_EQ(k3, 1);
  KInt<2> k4("CA");
  EXPECT_EQ(k4, 5);
  KInt<2> k5("CC");
  EXPECT_EQ(k5, 6);
  KInt<2> k6("--");
  EXPECT_EQ(k6, 24);
}

TEST_F(MyLibraryTest, KInt_Simple3_Test) {
  KInt<3> k1;
  EXPECT_EQ(k1, 0);
  KInt<3> k4("GGC");
  EXPECT_EQ(k4, 61);
  KInt<3> k5("G-G");
  EXPECT_EQ(k5, 72);
  KInt<3> k6("-AG");
  EXPECT_EQ(k6, 102);
}

TEST_F(MyLibraryTest, LoadMultiFASTA_Test) {
  MultiFASTA mf = loadFromFASTA("10.ref");
  EXPECT_EQ(mf.size(), 1);
  EXPECT_STREQ(mf.begin()->first.c_str(), "testref");
  EXPECT_STREQ(BString2String(mf.begin()->second).c_str(), "CGACTATTCC");
}


TEST_F(MyLibraryTest, SplitBy1stSpace_Test) {
  std::string head = splitBy1stSpace("AAA WWW");
  EXPECT_STREQ(head.c_str(), "AAA");
  EXPECT_EQ(head.size(), 3);
  head = splitBy1stSpace("i;oarehwv wir efw;ioh");
  EXPECT_STREQ(head.c_str(), "i;oarehwv");
  EXPECT_EQ(head.size(), 9);
}

TEST_F(MyLibraryTest, ParseSAM_Test) {
  FastTSVParse ftp("10.sam");
  ASSERT_FALSE(!ftp);
  while(true) {
    ASSERT_TRUE(ftp.readNextLine());
    if(ftp.c_str()[0] != '@') break;
  }
  SAMRecord r;
  //cout << ftp.c_str() << endl;
  //cout << ftp.size() << endl;
  ASSERT_TRUE(r.fill(ftp));
  EXPECT_STREQ(r.qname.c_str(), "testread");
  EXPECT_STREQ(r.rname.c_str(), "testref");
  EXPECT_EQ(r.flag, 0);
  EXPECT_STREQ(r.cigar.c_str(), "10M");
  EXPECT_STREQ(r.seq.c_str(), "CGTCTATTCC");
  EXPECT_STREQ(r.qual.c_str(), "*");
  EXPECT_STREQ(r.rnext.c_str(), "*");
  EXPECT_EQ(r.tlen, 0);
  EXPECT_EQ(r.pnext, 0);
  EXPECT_EQ(r.mapQ, 60);
  EXPECT_EQ(r.pos, 0);
}

TEST_F(MyLibraryTest, BString_Test) {
  const string s = "ACCGGGTGGTC--G-";
  const BString bs = String2BString(s);
  const string rev_s = BString2String(bs);
  EXPECT_STREQ(s.c_str(), rev_s.c_str());
  EXPECT_EQ(bs[0], 0);
  EXPECT_EQ(bs[1], 1);
  EXPECT_EQ(bs.back(), 4);
}

TEST_F(MyLibraryTest, ParseCIGAR_Simple1_Test) {
  CIGAROPS ops = parseCIGARString("1M");
  ASSERT_EQ(ops.size(), 1);
  EXPECT_EQ(ops[0].op, 'M');
  EXPECT_EQ(ops[0].len, 1);
}

TEST_F(MyLibraryTest, ParseCIGAR_Simple2_Test) {
  CIGAROPS ops = parseCIGARString("4S21M412D");
  ASSERT_EQ(ops.size(), 3);
  EXPECT_EQ(ops[0].op, 'S');
  EXPECT_EQ(ops[0].len, 4);
  EXPECT_EQ(ops[1].op, 'M');
  EXPECT_EQ(ops[1].len, 21);
  EXPECT_EQ(ops[2].op, 'D');
  EXPECT_EQ(ops[2].len, 412);
}

TEST_F(MyLibraryTest, GenerateAlignmentFromCIGAR_Test) {
  const string refSeq     = "ACGTGCGT";
  const string querySeq   = "CGGCAG";
  const BString refBSeq   = String2BString(refSeq);
  const BString queryBSeq = String2BString(querySeq);
  // ACGTGC-GT
  //  || || |
  //  CG-GCAG
  CIGAROPS ops = parseCIGARString("1H2M1D2M1I1M1H");
  BString ras, qas;
  generateAlignmentSequencesFromCIGARAndSeqs(refBSeq, queryBSeq, ops, 1, 0, ras, qas);
  const string rass = BString2String(ras);
  const string qass = BString2String(qas);
  EXPECT_STREQ(rass.c_str(), "CGTGC-G");
  EXPECT_STREQ(qass.c_str(), "CG-GCAG");
}

TEST_F(MyLibraryTest,revComp_Test){
  string preSeq   = "AAAC";
  BString preBSeq = String2BString(preSeq);
  BString revBSeq = revCompBString(preBSeq);
  string ansSeq   = BString2String(revBSeq);
  EXPECT_STREQ(ansSeq.c_str(), "GTTT");

  preSeq   = "A";
  preBSeq  = String2BString(preSeq);
  revBSeq  = revCompBString(preBSeq);
  ansSeq   = BString2String(revBSeq);
  EXPECT_STREQ(ansSeq.c_str(), "T");

  preSeq   = "CGTGAAA";
  preBSeq  = String2BString(preSeq);
  revBSeq  = revCompBString(preBSeq);
  ansSeq   = BString2String(revBSeq);
  EXPECT_STREQ(ansSeq.c_str(), "TTTCACG");
}

int main(int argc, char **argv) {
  GDB_On_SEGV g(argv[0]);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
