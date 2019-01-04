#include "gtest/gtest.h"
#include "kmer_library.h"

class MyLibraryTest : public ::testing::Test {
  protected:
        // データメンバーの初期化
  virtual void SetUp() {
  }
};


TEST_F(MyLibraryTest, KInt_Simple1_Test) {
  KInt<1> k1;
  EXPECT_EQ(k1, 0);
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
