/*
 * (c) Copyright 2020 CORSIKA Project, corsika-project@lists.kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#ifdef STANDALONE
#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this
                          // in one cpp file
#endif
#include <catch2/catch_all.hpp>

#include <corsika_data/Interface.h>

#include <array>
#include <iostream>
#include <string>
using namespace std;
using namespace corsika_data;

const array<double, 3 * 9> testData = {
    0.00000000,     0.00000000,     0.00000000,     0.00000000,
    0.00000000,     0.00000000,     0.00000000,     0.00000000,
    0.00000000,     7.74900036E-06, 1.00000000,     1.47967212E-05,
    3.76063472E-05, 3.91302201E-05, 5.68050455E-05, 5.91166172E-05,
    6.64338659E-05, 8.79753206E-05, 9.61698461E-05, 1.02615748E-04,
    1.12829475E-04, 1.28208514E-04, 1.50734515E-04, 1.79381183E-04,
    2.09064019E-04, 2.32736333E-04, 2.50067125E-04};

const std::string fileName = std::string(TESTDATA);
const std::string fileNameBZ2 = std::string(TESTDATA) + std::string(".bz2");

TEST_CASE("Data", "[data]") {

  auto file = GENERATE(as<std::string>{}, fileName, fileNameBZ2);

  SECTION(std::string("c++,") + file) {
    bool b = CorDataCanDeCompress();
    CHECK(b == true);
    cout << "Open/Reading: " << file << endl;
    CorDataOpenFile(file);
    std::string str;
    CorDataNextText(str);
    CHECK(str == "sibyll20");
    CHECK(CorDataNextNumber() == 0.10000E-02);
    CHECK(CorDataNextNumber() == 91);
    CHECK(CorDataNextNumber() == 261);
    CHECK(CorDataNextNumber() == 5);
    CHECK(CorDataNextNumber() == 11);
    CHECK(CorDataNextNumber() == 20);

    double aData[3 * 9];
    const int length = 3 * 9;
    CorDataFillArray(aData, length);

    for (int i = 0; i < length; ++i) {
      CHECK(aData[i] == testData[i]);
    }
    CorDataCloseFile();
  }

  SECTION(std::string("fortran, ") + file) {
    int i = 0;
    cordatacandecompress_(i);
    CHECK(i == 1);
    cout << "Open/Reading: " << file << endl;
    cordataopenfile_(file.c_str(), file.length());
    char str[10];
    cordatanexttext_(str, 10);
    CHECK(std::string(str) == "sibyll20");
    CHECK(cordatanextnumber_() == 0.10000E-02);
    CHECK(cordatanextnumber_() == 91);
    CHECK(cordatanextnumber_() == 261);
    CHECK(cordatanextnumber_() == 5);
    CHECK(cordatanextnumber_() == 11);
    CHECK(cordatanextnumber_() == 20);

    double aData[3 * 9];
    const int length = 3 * 9;
    cordatafillarray_(aData, length);

    for (int i = 0; i < length; ++i) {
      CHECK(aData[i] == testData[i]);
    }
    cordataclosefile_();
  }

  // now read file two times in a row
  SECTION(std::string("twice, c++,") + file) {
    bool b = CorDataCanDeCompress();
    CHECK(b == true);
    cout << "Open/Close/Open/Reading: " << file << endl;
    CorDataOpenFile(file);
    CorDataCloseFile();
    CorDataOpenFile(file);
    std::string str;
    CorDataNextText(str);
    CHECK(str == "sibyll20");
    CHECK(CorDataNextNumber() == 0.10000E-02);
    CHECK(CorDataNextNumber() == 91);
    CHECK(CorDataNextNumber() == 261);
    CHECK(CorDataNextNumber() == 5);
    CHECK(CorDataNextNumber() == 11);
    CHECK(CorDataNextNumber() == 20);

    double aData[3 * 9];
    const int length = 3 * 9;
    CorDataFillArray(aData, length);

    for (int i = 0; i < length; ++i) {
      CHECK(aData[i] == testData[i]);
    }
    CorDataCloseFile();
  }

  SECTION(std::string("twice, fortran, ") + file) {
    int i = 0;
    cordatacandecompress_(i);
    CHECK(i == 1);
    cout << "Open/Close/Open/Reading: " << file << endl;
    cordataopenfile_(file.c_str(), file.length());
    cordataclosefile_();
    cordataopenfile_(file.c_str(), file.length());
    char str[10];
    cordatanexttext_(str, 10);
    CHECK(std::string(str) == "sibyll20");
    CHECK(cordatanextnumber_() == 0.10000E-02);
    CHECK(cordatanextnumber_() == 91);
    CHECK(cordatanextnumber_() == 261);
    CHECK(cordatanextnumber_() == 5);
    CHECK(cordatanextnumber_() == 11);
    CHECK(cordatanextnumber_() == 20);

    double aData[3 * 9];
    const int length = 3 * 9;
    cordatafillarray_(aData, length);

    for (int i = 0; i < length; ++i) {
      CHECK(aData[i] == testData[i]);
    }
    cordataclosefile_();
  }
}
