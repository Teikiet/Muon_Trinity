/*
 * (c) Copyright 2020 CORSIKA Project, corsika-project@lists.kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <stdexcept>
#include <string>

namespace corsika_data {

  void CorDataOpenFile(const std::string&) {
    throw std::runtime_error(
        "CorDataOpenFile: Cannot read compressed data files with dummy library.");
  }
  void CorDataFillArray(double*, const int&) {
    throw std::runtime_error(
        "CorDataFillArray: Cannot read compressed data files with dummy library.");
  }
  void CorDataCloseFile() {
    throw std::runtime_error(
        "CorDataCloseFile: Cannot read compressed data files with dummy library.");
  }
  double CorDataNextNumber() {
    throw std::runtime_error(
        "CorDataNextNumber: Cannot read compressed data files with dummy library.");
    return 0;
  }
  void CorDataNextText(std::string&) {
    throw std::runtime_error(
        "CorDataNextText(string&): Cannot read compressed data files with dummy "
        "library.");
  }
  void CorDataNextText(char*, int) {
    throw std::runtime_error(
        "CorDataNextText(char*): Cannot read compressed data files with dummy library.");
  }
  bool CorDataCanDeCompress() { return false; }

  // the fortran interface functions
  extern "C" {
  void cordataopenfile_(const char* name, const int) { CorDataOpenFile(name); }
  void cordatafillarray_(double* data, const int& length) {
    CorDataFillArray(data, length);
  }
  void cordataclosefile_() { CorDataCloseFile(); }
  double cordatanextnumber_() { return CorDataNextNumber(); }
  void cordatanexttext_(char*, int) {
    throw std::runtime_error(
        "cordatanexttext_: Cannot read compressed data files with dummy library.");
  }
  void cordatacandecompress_(int& v) { v = 0; }
  }
} // namespace corsika_data
