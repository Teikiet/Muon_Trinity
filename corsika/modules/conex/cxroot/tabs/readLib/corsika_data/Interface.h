/*
 * (c) Copyright 2020 CORSIKA Project, corsika-project@lists.kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <string>

namespace corsika_data {

  // the c++ interface functions
  // void CorDataOpenFile(const char* name);
  void CorDataOpenFile(const std::string& name);
  void CorDataFillArray(double* data, const int& length);
  void CorDataCloseFile();
  double CorDataNextNumber();
  void CorDataNextText(std::string& data);
  void CorDataNextText(char* data, int maxlength);
  bool CorDataCanDeCompress();

  // the fortran interface functions
  extern "C" {
  void cordataopenfile_(const char*, const int stringlength);
  void cordatafillarray_(double*, const int&);
  void cordataclosefile_();
  double cordatanextnumber_();
  void cordatanexttext_(char*, int);
  void cordatacandecompress_(int&);
  }
} // namespace corsika_data
