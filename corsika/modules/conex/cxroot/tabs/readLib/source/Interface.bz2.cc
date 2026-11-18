/*
 * (c) Copyright 2020 CORSIKA Project, corsika-project@lists.kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika_data/Interface.h>

#include <bzstream.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

namespace corsika_data {

struct DataFile {

  std::istream *file_in_ = 0;
  // ibzstream file_in_;
  std::string file_name_{};
  DataFile() {}
  ~DataFile() {}
  void Open(const char *name) {
    std::ifstream *if_file = dynamic_cast<std::ifstream *>(file_in_);
    if (if_file && if_file->is_open()) {
      std::cout << "DataFile is still open! Closing it now..." << std::endl;
      if_file->close();
      delete file_in_;
      file_in_ = 0;
    }
    std::string s(name);
    auto i1 = s.find_first_not_of(" \t\r\n"); // trimm right
    auto i2 = s.find_first_of(" \t\r\n", i1); // trimm left
    std::string trimmed(s.substr(i1, i2 - i1));
    file_name_ = trimmed;
    if (trimmed.rfind(".bz2") == trimmed.length() - 4) {
      file_in_ = new ibzstream(trimmed, ios_base::in | ios_base::binary);
    } else {
      file_in_ = new ifstream(trimmed, ios_base::in | ios_base::binary);
    }
    if (!file_in_ || !file_in_->good()) {
      std::cerr << std::string("Cannot open file: \'") + file_name_ << "\'"
                << std::endl;
      throw std::runtime_error(std::string("Cannot open file: ") + file_name_);
    }
  }
  void Close() {
    if (file_in_) {
      delete file_in_;
      file_in_ = 0;
    }
  }
  void FillArray(double *data, const int &length) {
    for (int i = 0; i < length; i++) {
      std::string token;
      (*file_in_) >> token;
      if(token == "NaN") std::cout << "readLib: warning! NaN in table!" << std::endl;
      std:stringstream tokenStream{token};
      tokenStream >> data[i];
      if (!file_in_->good()) {
        throw std::runtime_error(
            std::string("Cannot FillArray i/length=") + std::to_string(i) +
            "/" + std::to_string(length) + " in file " + file_name_);
      }
    }
  }
  double ReadNextNumber() {
    double data;
    (*file_in_) >> data;
    if (!file_in_->good()) {
      throw std::runtime_error(
          std::string("Cannot ReadNextNumber in file " + file_name_));
    }
    return data;
  }
  int ReadNextText(std::string &data) {
    std::getline(*file_in_, data, ' ');
    if (!file_in_->good()) {
      throw std::runtime_error(
          std::string("Cannot ReadNextText in file " + file_name_));
    }
    return 0;
  }
};

static DataFile global_DataFile{};

void CorDataOpenFile(const std::string &name) {
  global_DataFile.Open(name.c_str());
}
void CorDataFillArray(double *data, const int &length) {
  global_DataFile.FillArray(data, length);
}
void CorDataCloseFile() { global_DataFile.Close(); }
double CorDataNextNumber() { return global_DataFile.ReadNextNumber(); }
void CorDataNextText(std::string &data) { global_DataFile.ReadNextText(data); }
void CorDataNextText(char *data, const int maxlength) {
  std::string STR;
  CorDataNextText(STR);
  int i = 0;
  for (; i < maxlength && i < (int)STR.size(); ++i)
    data[i] = STR[i];
  for (; i < maxlength; ++i)
    data[i] = 0;
}
bool CorDataCanDeCompress() { return true; }

extern "C" {
void cordataopenfile_(const char *name, const int stringlength) {
  // make a 'real' string out of fortran one (if \0 is missing)
  char *buffer = new char[stringlength + 1];
  for (int i = 0; i < stringlength; ++i)
    buffer[i] = name[i];
  // if not already terminated by 0, add it
  if (buffer[stringlength - 1] != '\0')
    buffer[stringlength] = '\0';
  CorDataOpenFile(buffer);
  delete[] buffer;
}
void cordatafillarray_(double *data, const int &length) {
  global_DataFile.FillArray(data, length);
}
double cordatanextnumber_() { return global_DataFile.ReadNextNumber(); }
void cordatanexttext_(char *data, const int maxlength) {
  return CorDataNextText(data, maxlength);
}
void cordataclosefile_() { global_DataFile.Close(); }
void cordatacandecompress_(int &v) { v = 1; }
}
} // namespace corsika_data
