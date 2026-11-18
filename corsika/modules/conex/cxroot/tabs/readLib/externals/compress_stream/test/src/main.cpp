// ============================================================================
// compress_stream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2018 Pierre Marijon <pierre@marijon.fr>
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================

/* standard include */

/* project include */
#include "gzstream.hpp"
#include "bzstream.hpp"

#include <algorithm>
#include <iterator>
#include <string>
#include <fstream>

template <typename InputIterator1, typename InputIterator2>
bool range_equal(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2,
                 InputIterator2 last2) {
  while (first1 != last1 && first2 != last2) {
    if (*first1 != *first2) return false;
    ++first1;
    ++first2;
  }
  return (first1 == last1) && (first2 == last2);
}

bool compare_files(const std::string& filename1, const std::string& filename2) {
  std::ifstream file1(filename1);
  std::ifstream file2(filename2);

  std::istreambuf_iterator<char> begin1(file1);
  std::istreambuf_iterator<char> begin2(file2);

  std::istreambuf_iterator<char> end;

  return range_equal(begin1, end, begin2, end);
}

void test_out(bool result, std::string message) {
  if (result) {
    std::cout << "\033[32m" << message << ": Passed"
              << "\033[0m" << std::endl;
  } else {
    std::cout << "\033[31m" << message << ": Failled"
              << "\033[0m" << std::endl;
  }
}

int main(int argc, char** argv) {
  std::string message =
      "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Etiam vitae erat sit "
      "amet lorem vehicula malesuada. Praesent tempus, tortor sed dapibus molestie, "
      "tortor sem sagittis velit, eget dapibus est odio aliquet magna. Morbi sed "
      "malesuada quam, nec rhoncus ligula. Ut non fermentum metus. Suspendisse potenti. "
      "Nam auctor facilisis lorem fermentum.";

  /* write a message */
  ogzstream out_gz("out.gz");
  obzstream out_bz("out.bz");

  out_gz << message << std::endl;
  out_bz << message << std::endl;

  out_gz.close();
  out_bz.close();

  /* read message */
  igzstream in_gz("out.gz");
  ibzstream in_bz("out.bz");

  std::string gzread_message;
  std::getline(in_gz, gzread_message);

  std::string bzread_message;
  std::getline(in_bz, bzread_message);

  /* check */
  /* created file to reference file */
  test_out(compare_files("out.gz", "reference.gz"), "gzip file same as reference");
  test_out(compare_files("out.bz", "reference.bz"), "bzip file same as reference");

  /* read message to write message */
  test_out(gzread_message == message, "gz message read same as original");
  test_out(bzread_message == message, "bz message read same as original");
}
