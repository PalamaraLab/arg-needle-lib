/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023-2025 ARG-Needle Developers.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// The below code was copied and modified from the Eagle software file
// https://github.com/poruloh/Eagle/blob/master/src/FileUtils.hpp
// developed by Po-Ru Loh and released under the GNU General Public
// License v3.0 (GPLv3).
//
// The license file can be found at 3rd_party/Eagle/COPYING from the
// root of this repository.

#ifndef __FILE_UTILS_HPP_
#define __FILE_UTILS_HPP_

#include <fstream>
#include <string>

#include <boost/iostreams/filtering_stream.hpp>

namespace file_utils {

class AutoGzOfstream {

  boost::iostreams::filtering_ostream boost_out;
  std::ofstream fout;

public:
  void openOrExit(const std::string& file, std::ios_base::openmode mode = std::ios::out);
  void close();
  template <class T> AutoGzOfstream& operator<<(const T& x) {
    boost_out << x;
    return *this;
  }

  AutoGzOfstream& operator<<(std::ostream& (*manip)(std::ostream&) );
  void unsetf(std::ios_base::fmtflags);
  operator bool() const;
};

} // namespace file_utils

#endif // __FILE_UTILS_HPP_
