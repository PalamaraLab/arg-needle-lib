/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023-2024 ARG-Needle Developers.

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

#include <iostream>

#include "file_utils.hpp"

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

using std::cerr;
using std::endl;

namespace file_utils {

/***** AutoGzOfstream class implementation *****/

void AutoGzOfstream::openOrExit(const std::string& file, std::ios_base::openmode mode) {
  fout.open(file.c_str(), mode);
  if (!fout) {
    cerr << "ERROR: Unable to open file: " << file << endl;
    exit(1);
  }
  if ((int) file.length() > 3 && file.substr(file.length() - 3) == ".gz") {
    boost_out.push(boost::iostreams::gzip_compressor());
  }
  boost_out.push(fout);
}

void AutoGzOfstream::close() {
  boost_out.reset();
}

AutoGzOfstream& AutoGzOfstream::operator<<(std::ostream& (*manip)(std::ostream&) ) {
  manip(boost_out);
  return *this;
}

void AutoGzOfstream::unsetf(std::ios_base::fmtflags mask) {
  boost_out.unsetf(mask);
}

AutoGzOfstream::operator bool() const {
  return !boost_out.fail();
}

} // namespace file_utils
