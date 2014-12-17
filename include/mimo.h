/*
 * Copyright (c) 2014 Manu T S
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include <math.h>
#include <complex>
#include <stdlib.h>
#include <liquid/liquid.h>

namespace liquid {
  namespace mimo {
    class framegen
    {
      private:
        msequence ms1;
        msequence ms2;
        unsigned int pn_len_exp;
        unsigned int pn_len;
        std::complex<float> * pn1;
        std::complex<float> * pn2;
      public:
        framegen();
        ~framegen();
        unsigned int get_pn_len();
        void work(std::complex<float> *, std::complex<float> *, unsigned int num_output);
    };
  }
}
