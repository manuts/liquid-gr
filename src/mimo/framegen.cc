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

#include "mimo.h"

namespace liquid {
  namespace mimo {
    framegen::framegen() {
      pn_len_exp = 7;
      pn_len = (unsigned int)(pow(2, pn_len_exp)) - 1;
      pn1 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*pn_len);
      pn2 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*pn_len);
      ms1 = msequence_create(pn_len_exp, 0x0089, 1);
      ms2 = msequence_create(pn_len_exp, 0x00FD, 1);
    }
    unsigned int framegen::get_pn_len() {
      return pn_len;
    }
    framegen::~framegen() {
      free(pn1);
      free(pn2);
      msequence_destroy(ms1);
      msequence_destroy(ms2);
    }
    void framegen::work(std::complex<float> * chan1_symb,
                        std::complex<float> * chan2_symb,
                        unsigned int num_output)
    {
      for(unsigned int i = 0; i < num_output; i++)
      {
        chan1_symb[i] = (msequence_advance(ms1) ? 1.0f : -1.0f);
        chan2_symb[i] = (msequence_advance(ms2) ? 1.0f : -1.0f);
      }
    }
  }
}
