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
    framegen::framegen(unsigned int _k,
                       unsigned int _m,
                       float _beta,
                       float _gain) {
      seq_len_exp = 6;
      seq_len = (unsigned int)(pow(2, seq_len_exp)) - 1;
      pn1 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*seq_len);
      pn2 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*seq_len);
      msequence ms1 = msequence_create(seq_len_exp, 0x005b, 1);
      msequence ms2 = msequence_create(seq_len_exp, 0x0043, 1);
      for (unsigned int i = 0; i < seq_len; i++){
        pn1[i] = (msequence_advance(ms1)) ? 1.0f : -1.0f;
        pn2[i] = (msequence_advance(ms2)) ? 1.0f : -1.0f;
      }
      msequence_destroy(ms1);
      msequence_destroy(ms2);
      k = _k;
      m = _m;
      beta = _beta;
      gain = _gain;
      interp = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);

      reset();
    }

    unsigned int framegen::get_pn_len() {
      return seq_len;
    }

    framegen::~framegen() {
      free(pn1);
      free(pn2);
      firinterp_crcf_destroy(interp);
    }

    unsigned int framegen::work(std::vector<std::complex<float>*> tx_sig,
                                unsigned int num_output)
    {
      unsigned int count = 0;
      std::complex<float> * samps = 
        (std::complex<float> *)malloc(sizeof(std::complex<float>)*k);
      while(count + k < num_output)
      {
        switch(state) {
          case(STATE_TXPN1):
            firinterp_crcf_execute(interp, pn1[pn_count], samps);
            break;
          case(STATE_TXPN2):
            firinterp_crcf_execute(interp, pn2[pn_count], samps);
            break;
        }
        pn_count += 1;
        for(unsigned int sps_count = 0; sps_count < k; sps_count++)
        {
          tx_sig[0][count] = gain*(samps[sps_count].real()) + liquid::math::Z;
          tx_sig[1][count] = liquid::math::Z + liquid::math::I*(gain*(samps[sps_count].imag()));
          count++;
        }
        if(pn_count == seq_len)
        {
          switch(state) {
            case(STATE_TXPN1):
              state = STATE_TXPN2;
              break;
            case(STATE_TXPN2):
              frame_count += 1;
              state = STATE_TXPN1;
              break;
          }
          pn_count = 0;
        }
      }
      free(samps);
      return count;
    }

    void framegen::reset() {
      frame_count = 0;
      pn_count = 0;
      state = STATE_TXPN1;
      firinterp_crcf_reset(interp);
    }

    unsigned int framegen::get_num_frames()
    {
      return frame_count;
    }
  }
}
