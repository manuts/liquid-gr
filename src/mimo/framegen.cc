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

#define __SAVE__ 1
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
      sps = (std::complex<float> *)malloc(sizeof(std::complex<float>)*k);

      if(__SAVE__){
        std::complex<float> z;
        FILE * f_pn1;
        FILE * f_pn2;
        FILE * f_pn3;
        f_pn1 = fopen("/tmp/tx_pn1", "wb");
        f_pn2 = fopen("/tmp/tx_pn2", "wb");
        f_pn3 = fopen("/tmp/tx_pn3", "wb");
        for(unsigned int i = 0; i < seq_len; i++){
          z = pn1[i] + liquid::math::I*pn2[i];
          fwrite((void *)(pn1 + i), sizeof(std::complex<float>), 1, f_pn1);
          fwrite((void *)(pn2 + i), sizeof(std::complex<float>), 1, f_pn2);
          fwrite((void *)&z, sizeof(std::complex<float>), 1, f_pn3);
        }
        fclose(f_pn1);
        fclose(f_pn2);
        fclose(f_pn3);
      }
      reset();
    }

    unsigned int framegen::get_pn_len() {
      return seq_len;
    }

    framegen::~framegen() {
      free(pn1);
      free(pn2);
      free(sps);
      firinterp_crcf_destroy(interp);
    }

    unsigned int framegen::work(std::complex<float> ** tx_sig,
                                unsigned int num_output)
    {
      unsigned int count = 0;
      // flush out any remaining samples from the prevous symbol
      for(; (sps_count < k)
          && (count < num_output); sps_count++)
      {
        tx_sig[0][count] = gain*(sps[sps_count].real()) + liquid::math::Z;
        tx_sig[1][count] = liquid::math::Z + liquid::math::I*(gain*(sps[sps_count].imag()));
        count++;
      }
      // tx new symbols
      while(count < num_output)
      {
        switch(state) {
          case(STATE_TXPN1):
            firinterp_crcf_execute(interp, pn1[pn_count], sps);
            break;
          case(STATE_PNZ1):
            firinterp_crcf_execute(interp, liquid::math::Z, sps);
            break;
          case(STATE_TXPN2):
            firinterp_crcf_execute(interp, pn2[pn_count], sps);
            break;
          case(STATE_PNZ2):
            firinterp_crcf_execute(interp, liquid::math::Z, sps);
            break;
          case(STATE_TXPN3):
            firinterp_crcf_execute(interp, pn1[pn_count] + 
                liquid::math::I*pn2[pn_count], sps);
            break;
          case(STATE_PNZ3):
            firinterp_crcf_execute(interp, liquid::math::Z, sps);
            break;
        }
        pn_count += 1;
        for(sps_count = 0; (sps_count < k)
            && (count < num_output); sps_count++)
        {
          tx_sig[0][count] = gain*(sps[sps_count].real()) + liquid::math::Z;
          tx_sig[1][count] = liquid::math::Z + liquid::math::I*(gain*(sps[sps_count].imag()));
          count++;
        }
        if(pn_count == seq_len)
        {
          frame_count += 1;
          switch(state) {
            case(STATE_TXPN1):
              state = STATE_PNZ1;
              break;
            case(STATE_PNZ1):
              state = STATE_TXPN2;
              break;
            case(STATE_TXPN2):
              state = STATE_PNZ2;
              break;
            case(STATE_PNZ2):
              state = STATE_TXPN3;
              break;
            case(STATE_TXPN3):
              state = STATE_PNZ3;
              break;
            case(STATE_PNZ3):
              state = STATE_TXPN1;
              break;
          }
          pn_count = 0;
        }
      }
      return count;
    }

    void framegen::reset() {
      frame_count = 0;
      pn_count = 0;
      state = STATE_TXPN1;
      firinterp_crcf_reset(interp);
      sps_count = k;
    }

    unsigned int framegen::get_num_frames()
    {
      return frame_count;
    }
  }
}
