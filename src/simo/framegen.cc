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

#include "simo.h"

namespace liquid {
  namespace simo {
    framegen::framegen(unsigned int _k,
                       unsigned int _m,
                       float _beta) {
      training_seq_len = 63;
      payload_len = 128;
      pn = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len);
      payload = (std::complex<float> *)malloc(sizeof(std::complex<float>)*payload_len);
      msequence ms1 = msequence_create(6, 0x005b, 1);
      msequence ms3 = msequence_create(7, 0x0089, 1);
      for (unsigned int i = 0; i < training_seq_len; i++){
        pn[i] = (msequence_advance(ms1)) ? 1.0f : -1.0f;
      }
      for (unsigned int i = 0; i < payload_len; i++) {
        payload[i] = (msequence_advance(ms3)) ? 1.0f : -1.0f;
      }
      msequence_destroy(ms1);
      msequence_destroy(ms3);

      k = _k;
      m = _m;
      beta = _beta;
      gain = 0.25;

      frame_len = k*(training_seq_len + payload_len + m);

      interp = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);
      sps = (std::complex<float> *)malloc(sizeof(std::complex<float>)*k);

      assert(payload_len%2 == 0);

      reset();
    }

    void framegen::set_gains(float g)
    {
      gain = g;
    }

    unsigned int framegen::get_training_seq_len() {
      return training_seq_len;
    }

    unsigned int framegen::get_frame_len() {
      return frame_len;
    }

    framegen::~framegen() {
      free(pn);
      free(sps);
      firinterp_crcf_destroy(interp);
    }

    unsigned int framegen::work(std::complex<float> * tx_sig)
    {
      unsigned int count = 0;

      // transmit pn
      for(unsigned int i = 0; i < training_seq_len; i++)
      {
        firinterp_crcf_execute(interp, pn[i], sps);
        for(unsigned int sps_count = 0; (sps_count < k); sps_count++) {
          tx_sig[count] = sps[sps_count]*gain;
          count++;
        }
      }

      // transmit symbols
      for(unsigned int i = 0; i < payload_len; i++)
      {
        // send u1 u2
        firinterp_crcf_execute(interp, payload[i], sps);
        for(unsigned int sps_count = 0; (sps_count < k); sps_count++) {
          tx_sig[count] = sps[sps_count]*gain;
          count++;
        }
      }

      // send m zeros to settle
      for(unsigned int i = 0; i < m; i++)
      {
        firinterp_crcf_execute(interp, liquid::math::Z, sps);
        for(unsigned int sps_count = 0; (sps_count < k); sps_count++) {
          tx_sig[count] = sps[sps_count]*gain;
          count++;
        }
      }
      return count;
    }

    void framegen::reset() {
      firinterp_crcf_reset(interp);
    }
  }
}
