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

#include "alamouti.h"

namespace liquid {
  namespace alamouti {
    framegen::framegen(unsigned int _k,
                       unsigned int _m,
                       float _beta) {
      training_seq_len = 63;
      payload_len = 128;
      pn1 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len);
      pn2 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len);
      payload = (std::complex<float> *)malloc(sizeof(std::complex<float>)*payload_len);
      msequence ms1 = msequence_create(6, 0x005b, 1);
      msequence ms2 = msequence_create(6, 0x0043, 1);
      msequence ms3 = msequence_create(7, 0x0089, 1);
      for (unsigned int i = 0; i < training_seq_len; i++){
        pn1[i] = (msequence_advance(ms1)) ? 1.0f : -1.0f;
        pn2[i] = (msequence_advance(ms2)) ? 1.0f : -1.0f;
      }
      for (unsigned int i = 0; i < payload_len; i++) {
        payload[i] = (msequence_advance(ms3)) ? 1.0f : -1.0f;
      }
      msequence_destroy(ms1);
      msequence_destroy(ms2);
      msequence_destroy(ms3);

      k = _k;
      m = _m;
      beta = _beta;
      gain1 = 0.25;
      gain2 = 0.25;

      frame_len = k*(2*training_seq_len + payload_len + m);

      interp1 = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);
      interp2 = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);
      sps1 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*k);
      sps2 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*k);

      assert(payload_len%2 == 0);

      reset();
    }

    void framegen::set_gains(float g1, float g2)
    {
      gain1 = g1;
      gain2 = g2;
    }

    unsigned int framegen::get_training_seq_len() {
      return training_seq_len;
    }

    unsigned int framegen::get_frame_len() {
      return frame_len;
    }

    framegen::~framegen() {
      free(pn1);
      free(pn2);
      free(sps1);
      free(sps2);
      firinterp_crcf_destroy(interp1);
      firinterp_crcf_destroy(interp2);
    }

    unsigned int framegen::work(std::complex<float> ** tx_sig)
    {
      unsigned int count = 0;

      // transmit pn1 on channel 1
      for(unsigned int i = 0; i < training_seq_len; i++)
      {
        firinterp_crcf_execute(interp1, pn1[i], sps1);
        firinterp_crcf_execute(interp2, liquid::math::Z, sps2);
        for(unsigned int sps_count = 0; (sps_count < k); sps_count++) {
          tx_sig[0][count] = sps1[sps_count];
          tx_sig[1][count] = sps2[sps_count];
          count++;
        }
      }

      // transmit pn2 on channel 2
      for(unsigned int i = 0; i < training_seq_len; i++)
      {
        firinterp_crcf_execute(interp1, liquid::math::Z, sps1);
        firinterp_crcf_execute(interp2, pn2[i], sps2);
        for(unsigned int sps_count = 0; (sps_count < k); sps_count++) {
          tx_sig[0][count] = sps1[sps_count]*gain1;
          tx_sig[1][count] = sps2[sps_count]*gain2;
          count++;
        }
      }

      // transmit alamouti encoded symbols
      for(unsigned int i = 0; i < payload_len; i += 2)
      {
        // send u1 u2
        firinterp_crcf_execute(interp1, payload[i + 0], sps1);
        firinterp_crcf_execute(interp2, payload[i + 1], sps2);
        for(unsigned int sps_count = 0; (sps_count < k); sps_count++) {
          tx_sig[0][count] = sps1[sps_count]*gain1;
          tx_sig[1][count] = sps2[sps_count]*gain2;
          count++;
        }
        // send -conj(u2) + conj(u1)
        firinterp_crcf_execute(interp1, -std::conj(payload[i + 1]), sps1);
        firinterp_crcf_execute(interp2, std::conj(payload[i + 0]), sps2);
        for(unsigned int sps_count = 0; (sps_count < k); sps_count++) {
          tx_sig[0][count] = sps1[sps_count]*gain1;
          tx_sig[1][count] = sps2[sps_count]*gain2;
          count++;
        }
      }

      // send m zeros to settle
      for(unsigned int i = 0; i < m; i++)
      {
        firinterp_crcf_execute(interp1, liquid::math::Z, sps1);
        firinterp_crcf_execute(interp2, liquid::math::Z, sps2);
        for(unsigned int sps_count = 0; (sps_count < k); sps_count++) {
          tx_sig[0][count] = sps1[sps_count]*gain1;
          tx_sig[1][count] = sps2[sps_count]*gain2;
          count++;
        }
      }
      return count;
    }

    void framegen::reset() {
      firinterp_crcf_reset(interp1);
      firinterp_crcf_reset(interp2);
    }
  }
}
