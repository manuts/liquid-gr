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

#define PRINT_STATS_TX1RX1 true
#define PRINT_STATS_TX1RX2 true
#define PRINT_STATS_TX2RX1 true
#define PRINT_STATS_TX2RX2 true
#define PRINT_STAT_DIFFS false

namespace liquid {
  namespace simo {
    framesync::framesync(unsigned int _k,
                         unsigned int _m,
                         float _beta)
    {
      training_seq_len = 63;
      payload_len = 128;

      msequence ms1 = msequence_create(6, 0x005b, 1);
      preamble_pn = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len); 
      preamble_rx[0] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len); 
      preamble_rx[1] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len); 
      for (unsigned int i = 0; i < training_seq_len; i++){
        preamble_pn[i] = (msequence_advance(ms1)) ? 1.0f : -1.0f;
      }
      msequence_destroy(ms1);

      // interpolate p/n sequence with matched filter
      k     = _k;        // samples/symbol
      m     = _m;        // filter delay (symbols)
      beta  = _beta;    // excess bandwidth factor
      frame_len = k*(training_seq_len + payload_len + m);

      std::complex<float> seq[k*training_seq_len];
      firinterp_crcf interp = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);

      for (unsigned int i = 0; i < training_seq_len; i++) {
        firinterp_crcf_execute(interp, preamble_pn[i%training_seq_len], &seq[k*i]);
      }
      firinterp_crcf_destroy(interp);

      // create frame detector
      float threshold = 0.4f;     // detection threshold
      float dphi_max  = 0.05f;    // maximum carrier offset allowable
      frame_detector[0] = detector_cccf_create(seq, k*training_seq_len, threshold, dphi_max);
      frame_detector[1] = detector_cccf_create(seq, k*training_seq_len, threshold, dphi_max);
      buffer[0] = windowcf_create(frame_len);
      buffer[1] = windowcf_create(frame_len);

      npfb = 32;
      mf[0]   = firpfb_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER, npfb,k,m,beta);
      dmf[0]  = firpfb_crcf_create_drnyquist(LIQUID_FIRFILT_ARKAISER,npfb,k,m,beta);
      mf[1]   = firpfb_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER, npfb,k,m,beta);
      dmf[1]  = firpfb_crcf_create_drnyquist(LIQUID_FIRFILT_ARKAISER,npfb,k,m,beta);

      nco_coarse_freq = 0.0f;
      nco_coarse[0] = nco_crcf_create(LIQUID_NCO);
      nco_coarse[1] = nco_crcf_create(LIQUID_NCO);
      nco_fine[0]   = nco_crcf_create(LIQUID_VCO);
      nco_fine[1]   = nco_crcf_create(LIQUID_VCO);
      nco_crcf_pll_set_bandwidth(nco_fine[0], 0.05f);
      nco_crcf_pll_set_bandwidth(nco_fine[1], 0.05f);

      dphi_hat[0] = 0;
      dphi_hat[1] = 0;

      nco_crcf_reset(nco_coarse[0]);
      nco_crcf_reset(nco_coarse[1]);
      nco_crcf_reset(nco_fine[0]);
      nco_crcf_reset(nco_fine[1]);

      training_seq_n[0] = 0;
      training_seq_n[1] = 0;

      reset();
    }

    void framesync::reset()
    {
      frame_detect_flag[0] = false;
      frame_detect_flag[1] = false;
      state = STATE_DETECTFRAME;

      detector_cccf_reset(frame_detector[0]);
      detector_cccf_reset(frame_detector[1]);

      windowcf_clear(buffer[0]);
      windowcf_clear(buffer[1]);

      // reset symbol timing recovery state
      firpfb_crcf_reset(mf[0]);
      firpfb_crcf_reset(mf[1]);
      firpfb_crcf_reset(dmf[0]);
      firpfb_crcf_reset(dmf[1]);
      update_training_sequence_n();
    }

    framesync::~framesync()
    {
      free(preamble_pn);
      free(preamble_rx[0]);
      free(preamble_rx[1]);

      detector_cccf_destroy(frame_detector[0]);
      detector_cccf_destroy(frame_detector[1]);

      windowcf_destroy(buffer[0]);
      windowcf_destroy(buffer[1]);

      firpfb_crcf_destroy(mf[0]);
      firpfb_crcf_destroy(mf[1]);
      firpfb_crcf_destroy(dmf[0]);
      firpfb_crcf_destroy(dmf[1]);

      nco_crcf_destroy(nco_coarse[0]);
      nco_crcf_destroy(nco_coarse[1]);
      nco_crcf_destroy(nco_fine[0]);
      nco_crcf_destroy(nco_fine[1]);
    }

    void framesync::work(std::complex<float> ** _x,
                         unsigned int    _n)
    {
      for (unsigned int i = 0; i < _n; i++) {
        switch (state) {
          case STATE_DETECTFRAME:
            // detect frame (look for p/n sequence)
            execute_seekpn(_x, i);
            break;
          case STATE_RXPN:
            // detect frame (look for p/n sequence)
            execute_rxpn(_x, i);
            break;
          case STATE_RXPAYLOAD:
            // detect frame (look for p/n sequence)
            execute_rxpayload(_x, i);
            break;
        default:
            fprintf(stderr,"error: framesync64_exeucte(), unknown/unsupported state\n");
            exit(1);
        }
      }
    }

    void framesync::update_training_sequence_n() {
      while(training_seq_n[0] < training_seq_n[1])
        training_seq_n[1]--;
      while(training_seq_n[1] < training_seq_n[0])
        training_seq_n[0]--;
    }

    void framesync::training_match_index_diffs()
    {
      ;
    }

    void framesync::execute_seekpn(std::complex<float> ** _x, unsigned int i)
    {
      nco_crcf_mix_down(nco_coarse[0], _x[0][i], &_x[0][i]);
      nco_crcf_step(nco_coarse[0]);
      nco_crcf_mix_down(nco_coarse[1], _x[1][i], &_x[1][i]);
      nco_crcf_step(nco_coarse[1]);

      windowcf_push(buffer[0], _x[0][i]);
      windowcf_push(buffer[1], _x[1][i]);

      int ch1_flag = detector_cccf_correlate(frame_detector[0],
                                             _x[0][i],
                                             &tau_hat[0],
                                             &dphi_hat[0],
                                             &gamma_hat[0]);
      int ch2_flag = detector_cccf_correlate(frame_detector[1],
                                             _x[1][i],
                                             &tau_hat[1],
                                             &dphi_hat[1],
                                             &gamma_hat[1]);
 
      if(ch1_flag) {
        frame_detect_flag[0] = true;
        training_match_index[0] = i;
        training_seq_n[0]++;
      }
      if(ch2_flag) {
        frame_detect_flag[1] = true;
        training_match_index[1] = i;
        training_seq_n[1]++;
      }

      if(frame_detect_flag[0] and frame_detect_flag[1]) {
        if (training_seq_n[0] != training_seq_n[1]) {
          reset();
          return;
        }

        float delta_f = (dphi_hat[0] + dphi_hat[1])/2.0;
        if(fabs(delta_f) > 0.001) {
          nco_coarse_freq += delta_f;
          nco_crcf_set_frequency(nco_coarse[0], nco_coarse_freq);
          nco_crcf_set_frequency(nco_coarse[1], nco_coarse_freq);
        }
        if(PRINT_STATS_TX1RX1) {
          printf("***** TX1 RX1 Statistics! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB, index:%u, seq_n:%lu\n",
                 tau_hat[0], dphi_hat[0], 20*log10f(gamma_hat[0]), training_match_index[0],
                 training_seq_n[0]);
        }
        if(PRINT_STATS_TX1RX2) {
          printf("***** TX1 RX2 Statistics! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB, index:%u, seq_n:%lu\n",
                 tau_hat[1], dphi_hat[1], 20*log10f(gamma_hat[1]), training_match_index[1],
                 training_seq_n[1]);
        }

        if(PRINT_STAT_DIFFS) {
          printf("***** RX1 RX2 Defferences! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB\n",
                 tau_hat[0]-tau_hat[1], 
                 dphi_hat[0]-dphi_hat[1], 
                 20*log10f(gamma_hat[0]/gamma_hat[1]));
        }

        pushpn();
      }
    }

    void framesync::pushpn()
    {
      int delay = find_index_of_corr();
      reset();
    }

    int framesync::find_index_of_corr() {
      return 0;
    }

    void framesync::execute_rxpn(std::complex<float> ** _x, unsigned int i)
    {
    }

    int framesync::update_symsync(std::complex<float>   _x,
                           std::complex<float> * _y)
    {
      return 0;
    }

    void framesync::syncpn()
    {
      ;
    }

    void framesync::execute_rxpayload(std::complex<float> ** _x, unsigned int i)
    {
      ;
    }

  }   // namespace mimo
}     // namespace liquid
