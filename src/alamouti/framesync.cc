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

#define PRINT_STATS_TX1RX1  true 
#define PRINT_STATS_TX1RX2  true 
#define PRINT_STATS_TX2RX1  true 
#define PRINT_STATS_TX2RX2  true 
#define PRINT_STAT_DIFFS    false
#define PRINT_PHASES        false
#define SHOW_OTHER_DOTS     true
#define PRINT_PAYLOAD       false
#define USE_AVERAGE_METRIC  false

namespace liquid {
  namespace alamouti {
    framesync::framesync(unsigned int _k,
                         unsigned int _m,
                         float _beta)
    {
      training_seq_len = 63;
      payload_len = 1024;

      msequence ms1 = msequence_create(6, 0x005B, 1);
      msequence ms2 = msequence_create(6, 0x0067, 1);
      msequence ms3 = msequence_create(10, 0x0409, 1);

      preamble_pn[0] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len); 
      preamble_pn[1] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len); 
      preamble_rx[0][0] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len); 
      preamble_rx[0][1] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len); 
      preamble_rx[1][0] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len); 
      preamble_rx[1][1] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len); 

      expected_payload = (unsigned char *)malloc(sizeof(unsigned char)*payload_len);
      for (unsigned int i = 0; i < training_seq_len; i++){
        preamble_pn[0][i] = (msequence_advance(ms1)) ? 1.0f : -1.0f;
        preamble_pn[1][i] = (msequence_advance(ms2)) ? 1.0f : -1.0f;
      }

      for (unsigned int i = 0; i < payload_len; i++) {
        expected_payload[i] = msequence_advance(ms3);
      }

      msequence_destroy(ms1);
      msequence_destroy(ms2);

      // interpolate p/n sequence with matched filter
      k     = _k;        // samples/symbol
      m     = _m;        // filter delay (symbols)
      beta  = _beta;    // excess bandwidth factor
      frame_len = k*(3*training_seq_len + payload_len + 3*m);

      std::complex<float> seq1[k*training_seq_len];
      std::complex<float> seq2[k*training_seq_len];
      firinterp_crcf interp1 = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);
      firinterp_crcf interp2 = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);

/*      for (unsigned int i = 0; i < training_seq_len; i++) {
        firinterp_crcf_execute(interp1, preamble_pn[0][i%training_seq_len], &seq1[k*i]);
        firinterp_crcf_execute(interp2, preamble_pn[1][i%training_seq_len], &seq2[k*i]);
      }*/     
      for (unsigned int i = 0; i < training_seq_len + m; i++) {
        // compensate for filter delay
        if (i < m) {
            firinterp_crcf_execute(interp1, preamble_pn[0][i], &seq1[0]);
            firinterp_crcf_execute(interp2, preamble_pn[1][i], &seq2[0]);
        }
        else {
            firinterp_crcf_execute(interp1, preamble_pn[0][i%training_seq_len], &seq1[k*(i - m)]);
            firinterp_crcf_execute(interp2, preamble_pn[1][i%training_seq_len], &seq2[k*(i - m)]);
        }
      }
      firinterp_crcf_destroy(interp1);
      firinterp_crcf_destroy(interp2);

      // create frame detector
      float threshold = 0.5f;     // detection threshold
      float dphi_max  = 0.05f;    // maximum carrier offset allowable
      frame_detector[0][0] = detector_cccf_create(seq1, k*training_seq_len, threshold, dphi_max);
      frame_detector[0][1] = detector_cccf_create(seq1, k*training_seq_len, threshold, dphi_max);
      frame_detector[1][0] = detector_cccf_create(seq2, k*training_seq_len, threshold, dphi_max);
      frame_detector[1][1] = detector_cccf_create(seq2, k*training_seq_len, threshold, dphi_max);

      rx_payload[0] = (unsigned char *)malloc(sizeof(unsigned char)*payload_len);
      rx_payload[1] = (unsigned char *)malloc(sizeof(unsigned char)*payload_len);
      demod[0] = modem_create(LIQUID_MODEM_BPSK);
      demod[1] = modem_create(LIQUID_MODEM_BPSK);
      payload_window[0] = windowcf_create(payload_len);
      payload_window[1] = windowcf_create(payload_len);
      pn_window[0] = windowcf_create(3*training_seq_len + 2*m);
      pn_window[1] = windowcf_create(3*training_seq_len + 2*m);
      pn_window_unscaled[0] = windowcf_create(3*training_seq_len + 2*m);
      pn_window_unscaled[1] = windowcf_create(3*training_seq_len + 2*m);
      pn_dotprods[0] = dotprod_cccf_create(preamble_pn[0], training_seq_len);
      pn_dotprods[1] = dotprod_cccf_create(preamble_pn[1], training_seq_len);

      rx_sig[0] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*frame_len*2);
      rx_sig[1] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*frame_len*2);
      last_rx_sig[0] = rx_sig[0];
      last_rx_sig[1] = rx_sig[1];
      curr_rx_sig[0] = rx_sig[0] + frame_len;
      curr_rx_sig[1] = rx_sig[1] + frame_len;

      npfb = 32;
      mf[0]   = firpfb_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER, npfb,k,m,beta);
      dmf[0]  = firpfb_crcf_create_drnyquist(LIQUID_FIRFILT_ARKAISER,npfb,k,m,beta);
      mf[1]   = firpfb_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER, npfb,k,m,beta);
      dmf[1]  = firpfb_crcf_create_drnyquist(LIQUID_FIRFILT_ARKAISER,npfb,k,m,beta);

      nco_coarse_freq = 0.0f;
      nco_fine_freq[0] = 0.0f;
      nco_fine_freq[1] = 0.0f;
      nco_fine_phase[0] = 0.0f;
      nco_fine_phase[1] = 0.0f;
      nco_coarse[0] = nco_crcf_create(LIQUID_NCO);
      nco_coarse[1] = nco_crcf_create(LIQUID_NCO);
      nco_fine[0]   = nco_crcf_create(LIQUID_VCO);
      nco_fine[1]   = nco_crcf_create(LIQUID_VCO);
      nco_crcf_pll_set_bandwidth(nco_fine[0], 0.05f);
      nco_crcf_pll_set_bandwidth(nco_fine[1], 0.05f);

      dphi_hat[0][0] = 0;
      dphi_hat[0][1] = 0;
      dphi_hat[1][0] = 0;
      dphi_hat[1][1] = 0;

      nco_crcf_reset(nco_coarse[0]);
      nco_crcf_reset(nco_coarse[1]);
      nco_crcf_reset(nco_fine[0]);
      nco_crcf_reset(nco_fine[1]);

      // reset symbol timing recovery state
      firpfb_crcf_reset(mf[0]);
      firpfb_crcf_reset(mf[1]);
      firpfb_crcf_reset(dmf[0]);
      firpfb_crcf_reset(dmf[1]);
      pfb_timer = k;
      pfb_soft = 0.0;
      pfb_index = 0;

      training_seq_n[0][0] = 0;
      training_seq_n[0][1] = 0;
      training_seq_n[1][0] = 0;
      training_seq_n[1][1] = 0;

      estimator = new channel_estimator(seq1, k*training_seq_len, threshold, dphi_max);
      work_exec_cycle = 0;

      out_files[0] = fopen("/tmp/payload_data1.txt", "w");
      out_files[1] = fopen("/tmp/payload_data2.txt", "w");

      num_frames_detected = 0;
      errors[0] = 0;
      errors[1] = 0;

      reset();
    }

    void framesync::print_frame_struct()
    {
      printf("*****                     Frame Structure                             *****\n");
      printf("Duration : || 63 Symb |  3 Symb | 63 Symb | 63 Symb | 1024 Symb | 3 Symb ||\n");
      printf("Channel 1: || PN1 Seq |  Zeros  |  Zeros  | Phasing |  Payload  |  Zeros ||\n");
      printf("Channel 2: ||  Zeros  |  Zeros  | PN2 Seq | Phasing |  Payload  |  Zeros ||\n");
      printf("Pulse Shape : k = %u, m = %u, beta = %f\n", k, m, beta);
      printf("Frame length : %u\n", frame_len);
      printf("***************************************************************************\n");
    }

    void framesync::reset()
    {
      frame_detect_flag[0][0] = false;
      frame_detect_flag[0][1] = false;
      frame_detect_flag[1][0] = false;
      frame_detect_flag[1][1] = false;
      state = STATE_DETECTFRAME1;

      detector_cccf_reset(frame_detector[0][0]);
      detector_cccf_reset(frame_detector[0][1]);
      detector_cccf_reset(frame_detector[1][0]);
      detector_cccf_reset(frame_detector[1][1]);

      windowcf_clear(pn_window[0]);
      windowcf_clear(pn_window[1]);
      windowcf_clear(pn_window_unscaled[0]);
      windowcf_clear(pn_window_unscaled[1]);
      windowcf_clear(payload_window[0]);
      windowcf_clear(payload_window[1]);

      pn_count = 0;
      payload_count = 0;
//      update_training_sequence_n();
    }

    framesync::~framesync()
    {
      free(preamble_pn[0]);
      free(preamble_pn[1]);
      free(preamble_rx[0][0]);
      free(preamble_rx[0][1]);
      free(preamble_rx[1][0]);
      free(preamble_rx[1][1]);
      free(rx_sig[0]);
      free(rx_sig[1]);
      free(rx_payload[0]);
      free(rx_payload[1]);
      free(expected_payload);
      delete estimator;

      detector_cccf_destroy(frame_detector[0][0]);
      detector_cccf_destroy(frame_detector[0][1]);
      detector_cccf_destroy(frame_detector[1][0]);
      detector_cccf_destroy(frame_detector[1][1]);

      windowcf_destroy(pn_window[0]);
      windowcf_destroy(pn_window[1]);
      windowcf_destroy(payload_window[0]);
      windowcf_destroy(payload_window[1]);

      firpfb_crcf_destroy(mf[0]);
      firpfb_crcf_destroy(mf[1]);
      firpfb_crcf_destroy(dmf[0]);
      firpfb_crcf_destroy(dmf[1]);

      nco_crcf_destroy(nco_coarse[0]);
      nco_crcf_destroy(nco_coarse[1]);
      nco_crcf_destroy(nco_fine[0]);
      nco_crcf_destroy(nco_fine[1]);

      modem_destroy(demod[0]);
      modem_destroy(demod[1]);

      dotprod_cccf_destroy(pn_dotprods[0]);
      dotprod_cccf_destroy(pn_dotprods[1]);

      fclose(out_files[0]);
      fclose(out_files[1]);
    }

    unsigned int framesync::get_frame_len() {
      return frame_len;
    }

    void framesync::work(std::complex<float> ** _x,
                         unsigned int    _n)
    {
      assert(_n == frame_len);
      rx_sig_index = 0;
      for (unsigned int i = 0; i < _n; i++) {
        switch (state) {
          case STATE_DETECTFRAME1:
            // detect frame (look for p/n sequence)
            execute_seekpn1(_x, i);
            break;
          case STATE_DETECTFRAME2:
            // detect frame (look for p/n sequence)
            execute_seekpn2(_x, i);
            break;
          case STATE_RXPN:
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
      memmove(last_rx_sig[0], curr_rx_sig[0], sizeof(std::complex<float>)*frame_len);
      memmove(last_rx_sig[1], curr_rx_sig[1], sizeof(std::complex<float>)*frame_len);
      printf("Synchronizer Execution cycle :%u\n",work_exec_cycle++);
      training_match_index_diffs1();
      training_match_index_diffs2();
    }

    void framesync::execute_seekpn1(std::complex<float> ** _x, unsigned int i)
    {
      std::complex<float> y0, y1;

      nco_crcf_mix_down(nco_coarse[0], _x[0][i], &y0);
      nco_crcf_step(nco_coarse[0]);
      nco_crcf_mix_down(nco_coarse[1], _x[1][i], &y1);
      nco_crcf_step(nco_coarse[1]);

      curr_rx_sig[0][rx_sig_index] = y0;
      curr_rx_sig[1][rx_sig_index] = y1;
      rx_sig_index++;

      if(PRINT_PHASES) {
        float phase_1 = std::arg(y0);
        float phase_2 = std::arg(y1);
        printf("Channel 1 Phase:%8.4f, Channel 2 Phase:%8.4f, Phase Diffs:%8.4f\n",
            phase_1, phase_2, phase_2 - phase_1);
      }

      int ch1_flag = detector_cccf_correlate(frame_detector[0][0],
                                             y0,
                                             &tau_hat[0][0],
                                             &dphi_hat[0][0],
                                             &gamma_hat[0][0]);
      int ch2_flag = detector_cccf_correlate(frame_detector[0][1],
                                             y1,
                                             &tau_hat[0][1],
                                             &dphi_hat[0][1],
                                             &gamma_hat[0][1]);
 
      if(ch1_flag) {
        frame_detect_flag[0][0] = true;
        training_match_index[0][0] = i;
      }
      if(ch2_flag) {
        frame_detect_flag[0][1] = true;
        training_match_index[0][1] = i;
      }

      if(frame_detect_flag[0][0] and frame_detect_flag[0][1]) {
        float delta_f = (dphi_hat[0][0] + dphi_hat[0][1])/2.0;
        if(fabs(delta_f) > 0.001) {
          nco_coarse_freq += delta_f;
          nco_crcf_set_frequency(nco_coarse[0], nco_coarse_freq);
          nco_crcf_set_frequency(nco_coarse[1], nco_coarse_freq);
        }

        if(PRINT_STATS_TX1RX1) {
          printf("***** TX1 RX1 Statistics! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB, index:%u, seq_n:%lu\n",
                 tau_hat[0][0], dphi_hat[0][0], 20*log10f(gamma_hat[0][0]), training_match_index[0][0],
                 training_seq_n[0][0]++);
        }
        if(PRINT_STATS_TX1RX2) {
          printf("***** TX1 RX2 Statistics! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB, index:%u, seq_n:%lu\n",
                 tau_hat[0][1], dphi_hat[0][1], 20*log10f(gamma_hat[0][1]), training_match_index[0][1],
                 training_seq_n[0][1]++);
        }
        if(PRINT_STAT_DIFFS) {
          printf("***** RX1 RX2 Defferences! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB\n",
                 tau_hat[0][0]-tau_hat[0][1], 
                 dphi_hat[0][0]-dphi_hat[0][1], 
                 20*log10f(gamma_hat[0][0]/gamma_hat[0][1]));
        }

        wait_for_pn2 = 0;
        state = STATE_DETECTFRAME2;
      }
    }

    void framesync::execute_seekpn2(std::complex<float> ** _x, unsigned int i)
    {
      if(wait_for_pn2 > 4*training_seq_len*k)
      {
        printf("Well, I think I dropped a frame in between, Going to reset\n");
        reset();
      }
      wait_for_pn2++;
      std::complex<float> y0, y1;

      nco_crcf_mix_down(nco_coarse[0], _x[0][i], &y0);
      nco_crcf_step(nco_coarse[0]);
      nco_crcf_mix_down(nco_coarse[1], _x[1][i], &y1);
      nco_crcf_step(nco_coarse[1]);

      curr_rx_sig[0][rx_sig_index] = y0;
      curr_rx_sig[1][rx_sig_index] = y1;
      rx_sig_index++;

      if(PRINT_PHASES) {
        float phase_1 = std::arg(y0);
        float phase_2 = std::arg(y1);
        printf("Channel 1 Phase:%8.4f, Channel 2 Phase:%8.4f, Phase Diffs:%8.4f\n",
            phase_1, phase_2, phase_2 - phase_1);
      }

      int ch1_flag = detector_cccf_correlate(frame_detector[1][0],
                                             y0,
                                             &tau_hat[1][0],
                                             &dphi_hat[1][0],
                                             &gamma_hat[1][0]);
      int ch2_flag = detector_cccf_correlate(frame_detector[1][1],
                                             y1,
                                             &tau_hat[1][1],
                                             &dphi_hat[1][1],
                                             &gamma_hat[1][1]);
 
      if(ch1_flag) {
        frame_detect_flag[1][0] = true;
        training_match_index[1][0] = i;
      }
      if(ch2_flag) {
        frame_detect_flag[1][1] = true;
        training_match_index[1][1] = i;
      }

      if(frame_detect_flag[1][0] and frame_detect_flag[1][1]) {

        float delta_f = (dphi_hat[0][0] + dphi_hat[0][1])/2.0;
        if(fabs(delta_f) > 0.001) {
          nco_coarse_freq += (dphi_hat[1][0] + dphi_hat[1][1])/2.0;
          nco_crcf_set_frequency(nco_coarse[0], nco_coarse_freq);
          nco_crcf_set_frequency(nco_coarse[1], nco_coarse_freq);
        }
        
        if(PRINT_STATS_TX2RX1) {
          printf("***** TX2 RX1 Statistics! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB, index:%u, seq_n:%lu\n",
                 tau_hat[1][0], dphi_hat[1][0], 20*log10f(gamma_hat[1][0]), training_match_index[1][0],
                 training_seq_n[1][0]++);
        }
        if(PRINT_STATS_TX2RX2) {
          printf("***** TX2 RX2 Statistics! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB, index:%u, seq_n:%lu\n",
                 tau_hat[1][1], dphi_hat[1][1], 20*log10f(gamma_hat[1][1]), training_match_index[1][1],
                 training_seq_n[1][1]++);
        }
        if(PRINT_STAT_DIFFS) {
          printf("***** RX1 RX2 Defferences! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB\n",
                 tau_hat[1][0]-tau_hat[1][1], 
                 dphi_hat[1][0]-dphi_hat[1][1], 
                 20*log10f(gamma_hat[1][0]/gamma_hat[1][1]));
        }

        pushpn();
      }
    }

    // reads the already stored data and pushes through the matche filter
    // should check if we should proceed to execute_rxpn / sync_pn.
    void framesync::pushpn()
    {
      std::complex<float> mf_in[2];
      std::complex<float> mf_out[2];
      int start; 

      printf("training_match_index[0][0] :%6u\n", training_match_index[0][0]);
      printf("training_match_index[0][1] :%6u\n", training_match_index[0][1]);
      printf("training_match_index[1][0] :%6u\n", training_match_index[1][0]);
      printf("training_match_index[1][1] :%6u\n", training_match_index[1][1]);
      printf("              rx_sig_index :%6u\n", rx_sig_index);

      if(training_match_index[0][0] > training_match_index[1][0] and
         training_match_index[0][1] > training_match_index[1][1])
      {
        start = -training_seq_len*k - frame_len +
                ((training_match_index[0][0] < training_match_index[0][1])?
                  training_match_index[0][0] : training_match_index[0][1]);
      }
      else if(training_match_index[0][0] < training_match_index[1][0] and
              training_match_index[0][1] < training_match_index[1][1])
      {
        start = -training_seq_len*k +
                ((training_match_index[0][0] < training_match_index[0][1])?
                  training_match_index[0][0] : training_match_index[0][1]);
      }
      else {
        printf("***** I smell somthing fishy\n");
        assert(false);
      }
      printf("               start_index :%6u\n", start);

      for(int i = start; i < rx_sig_index; i++)
      {
        mf_in[0] = curr_rx_sig[0][i];
        mf_in[1] = curr_rx_sig[1][i];

        if(update_symsync(mf_in, mf_out))
        {
          windowcf_push(pn_window_unscaled[0], mf_out[0]);
          windowcf_push(pn_window_unscaled[1], mf_out[1]);
          windowcf_push(pn_window[0], std::polar(1.0f, std::arg(mf_out[0])));
          windowcf_push(pn_window[1], std::polar(1.0f, std::arg(mf_out[1])));
          pn_count++;
          if(pn_count == 3*training_seq_len + 2*m)
            sync_pn();
        }
      }
      if(state == STATE_DETECTFRAME2)
      {
        printf("Setting state to RXPN\n");
        state = STATE_RXPN;
      }
    }

    // this checks if we have actually received the pn sequence and the
    // phasing sequence. This is fired when pn_count reaches full. It should 
    // check if the values stored in the pn_window actually corresponds to 
    // pn sequences. If not decrement the pn_count. If they match, then 
    // proceed to rx_payload.
    // This should find out the nco_fine offset and set the values accordingly.
    void framesync::sync_pn()
    {
      std::complex<float> dotp[2][2];
      std::complex<float> * pn_window_ptr[2];
      std::complex<float> * pn_window_unscaled_ptr[2];
      bool sync;
      windowcf_read(pn_window[0], &pn_window_ptr[0]);
      windowcf_read(pn_window[1], &pn_window_ptr[1]);
      windowcf_read(pn_window_unscaled[0], &pn_window_unscaled_ptr[0]);
      windowcf_read(pn_window_unscaled[1], &pn_window_unscaled_ptr[1]);

      dotprod_cccf_execute(pn_dotprods[0], pn_window_ptr[0], &dotp[0][0]);
      dotprod_cccf_execute(pn_dotprods[0], pn_window_ptr[1], &dotp[0][1]);
      dotprod_cccf_execute(pn_dotprods[1], pn_window_ptr[0] + 66, &dotp[1][0]);
      dotprod_cccf_execute(pn_dotprods[1], pn_window_ptr[1] + 66, &dotp[1][1]);

      sync = (std::abs(dotp[0][1]) > 50.0);

      if(sync){

        printf("dot products [0][0] :%8.4f + i%8.4f, mag :%8.4f\n",
            std::real(dotp[0][0]), std::imag(dotp[0][0]), std::abs(dotp[0][0]));
        printf("dot products [0][1] :%8.4f + i%8.4f, mag :%8.4f\n",
            std::real(dotp[0][1]), std::imag(dotp[0][1]), std::abs(dotp[0][1]));
        printf("dot products [1][0] :%8.4f + i%8.4f, mag :%8.4f\n",
            std::real(dotp[1][0]), std::imag(dotp[1][0]), std::abs(dotp[1][0]));
        printf("dot products [1][1] :%8.4f + i%8.4f, mag :%8.4f\n",
            std::real(dotp[1][1]), std::imag(dotp[1][1]), std::abs(dotp[1][1]));

        std::complex<float> phasing_symb;
        // #FIXME
        {
          std::complex<float> dphi_metric[2] = {0.0f, 0.0f};
          std::complex<float> theta_metric[2] = {0.0f, 0.0f};
          float dphi_hat[2];
          float theta_hat[2];
          std::complex<float> r0[2] = {0.0f, 0.0f};
          std::complex<float> r1[2] = {0.0f, 0.0f};
  
          for(unsigned int i = 0; i < training_seq_len; i++)
          {
            phasing_symb = (i % 2) ? 1.0f : -1.0f;
            r0[0] = r1[0];
            r0[1] = r1[1];
            r1[0] = pn_window_ptr[0][2*(m + training_seq_len) + i]*phasing_symb;
            r1[0] = pn_window_ptr[1][2*(m + training_seq_len) + i]*phasing_symb;
            dphi_metric[0] += r1[0]*std::conj(r0[0]);
            dphi_metric[1] += r1[1]*std::conj(r0[1]);
          }

          if(USE_AVERAGE_METRIC) {
            dphi_hat[0] = std::arg((dphi_metric[0] + dphi_metric[1])/2.0f);
            dphi_hat[1] = std::arg((dphi_metric[0] + dphi_metric[1])/2.0f);
          }
          else {
            dphi_hat[0] = std::arg(dphi_metric[0]);
            dphi_hat[1] = std::arg(dphi_metric[1]);
          }

          for(unsigned int i = 0; i < training_seq_len; i++)
          {
            phasing_symb = (i % 2) ? 1.0f : -1.0f;
            theta_metric[0] += (*(pn_window_ptr[0] + 132 + i))*
              std::exp(-dphi_hat[0]*i)*phasing_symb;
            theta_metric[1] += (*(pn_window_ptr[1] + 132 + i))*
              std::exp(-dphi_hat[1]*i)*phasing_symb;
          }

          if(USE_AVERAGE_METRIC) {
            theta_hat[0] = std::arg((theta_metric[0] + theta_metric[1])/2.0f);
            theta_hat[1] = std::arg((theta_metric[0] + theta_metric[1])/2.0f);
          }
          else {
            theta_hat[0] = std::arg(theta_metric[0]);
            theta_hat[1] = std::arg(theta_metric[1]);
          }

          printf("dphi_metric[0] :%8.4f + %8.4fi\n",
                 std::real(dphi_metric[0]), std::imag(dphi_metric[0]));
  
          printf("dphi_metric[1] :%8.4f + %8.4fi\n",
                 std::real(dphi_metric[1]), std::imag(dphi_metric[1]));
  
          printf("dphi_hat[0] :%8.4f\n", dphi_hat[0]);
          printf("dphi_hat[1] :%8.4f\n", dphi_hat[1]);
          
          printf("theta_metric[0] :%8.4f + %8.4fi\n",
                 std::real(theta_metric[0]), std::imag(theta_metric[0]));
          
          printf("theta_metric[1] :%8.4f + %8.4fi\n",
                 std::real(theta_metric[1]), std::imag(theta_metric[1]));
          
          printf("theta_hat[0] :%8.4f\n", theta_hat[0]);
          printf("theta_hat[1] :%8.4f\n", theta_hat[1]);

          // initialize fine-tuned nco
          nco_fine_freq[0] = dphi_hat[0];
          nco_fine_freq[1] = dphi_hat[1];
          nco_fine_phase[0] = theta_hat[0];
          nco_fine_phase[1] = theta_hat[1];

          nco_crcf_set_frequency(nco_fine[0], nco_fine_freq[0]);
          nco_crcf_set_phase(nco_fine[0], nco_fine_phase[0]);
          nco_crcf_set_frequency(nco_fine[1], nco_fine_freq[1]);
          nco_crcf_set_phase(nco_fine[1], nco_fine_phase[1]);

          std::complex<float> temp[2];
          float sumsq[2] = {0.0f, 0.0f};
          float phase_error[2];
          for(unsigned int i = 0; i < training_seq_len; i++)
          {
            phasing_symb = (i % 2) ? 1.0f : -1.0f;
            nco_crcf_mix_down(nco_fine[0], pn_window_unscaled_ptr[0][132 + i], &temp[0]);
            nco_crcf_mix_down(nco_fine[1], pn_window_unscaled_ptr[1][132 + i], &temp[1]);

            phase_error[0] = std::imag((temp[0]*phasing_symb + temp[1]*phasing_symb)/2.0f);
            phase_error[1] = std::imag((temp[0]*phasing_symb + temp[1]*phasing_symb)/2.0f);

            nco_crcf_pll_step(nco_fine[0], phase_error[0]);
            nco_crcf_pll_step(nco_fine[1], phase_error[1]);

            nco_crcf_step(nco_fine[0]);
            nco_crcf_step(nco_fine[1]);

            sumsq[0] += std::real(temp[0] * std::conj(temp[0]));
            sumsq[1] += std::real(temp[1] * std::conj(temp[1]));
          }
          pow[0] = sqrtf(sumsq[0]/(float)training_seq_len);
          pow[1] = sqrtf(sumsq[1]/(float)training_seq_len);
          printf("pow[0] : %8.4f dB\n", 20*log10f(pow[0]));
          printf("pow[1] : %8.4f dB\n", 20*log10f(pow[1]));
        }

        state = STATE_RXPAYLOAD;
      }
      else
      {
        pn_count--;
      }
    }

    // this should do the maximum ratio combining of the values stored in the
    // payload buffer.
    void framesync::decode_payload()
    {
      std::complex<float> * payload_window_ptr[2];
      std::complex<float> symbols[2];
      unsigned int sym_out[2] = {0, 0};
      windowcf_read(payload_window[0], &payload_window_ptr[0]);
      windowcf_read(payload_window[1], &payload_window_ptr[1]);

      for(unsigned int i = 0; i < payload_len; i++)
      {
        nco_crcf_mix_down(nco_fine[0], payload_window_ptr[0][i]*0.5f/pow[0], &symbols[0]);
        nco_crcf_mix_down(nco_fine[1], payload_window_ptr[1][i]*0.5f/pow[1], &symbols[1]);

        modem_demodulate(demod[0], symbols[0], &sym_out[0]);
        modem_demodulate(demod[1], symbols[1], &sym_out[1]);

        float phase_error = (modem_get_demodulator_phase_error(demod[0]) +
                             modem_get_demodulator_phase_error(demod[1]))/2.0f;

        rx_payload[0][i] = sym_out[0];
        rx_payload[1][i] = sym_out[1];

        if(sym_out[0] != expected_payload[i])
          errors[0]++;
        if(sym_out[1] != expected_payload[i])
          errors[1]++;

        nco_crcf_pll_step(nco_fine[0], phase_error);
        nco_crcf_pll_step(nco_fine[1], phase_error);

        nco_crcf_step(nco_fine[0]);
        nco_crcf_step(nco_fine[1]);
      }

      if(PRINT_PAYLOAD) {
        for(unsigned int i = 0; i < payload_len; i++)
        {
          printf("%d", rx_payload[0][i]);
        }
        printf("\n");
        for(unsigned int i = 0; i < payload_len; i++)
        {
          printf("%d", rx_payload[1][i]);
        }
        printf("\n");
      }
      num_frames_detected++;
      reset();
    }

    void framesync::execute_rxpn(std::complex<float> ** _x, unsigned int i)
    {
      std::complex<float> mf_in[2];
      std::complex<float> mf_out[2];

      nco_crcf_mix_down(nco_coarse[0], _x[0][i], &mf_in[0]);
      nco_crcf_step(nco_coarse[0]);
      nco_crcf_mix_down(nco_coarse[1], _x[1][i], &mf_in[1]);
      nco_crcf_step(nco_coarse[1]);

      curr_rx_sig[0][rx_sig_index] = mf_in[0];
      curr_rx_sig[1][rx_sig_index] = mf_in[1];
      rx_sig_index++;

      if(update_symsync(mf_in, mf_out))
      {
        windowcf_push(pn_window_unscaled[0], mf_out[0]);
        windowcf_push(pn_window_unscaled[1], mf_out[1]);
        windowcf_push(pn_window[0], std::polar(1.0f, std::arg(mf_out[0])));
        windowcf_push(pn_window[1], std::polar(1.0f, std::arg(mf_out[1])));
        pn_count++;
        if(pn_count == 3*training_seq_len + 2*m)
          sync_pn();
      }
    }

    unsigned long int framesync::get_num_bits_detected()
    {
      return(num_frames_detected*payload_len);
    }

    unsigned long int framesync::get_num_errors1()
    {
      return(errors[0]);
    }

    unsigned long int framesync::get_num_errors2()
    {
      return(errors[1]);
    }

    void framesync::execute_rxpayload(std::complex<float> ** _x, unsigned int i)
    {
      std::complex<float> mf_in[2];
      std::complex<float> mf_out[2];

      nco_crcf_mix_down(nco_coarse[0], _x[0][i], &mf_in[0]);
      nco_crcf_step(nco_coarse[0]);
      nco_crcf_mix_down(nco_coarse[1], _x[1][i], &mf_in[1]);
      nco_crcf_step(nco_coarse[1]);

      curr_rx_sig[0][rx_sig_index] = mf_in[0];
      curr_rx_sig[1][rx_sig_index] = mf_in[1];
      rx_sig_index++;

      if(update_symsync(mf_in, mf_out))
      {
        windowcf_push(payload_window[0], mf_out[0]);
        windowcf_push(payload_window[1], mf_out[1]);
        payload_count++;
        if(payload_count == payload_len)
          decode_payload();
      }
    }

    int framesync::update_symsync(std::complex<float> * in,
                                  std::complex<float> * out)
    {
      firpfb_crcf_push(mf[0], in[0]);
      firpfb_crcf_push(dmf[0], in[0]);
      firpfb_crcf_push(mf[1], in[1]);
      firpfb_crcf_push(dmf[1], in[1]);

      std::complex<float> mf_out[2];
      std::complex<float> dmf_out[2];

      int sample_available = 0;

      if (pfb_timer <= 0) {
        sample_available = 1;
        pfb_timer = k;
        firpfb_crcf_execute(mf[0], pfb_index, &mf_out[0]);
        firpfb_crcf_execute(mf[1], pfb_index, &mf_out[1]);
        firpfb_crcf_execute(dmf[0], pfb_index, &dmf_out[0]);
        firpfb_crcf_execute(dmf[1], pfb_index, &dmf_out[1]);

        // update filtered timing error
        // hi  bandwidth parameters: {0.92, 1.20}, about 100 symbols settling time
        // med bandwidth parameters: {0.98, 0.20}, about 200 symbols settling time
        // lo  bandwidth parameters: {0.99, 0.05}, about 500 symbols settling time
        pfb_q = 0.99f*pfb_q + 0.025f*std::real(std::conj(mf_out[0])*dmf_out[0] + 
                                               std::conj(mf_out[1])*dmf_out[1]);

        pfb_soft += pfb_q;
        pfb_index = roundf(pfb_soft);

        while(pfb_index < 0) {
          pfb_index += npfb;
          pfb_soft  += npfb;
          pfb_timer--;
        }
        while(pfb_index > (npfb - 1)) {
          pfb_index -= npfb;
          pfb_soft  -= npfb;
          pfb_timer++;
        }
      }
      pfb_timer--;

      out[0] = mf_out[0];
      out[1] = mf_out[1];
      return sample_available;
    }

    void framesync::find_channel_diffs() {
      printf("Estimated index for channel 1 :%u\n",
             training_match_index[0][0] + 
             estimator->find_corr_index(curr_rx_sig[0], training_match_index[0][0]));
      printf("Estimated index for channel 2 :%u\n",
             training_match_index[0][1] + 
             estimator->find_corr_index(curr_rx_sig[1], training_match_index[0][1]));
    }

    void framesync::update_training_sequence_n() {
      unsigned long int * seq_n_ptr = (unsigned long int *)&training_seq_n;
      unsigned long int min = ULONG_MAX;
      for(unsigned int i = 0; i < 4; i++) {
        min = ((seq_n_ptr[i] < min) ? seq_n_ptr[i] : min);
      }
      for(unsigned int i = 0; i < 4; i++) {
        seq_n_ptr[i] = min;
      }
    }

    void framesync::training_match_index_diffs1()
    {
      if(training_match_index[0][0] > training_match_index[1][0] and
         training_match_index[0][1] > training_match_index[1][1])
      {
        printf("***** Index Diffs! Channel 1:%u, Channel 2:%u\n",
               frame_len - training_match_index[0][0] + training_match_index[1][0],
               frame_len - training_match_index[0][1] + training_match_index[1][1]);
      }
      else if(training_match_index[0][0] < training_match_index[1][0] and
              training_match_index[0][1] < training_match_index[1][1])
      {
        printf("***** Index Diffs! Channel 1:%u, Channel 2:%u\n",
               0 - training_match_index[0][0] + training_match_index[1][0],
               0 - training_match_index[0][1] + training_match_index[1][1]);
      }
      else {
        printf("***** I smell somthing fishy\n");
      }
    }

    void framesync::training_match_index_diffs2()
    {
      if(training_match_index[1][0] > training_match_index[0][0] and
         training_match_index[1][1] > training_match_index[0][1])
      {
        printf("***** Index Diffs! Channel 1:%u, Channel 2:%u\n",
               frame_len - training_match_index[1][0] + training_match_index[0][0],
               frame_len - training_match_index[1][1] + training_match_index[0][1]);
      }
      else if(training_match_index[1][0] < training_match_index[0][0] and
              training_match_index[1][1] < training_match_index[0][1])
      {
        printf("***** Index Diffs! Channel 1:%u, Channel 2:%u\n",
               0 - training_match_index[1][0] + training_match_index[0][0],
               0 - training_match_index[1][1] + training_match_index[0][1]);
      }
      else {
        printf("***** I smell somthing fishy\n");
        assert(false);
      }
    }
  }   // namespace mimo
}     // namespace liquid

/*
        if(SHOW_OTHER_DOTS) {
          dotprod_cccf_execute(pn_dotprods[0], pn_window_ptr[0] + 1, &dotp[0][0]);
          dotprod_cccf_execute(pn_dotprods[0], pn_window_ptr[1] + 1, &dotp[0][1]);
          dotprod_cccf_execute(pn_dotprods[1], pn_window_ptr[0] + 67, &dotp[1][0]);
          dotprod_cccf_execute(pn_dotprods[1], pn_window_ptr[1] + 67, &dotp[1][1]);
          printf("dot products [0][0] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[0][0]), std::imag(dotp[0][0]), std::abs(dotp[0][0]));
          printf("dot products [0][1] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[0][1]), std::imag(dotp[0][1]), std::abs(dotp[0][1]));
          printf("dot products [1][0] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[1][0]), std::imag(dotp[1][0]), std::abs(dotp[1][0]));
          printf("dot products [1][1] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[1][1]), std::imag(dotp[1][1]), std::abs(dotp[1][1]));
          dotprod_cccf_execute(pn_dotprods[0], pn_window_ptr[0] + 2, &dotp[0][0]);
          dotprod_cccf_execute(pn_dotprods[0], pn_window_ptr[1] + 2, &dotp[0][1]);
          dotprod_cccf_execute(pn_dotprods[1], pn_window_ptr[0] + 68, &dotp[1][0]);
          dotprod_cccf_execute(pn_dotprods[1], pn_window_ptr[1] + 68, &dotp[1][1]);
          printf("dot products [0][0] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[0][0]), std::imag(dotp[0][0]), std::abs(dotp[0][0]));
          printf("dot products [0][1] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[0][1]), std::imag(dotp[0][1]), std::abs(dotp[0][1]));
          printf("dot products [1][0] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[1][0]), std::imag(dotp[1][0]), std::abs(dotp[1][0]));
          printf("dot products [1][1] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[1][1]), std::imag(dotp[1][1]), std::abs(dotp[1][1]));
          dotprod_cccf_execute(pn_dotprods[0], pn_window_ptr[0] + 3, &dotp[0][0]);
          dotprod_cccf_execute(pn_dotprods[0], pn_window_ptr[1] + 3, &dotp[0][1]);
          dotprod_cccf_execute(pn_dotprods[1], pn_window_ptr[0] + 65, &dotp[1][0]);
          dotprod_cccf_execute(pn_dotprods[1], pn_window_ptr[1] + 65, &dotp[1][1]);
          printf("dot products [0][0] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[0][0]), std::imag(dotp[0][0]), std::abs(dotp[0][0]));
          printf("dot products [0][1] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[0][1]), std::imag(dotp[0][1]), std::abs(dotp[0][1]));
          printf("dot products [1][0] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[1][0]), std::imag(dotp[1][0]), std::abs(dotp[1][0]));
          printf("dot products [1][1] :%8.4f + i%8.4f, mag :%8.4f\n",
              std::real(dotp[1][1]), std::imag(dotp[1][1]), std::abs(dotp[1][1]));

        std::complex<float> phasing_symb;
        std::complex<float> temp[2][training_seq_len];
        std::complex<float> hell[2];

        // #FIXME
        {
          std::complex<float> dphi_metric[2][2] = {{0.0f, 0.0f}, {0.0f, 0.0f}};
          std::complex<float> r0[2][2] = {{0.0f, 0.0f}, {0.0f, 0.0f}};
          std::complex<float> r1[2][2] = {{0.0f, 0.0f}, {0.0f, 0.0f}};
  
          for(unsigned int i = 0; i < training_seq_len; i++)
          {
            r0[0][0] = r1[0][0];
            r0[0][1] = r1[0][1];
            r0[1][0] = r1[1][0];
            r0[1][1] = r1[1][1];
            r1[0][0] = (*(pn_window_ptr[0] +  0 + i))*preamble_pn[0][i];
            r1[0][1] = (*(pn_window_ptr[1] +  0 + i))*preamble_pn[0][i];
            r1[1][0] = (*(pn_window_ptr[0] + 66 + i))*preamble_pn[1][i];
            r1[1][1] = (*(pn_window_ptr[1] + 66 + i))*preamble_pn[1][i];
            dphi_metric[0][0] += r1[0][0]*std::conj(r0[0][0]);
            dphi_metric[0][1] += r1[0][1]*std::conj(r0[0][1]);
            dphi_metric[1][0] += r1[1][0]*std::conj(r0[1][0]);
            dphi_metric[1][1] += r1[1][1]*std::conj(r0[1][1]);
          }
  
          float dphi_hat[2][2] = {{std::arg(dphi_metric[0][0]),
                                   std::arg(dphi_metric[0][1])},
                                  {std::arg(dphi_metric[1][0]),
                                   std::arg(dphi_metric[1][1])}};
  
          std::complex<float> theta_metric[2][2] = {{0.0f, 0.0f}, {0.0f, 0.0f}};
  
          for(unsigned int i = 0; i < training_seq_len; i++)
          {
            theta_metric[0][0] += (*(pn_window_ptr[0] +  0 + i))*
              std::exp(-dphi_hat[0][0]*i)*preamble_pn[0][i];
            theta_metric[0][1] += (*(pn_window_ptr[1] +  0 + i))*
              std::exp(-dphi_hat[0][1]*i)*preamble_pn[0][i];
            theta_metric[1][0] += (*(pn_window_ptr[0] + 66 + i))*
              std::exp(-dphi_hat[1][0]*i)*preamble_pn[1][i];
            theta_metric[1][1] += (*(pn_window_ptr[1] + 66 + i))*
              std::exp(-dphi_hat[1][1]*i)*preamble_pn[1][i];
          }

          float theta_hat[2][2] = {{std::arg(theta_metric[0][0]),
                                    std::arg(theta_metric[0][1])},
                                   {std::arg(theta_metric[1][0]),
                                    std::arg(theta_metric[1][1])}};
  
          printf("dphi_metric[0][0] :%8.4f + %8.4fi\n",
                 std::real(dphi_metric[0][0]), std::imag(dphi_metric[0][0]));
  
          printf("dphi_metric[0][1] :%8.4f + %8.4fi\n",
                 std::real(dphi_metric[0][1]), std::imag(dphi_metric[0][1]));
          
          printf("dphi_metric[1][0] :%8.4f + %8.4fi\n",
                 std::real(dphi_metric[1][0]), std::imag(dphi_metric[1][0]));
          
          printf("dphi_metric[1][1] :%8.4f + %8.4fi\n",
                 std::real(dphi_metric[1][1]), std::imag(dphi_metric[1][1]));
  
          printf("dphi_hat[0][0] :%8.4f\n", dphi_hat[0][0]);
          printf("dphi_hat[0][1] :%8.4f\n", dphi_hat[0][1]);
          printf("dphi_hat[1][0] :%8.4f\n", dphi_hat[1][0]);
          printf("dphi_hat[1][1] :%8.4f\n", dphi_hat[1][1]);
          
          printf("theta_metric[0][0] :%8.4f + %8.4fi\n",
                 std::real(theta_metric[0][0]), std::imag(theta_metric[0][0]));
          
          printf("theta_metric[0][1] :%8.4f + %8.4fi\n",
                 std::real(theta_metric[0][1]), std::imag(theta_metric[0][1]));
          
          printf("theta_metric[1][0] :%8.4f + %8.4fi\n",
                 std::real(theta_metric[1][0]), std::imag(theta_metric[1][0]));
          
          printf("theta_metric[1][1] :%8.4f + %8.4fi\n",
                 std::real(theta_metric[1][1]), std::imag(theta_metric[1][1]));
          
          printf("theta_hat[0][0] :%8.4f\n", theta_hat[0][0]);
          printf("theta_hat[0][1] :%8.4f\n", theta_hat[0][1]);
          printf("theta_hat[1][0] :%8.4f\n", theta_hat[1][0]);
          printf("theta_hat[1][1] :%8.4f\n", theta_hat[1][1]);

          // initialize fine-tuned nco
          nco_fine_freq[0] = dphi_hat[0][0];
          nco_fine_freq[1] = dphi_hat[0][1];
          nco_fine_phase[0] = theta_hat[0][0];
          nco_fine_phase[1] = theta_hat[0][1];

          nco_crcf_set_frequency(nco_fine[0], nco_fine_freq[0]);
          nco_crcf_set_phase(nco_fine[0], nco_fine_phase[0]);
          nco_crcf_set_frequency(nco_fine[1], nco_fine_freq[1]);
          nco_crcf_set_phase(nco_fine[1], nco_fine_phase[1]);
        }

        for(unsigned int i = 0; i < training_seq_len; i++) {
          nco_crcf_mix_down(nco_fine[0], *(pn_window_ptr[0] + i), temp[0] + i);
          nco_crcf_mix_down(nco_fine[1], *(pn_window_ptr[1] + i), temp[1] + i);

          nco_crcf_pll_step(nco_fine[0], std::imag(temp[0][i]*preamble_pn[0][i]));
          nco_crcf_pll_step(nco_fine[1], std::imag(temp[1][i]*preamble_pn[0][i]));

          nco_crcf_step(nco_fine[0]);
          nco_crcf_step(nco_fine[1]);
        }

        for(unsigned int i = 0; i < training_seq_len; i++) {
          phasing_symb = (i % 2) ? 1.0f : -1.0f;

          nco_crcf_mix_down(nco_fine[0], *(pn_window_ptr[0] + 132 + i), temp[0] + i);
          nco_crcf_mix_down(nco_fine[1], *(pn_window_ptr[1] + 132 + i), temp[1] + i);

          nco_crcf_pll_step(nco_fine[0], std::imag(temp[0][i]*phasing_symb));
          nco_crcf_pll_step(nco_fine[1], std::imag(temp[1][i]*phasing_symb));

          nco_crcf_step(nco_fine[0]);
          nco_crcf_step(nco_fine[1]);
        }

        }*/
