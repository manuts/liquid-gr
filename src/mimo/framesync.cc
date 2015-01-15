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
    framesync::framesync(framesync_callback _callback,
                         void * _userdata,
                         unsigned int _k,
                         unsigned int _m,
                         float _beta)
    {
      callback = _callback;
      userdata = _userdata;
      seq_len_exp = 6;
      seq_len = (unsigned int)(pow(2, seq_len_exp) - 1);
      msequence ms1 = msequence_create(seq_len_exp, 0x005b, 1);
      msequence ms2 = msequence_create(seq_len_exp, 0x0043, 1);
      preamble_pn1 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*seq_len); 
      preamble_pn2 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*seq_len); 
      preamble_pn3 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*seq_len); 
      preamble_rx1 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*seq_len); 
      preamble_rx2 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*seq_len); 
      preamble_rx3 = (std::complex<float> *)malloc(sizeof(std::complex<float>)*seq_len); 
      for (unsigned int i = 0; i < seq_len; i++){
        preamble_pn1[i] = (msequence_advance(ms1)) ? 1.0f : -1.0f;
        preamble_pn2[i] = (msequence_advance(ms2)) ? 1.0f : -1.0f;
        preamble_pn3[i] = preamble_pn1[i] + liquid::math::I*preamble_pn2[i];
      }
      msequence_destroy(ms1);
      msequence_destroy(ms2);

      // interpolate p/n sequence with matched filter
      k     = _k;        // samples/symbol
      m     = _m;        // filter delay (symbols)
      beta  = _beta;    // excess bandwidth factor

      std::complex<float> seq1[k*seq_len];
      std::complex<float> seq2[k*seq_len];
      std::complex<float> seq3[k*seq_len];
      firinterp_crcf interp1 = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);
      firinterp_crcf interp2 = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);
      firinterp_crcf interp3 = firinterp_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER,k,m,beta,0);
      for (unsigned int i = 0; i < seq_len + m; i++) {
        // compensate for filter delay
        if (i < m) {
            firinterp_crcf_execute(interp1, preamble_pn1[i], &seq1[0]);
            firinterp_crcf_execute(interp2, preamble_pn2[i], &seq2[0]);
            firinterp_crcf_execute(interp3, preamble_pn3[i], &seq3[0]);
        }
        else {
            firinterp_crcf_execute(interp1, preamble_pn1[i%seq_len], &seq1[k*(i - m)]);
            firinterp_crcf_execute(interp2, preamble_pn2[i%seq_len], &seq2[k*(i - m)]);
            firinterp_crcf_execute(interp3, preamble_pn3[i%seq_len], &seq3[k*(i - m)]);
        }
      }
      firinterp_crcf_destroy(interp1);
      firinterp_crcf_destroy(interp2);
      firinterp_crcf_destroy(interp3);

      // create frame detector
      float threshold = 0.4f;     // detection threshold
      float dphi_max  = 0.05f;    // maximum carrier offset allowable
      frame1_detector = detector_cccf_create(seq1, k*seq_len, threshold, dphi_max);
      frame2_detector = detector_cccf_create(seq2, k*seq_len, threshold, dphi_max);
      frame3_detector = detector_cccf_create(seq3, k*seq_len, threshold, dphi_max);
      buffer1 = windowcf_create(k*(seq_len + m));
      buffer2 = windowcf_create(k*(seq_len + m));
      buffer3 = windowcf_create(k*(seq_len + m));

      // create symbol timing recovery filters
      npfb1 = 32;   // number of filters in the bank
      mf1   = firpfb_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER, npfb1,k,m,beta);
      dmf1  = firpfb_crcf_create_drnyquist(LIQUID_FIRFILT_ARKAISER,npfb1,k,m,beta);
      npfb2 = 32;   // number of filters in the bank
      mf2   = firpfb_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER, npfb2,k,m,beta);
      dmf2  = firpfb_crcf_create_drnyquist(LIQUID_FIRFILT_ARKAISER,npfb2,k,m,beta);
      npfb3 = 32;   // number of filters in the bank
      mf3   = firpfb_crcf_create_rnyquist(LIQUID_FIRFILT_ARKAISER, npfb3,k,m,beta);
      dmf3  = firpfb_crcf_create_drnyquist(LIQUID_FIRFILT_ARKAISER,npfb3,k,m,beta);

      // create down-coverters for carrier phase tracking
      nco_coarse1 = nco_crcf_create(LIQUID_NCO);
      nco_fine1   = nco_crcf_create(LIQUID_VCO);
      nco_crcf_pll_set_bandwidth(nco_fine1, 0.05f);
      nco_coarse2 = nco_crcf_create(LIQUID_NCO);
      nco_fine2   = nco_crcf_create(LIQUID_VCO);
      nco_crcf_pll_set_bandwidth(nco_fine2, 0.05f);
      nco_coarse3 = nco_crcf_create(LIQUID_NCO);
      nco_fine3   = nco_crcf_create(LIQUID_VCO);
      nco_crcf_pll_set_bandwidth(nco_fine3, 0.05f);

      frame1_count = 0;
      frame2_count = 0;
      frame3_count = 0;

      // reset state
      state = STATE_DETECTFRAME1;
      reset1();
      reset2();
      reset3();

      if(__SAVE__){
        f_pn1 = fopen("/tmp/rx_pn1", "wb");
        f_pn2 = fopen("/tmp/rx_pn2", "wb");
        f_pn3 = fopen("/tmp/rx_pn3", "wb");
      }
    }

    void framesync::reset1()
    {
      // reset binary pre-demod synchronizer
      detector_cccf_reset(frame1_detector);

      // clear pre-demod buffer
      windowcf_clear(buffer1);

      // reset carrier recovery objects
      nco_crcf_reset(nco_coarse1);
      nco_crcf_reset(nco_fine1);

      // reset symbol timing recovery state
      firpfb_crcf_reset(mf1);
      firpfb_crcf_reset(dmf1);
      pfb_q1 = 0.0f;   // filtered error signal

      pn1_counter      = 0;

      // reset frame statistics
      framestats1.evm = 0.0f;
    }

    void framesync::reset2()
    {
      // reset binary pre-demod synchronizer
      detector_cccf_reset(frame2_detector);

      // clear pre-demod buffer
      windowcf_clear(buffer2);

      // reset carrier recovery objects
      nco_crcf_reset(nco_coarse2);
      nco_crcf_reset(nco_fine2);

      // reset symbol timing recovery state
      firpfb_crcf_reset(mf2);
      firpfb_crcf_reset(dmf2);
      pfb_q2 = 0.0f;   // filtered error signal

      pn2_counter      = 0;

      // reset frame statistics
      framestats2.evm = 0.0f;
    }

    void framesync::reset3()
    {
      // reset binary pre-demod synchronizer
      detector_cccf_reset(frame3_detector);

      // clear pre-demod buffer
      windowcf_clear(buffer3);

      // reset carrier recovery objects
      nco_crcf_reset(nco_coarse3);
      nco_crcf_reset(nco_fine3);

      // reset symbol timing recovery state
      firpfb_crcf_reset(mf3);
      firpfb_crcf_reset(dmf3);
      pfb_q3 = 0.0f;   // filtered error signal

      pn3_counter      = 0;

      // reset frame statistics
      framestats3.evm = 0.0f;
    }

    framesync::~framesync()
    {
      if(__SAVE__){
        fclose(f_pn1);
        fclose(f_pn2);
        fclose(f_pn3);
      }
      // destroy synchronization objects
      free(preamble_pn1);
      free(preamble_pn2);
      free(preamble_pn3);
      free(preamble_rx1);
      free(preamble_rx2);
      free(preamble_rx3);
      detector_cccf_destroy(frame1_detector);  // frame detector
      detector_cccf_destroy(frame2_detector);  // frame detector
      detector_cccf_destroy(frame3_detector);  // frame detector
      windowcf_destroy(buffer1);               // p/n sample buffer
      windowcf_destroy(buffer2);               // p/n sample buffer
      windowcf_destroy(buffer3);               // p/n sample buffer
      firpfb_crcf_destroy(mf1);                // matched filter
      firpfb_crcf_destroy(dmf1);               // derivative matched filter
      nco_crcf_destroy(nco_coarse1);           // coarse NCO
      nco_crcf_destroy(nco_fine1);             // fine-tuned NCO
      firpfb_crcf_destroy(mf2);                // matched filter
      firpfb_crcf_destroy(dmf2);               // derivative matched filter
      nco_crcf_destroy(nco_coarse2);           // coarse NCO
      nco_crcf_destroy(nco_fine2);             // fine-tuned NCO
      firpfb_crcf_destroy(mf3);                // matched filter
      firpfb_crcf_destroy(dmf3);               // derivative matched filter
      nco_crcf_destroy(nco_coarse3);           // coarse NCO
      nco_crcf_destroy(nco_fine3);             // fine-tuned NCO
    }

    void framesync::work(std::complex<float> * _x,
                         unsigned int    _n)
    {
      for (unsigned int i = 0; i < _n; i++) {
        switch (state) {
          case STATE_DETECTFRAME1:
            // detect frame (look for p/n sequence)
            execute_seekpn1(_x[i]);
            break;
          case STATE_RXPN1:
            // detect frame (look for p/n sequence)
            execute_rxpn1(_x[i]);
            break;
          case STATE_DETECTFRAME2:
            // detect frame (look for p/n sequence)
            execute_seekpn2(_x[i]);
            break;
          case STATE_RXPN2:
            // detect frame (look for p/n sequence)
            execute_rxpn2(_x[i]);
            break;
          case STATE_DETECTFRAME3:
            // detect frame (look for p/n sequence)
            execute_seekpn3(_x[i]);
            break;
          case STATE_RXPN3:
            // detect frame (look for p/n sequence)
            execute_rxpn3(_x[i]);
            break;
        default:
            fprintf(stderr,"error: framesync64_exeucte(), unknown/unsupported state\n");
            exit(1);
        }
      }
    }

    void framesync::execute_seekpn1(std::complex<float> _x)
    {
      // push sample into pre-demod p/n sequence buffer
      windowcf_push(buffer1, _x);

      // push through pre-demod synchronizer
      int detected = detector_cccf_correlate(frame1_detector,
                                             _x,
                                             &tau_hat1,
                                             &dphi_hat1,
                                             &gamma_hat1);

      // check if frame has been detected
      if (detected) {
        printf("***** frame1 detected! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB\n",
               tau_hat1, dphi_hat1, 20*log10f(gamma_hat1));

        // push buffered samples through synchronizer
        // NOTE: this will set internal state appropriately
        //       to STATE_SEEKPN
        pushpn1();
      }
    }

    // push buffered p/n sequence through synchronizer
    void framesync::pushpn1()
    {
      // reset filterbanks
      firpfb_crcf_reset(mf1);
      firpfb_crcf_reset(dmf1);

      // read buffer
      std::complex<float> * rc;
      windowcf_read(buffer1, &rc);

      // compute delay and filterbank index
      //  tau_hat < 0 :   delay = 2*k*m-1, index = round(   tau_hat *npfb), flag = 0
      //  tau_hat > 0 :   delay = 2*k*m-2, index = round((1-tau_hat)*npfb), flag = 0
      assert(tau_hat1 < 0.5f && tau_hat1 > -0.5f);
      unsigned int delay = 2*k*m - 1; // samples to buffer before computing output
      pfb_soft1       = -tau_hat1*npfb1;
      pfb_index1      = (int) roundf(pfb_soft1);
      while (pfb_index1 < 0) {
        delay         -= 1;
        pfb_index1 += npfb1;
        pfb_soft1  += npfb1;
      }
      pfb_timer1 = 0;

      // set coarse carrier frequency offset
      nco_crcf_set_frequency(nco_coarse1, dphi_hat1);
    
      unsigned int buffer_len = (seq_len + m)*k;
      for (unsigned int i = 0; i < buffer_len; i++) {
        if (i < delay) {
          std::complex<float> y;
          nco_crcf_mix_down(nco_coarse1, rc[i]*0.5f/gamma_hat1, &y);
          nco_crcf_step(nco_coarse1);

          // push initial samples into filterbanks
          firpfb_crcf_push(mf1,  y);
          firpfb_crcf_push(dmf1, y);
        }
        else {
            // run remaining samples through p/n sequence recovery
          if(state == STATE_DETECTFRAME1)
          {
            state = STATE_RXPN1;
            execute_rxpn1(rc[i]);
          }
          else if(state == STATE_RXPN1)
          {
            execute_rxpn1(rc[i]);
          }
          else if(state == STATE_DETECTFRAME2)
          {
            return;
          }
        }
      }
    }

    void framesync::execute_rxpn1(std::complex<float> _x)
    {
      // validate input
      if(pn1_counter == seq_len)
      {
        fprintf(stderr, "Buff1 Full\n");
        return;
      }

      // mix signal down
      std::complex<float> y;
      nco_crcf_mix_down(nco_coarse1, _x*0.5f/gamma_hat1, &y);
      nco_crcf_step(nco_coarse1);

      // update symbol synchronizer
      std::complex<float> mf_out = 0.0f;
      int sample_available = update_symsync1(y, &mf_out);

      // compute output if timeout
      if (sample_available) {
        // save output in p/n symbols buffer
        preamble_rx1[pn1_counter] = mf_out;

        // update p/n counter
        pn1_counter++;

        if (pn1_counter == seq_len) {
          syncpn1();
          state = STATE_DETECTFRAME2;
          frame1_count++;
          if(__SAVE__)
            fwrite((void *)preamble_rx1, sizeof(std::complex<float>), seq_len, f_pn1);
          reset1();
        }
      }
    }

    int framesync::update_symsync1(std::complex<float>   _x,
                           std::complex<float> * _y)
    {
      // push sample into filterbanks
      firpfb_crcf_push(mf1,  _x);
      firpfb_crcf_push(dmf1, _x);

      //
      std::complex<float> mf_out  = 0.0f;    // matched-filter output
      std::complex<float> dmf_out = 0.0f;    // derivatived matched-filter output

      int sample_available = 0;

      // compute output if timeout
      if (pfb_timer1 <= 0) {
        sample_available = 1;

        // reset timer
        pfb_timer1 = k;  // k samples/symbol

        firpfb_crcf_execute(mf1, pfb_index1, &mf_out);
        firpfb_crcf_execute(dmf1, pfb_index1, &dmf_out);

        // update filtered timing error
        // hi  bandwidth parameters: {0.92, 1.20}, about 100 symbols settling time
        // med bandwidth parameters: {0.98, 0.20}, about 200 symbols settling time
        // lo  bandwidth parameters: {0.99, 0.05}, about 500 symbols settling time
        pfb_q1 = 0.99f*pfb_q1 + 0.05f*std::real(std::conj(mf_out)*dmf_out );

        // accumulate error into soft filterbank value
        pfb_soft1 += pfb_q1;

        // compute actual filterbank index
        pfb_index1 = roundf(pfb_soft1);

        // contrain index to be in [0, npfb-1]
        while (pfb_index1 < 0) {
          pfb_index1 += npfb1;
          pfb_soft1  += npfb1;

          // adjust pfb output timer
          pfb_timer1--;
        }
        while (pfb_index1 > npfb1 - 1) {
          pfb_index1 -= npfb1;
          pfb_soft1  -= npfb1;

          // adjust pfb output timer
          pfb_timer1++;
        }
      }

      // decrement symbol timer
      pfb_timer1--;

      // set output and return
      *_y = mf_out;
    
      return sample_available;
    }

    void framesync::syncpn1()
    {
      // estimate residual carrier frequency offset from p/n symbols
      std::complex<float> dphi_metric = 0.0f;
      std::complex<float> r0 = 0.0f;
      std::complex<float> r1 = 0.0f;
      for (unsigned int i = 0; i < seq_len; i++) {
        r0 = r1;
        r1 = preamble_rx1[i]*preamble_pn1[i];
        dphi_metric += r1 * std::conj(r0);
      }
      float dphi_hat = std::arg(dphi_metric);

      // estimate carrier phase offset from p/n symbols
      std::complex<float> theta_metric = 0.0f;
      for (unsigned int i = 0; i < seq_len; i++)
        theta_metric += preamble_rx1[i]*liquid::math::cexpjf(-dphi_hat*i)*preamble_pn1[i];
      float theta_hat = std::arg(theta_metric);
      // TODO: compute gain correction factor

      // initialize fine-tuned nco
      nco_crcf_set_frequency(nco_fine1, dphi_hat);
      nco_crcf_set_phase(nco_fine1, theta_hat);

      // correct for carrier offset, pushing through phase-locked loop
      for (unsigned int i = 0; i < seq_len; i++) {
        // mix signal down
        nco_crcf_mix_down(nco_fine1, preamble_rx1[i], &preamble_rx1[i]);
        
        // push through phase-locked loop
        float phase_error = std::imag(preamble_rx1[i]*preamble_pn1[i]);
        nco_crcf_pll_step(nco_fine1, phase_error);
        // update nco phase
        nco_crcf_step(nco_fine1);
      }
    }

    void framesync::execute_seekpn2(std::complex<float> _x)
    {
      // push sample into pre-demod p/n sequence buffer
      windowcf_push(buffer2, _x);

      // push through pre-demod synchronizer
      int detected = detector_cccf_correlate(frame2_detector,
                                             _x,
                                             &tau_hat2,
                                             &dphi_hat2,
                                             &gamma_hat2);

      // check if frame has been detected
      if (detected) {
        printf("***** frame2 detected! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB\n",
               tau_hat2, dphi_hat2, 20*log10f(gamma_hat2));
        printf("***** Estimate diffs!! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB\n",
               tau_hat1 - tau_hat2, dphi_hat1 - dphi_hat2, 20*log10f(gamma_hat1/gamma_hat2));

        // push buffered samples through synchronizer
        // NOTE: this will set internal state appropriately
        //       to STATE_SEEKPN
        pushpn2();
      }
    }

    // push buffered p/n sequence through synchronizer
    void framesync::pushpn2()
    {
      // reset filterbanks
      firpfb_crcf_reset(mf2);
      firpfb_crcf_reset(dmf2);

      // read buffer
      std::complex<float> * rc;
      windowcf_read(buffer2, &rc);

      // compute delay and filterbank index
      //  tau_hat < 0 :   delay = 2*k*m-1, index = round(   tau_hat *npfb), flag = 0
      //  tau_hat > 0 :   delay = 2*k*m-2, index = round((1-tau_hat)*npfb), flag = 0
      assert(tau_hat2 < 0.5f && tau_hat2 > -0.5f);
      unsigned int delay = 2*k*m - 1; // samples to buffer before computing output
      pfb_soft2       = -tau_hat2*npfb2;
      pfb_index2      = (int) roundf(pfb_soft2);
      while (pfb_index2 < 0) {
        delay         -= 1;
        pfb_index2 += npfb2;
        pfb_soft2  += npfb2;
      }
      pfb_timer2 = 0;

      // set coarse carrier frequency offset
      nco_crcf_set_frequency(nco_coarse2, dphi_hat2);
    
      unsigned int buffer_len = (seq_len + m)*k;
      for (unsigned int i = 0; i < buffer_len; i++) {
        if (i < delay) {
          std::complex<float> y;
          nco_crcf_mix_down(nco_coarse2, rc[i]*0.5f/gamma_hat2, &y);
          nco_crcf_step(nco_coarse2);

          // push initial samples into filterbanks
          firpfb_crcf_push(mf2,  y);
          firpfb_crcf_push(dmf2, y);
        }
        else {
            // run remaining samples through p/n sequence recovery
          if(state == STATE_DETECTFRAME2)
          {
            state = STATE_RXPN2;
            execute_rxpn1(rc[i]);
          }
          else if(state == STATE_RXPN2)
          {
            execute_rxpn2(rc[i]);
          }
          else if(state == STATE_DETECTFRAME3)
          {
            return;
          }
        }
      }
    }

    void framesync::execute_rxpn2(std::complex<float> _x)
    {
      if(pn2_counter == seq_len)
      {
        fprintf(stderr, "Buff2 Full\n");
        return;
      }

      // mix signal down
      std::complex<float> y;
      nco_crcf_mix_down(nco_coarse2, _x*0.5f/gamma_hat2, &y);
      nco_crcf_step(nco_coarse2);

      // update symbol synchronizer
      std::complex<float> mf_out = 0.0f;
      int sample_available = update_symsync2(y, &mf_out);

      // compute output if timeout
      if (sample_available) {
        // save output in p/n symbols buffer
        preamble_rx2[pn2_counter] = mf_out;

        // update p/n counter
        pn2_counter++;

        if (pn2_counter == seq_len) {
          syncpn2();
          state = STATE_DETECTFRAME3;
          frame2_count++;
          if(__SAVE__)
            fwrite((void *)preamble_rx2, sizeof(std::complex<float>), seq_len, f_pn2);
          reset2();
        }
      }
    }

    int framesync::update_symsync2(std::complex<float>   _x,
                           std::complex<float> * _y)
    {
      // push sample into filterbanks
      firpfb_crcf_push(mf2,  _x);
      firpfb_crcf_push(dmf2, _x);

      //
      std::complex<float> mf_out  = 0.0f;    // matched-filter output
      std::complex<float> dmf_out = 0.0f;    // derivatived matched-filter output

      int sample_available = 0;

      // compute output if timeout
      if (pfb_timer2 <= 0) {
        sample_available = 1;

        // reset timer
        pfb_timer2 = k;  // k samples/symbol

        firpfb_crcf_execute(mf2, pfb_index2, &mf_out);
        firpfb_crcf_execute(dmf2, pfb_index2, &dmf_out);

        // update filtered timing error
        // hi  bandwidth parameters: {0.92, 1.20}, about 100 symbols settling time
        // med bandwidth parameters: {0.98, 0.20}, about 200 symbols settling time
        // lo  bandwidth parameters: {0.99, 0.05}, about 500 symbols settling time
        pfb_q2 = 0.99f*pfb_q2 + 0.05f*std::real( std::conj(mf_out)*dmf_out );

        // accumulate error into soft filterbank value
        pfb_soft2 += pfb_q2;

        // compute actual filterbank index
        pfb_index2 = roundf(pfb_soft2);

        // contrain index to be in [0, npfb-1]
        while (pfb_index2 < 0) {
          pfb_index2 += npfb2;
          pfb_soft2  += npfb2;

          // adjust pfb output timer
          pfb_timer2--;
        }
        while (pfb_index2 > npfb2 - 1) {
          pfb_index2 -= npfb2;
          pfb_soft2  -= npfb2;

          // adjust pfb output timer
          pfb_timer2++;
        }
      }

      // decrement symbol timer
      pfb_timer2--;

      // set output and return
      *_y = mf_out;
    
      return sample_available;
    }

    void framesync::syncpn2()
    {
      // estimate residual carrier frequency offset from p/n symbols
      std::complex<float> dphi_metric = 0.0f;
      std::complex<float> r0 = 0.0f;
      std::complex<float> r1 = 0.0f;
      for (unsigned int i = 0; i < seq_len; i++) {
        r0 = r1;
        r1 = preamble_rx2[i]*preamble_pn2[i];
        dphi_metric += r1 * std::conj(r0);
      }
      float dphi_hat = std::arg(dphi_metric);

      // estimate carrier phase offset from p/n symbols
      std::complex<float> theta_metric = 0.0f;
      for (unsigned int i = 0; i < seq_len; i++)
        theta_metric += preamble_rx2[i]*liquid::math::cexpjf(-dphi_hat*i)*preamble_pn2[i];
      float theta_hat = std::arg(theta_metric);
      // TODO: compute gain correction factor

      // initialize fine-tuned nco
      nco_crcf_set_frequency(nco_fine2, dphi_hat);
      nco_crcf_set_phase(nco_fine2, theta_hat);

      // correct for carrier offset, pushing through phase-locked loop
      for (unsigned int i = 0; i < seq_len; i++) {
        // mix signal down
        nco_crcf_mix_down(nco_fine2, preamble_rx2[i], &preamble_rx2[i]);
        
        // push through phase-locked loop
        float phase_error = std::imag(preamble_rx2[i]*preamble_pn2[i]);
        nco_crcf_pll_step(nco_fine2, phase_error);
        // update nco phase
        nco_crcf_step(nco_fine2);
      }
    }

    void framesync::execute_seekpn3(std::complex<float> _x)
    {
      // push sample into pre-demod p/n sequence buffer
      windowcf_push(buffer3, _x);

      // push through pre-demod synchronizer
      int detected = detector_cccf_correlate(frame3_detector,
                                             _x,
                                             &tau_hat3,
                                             &dphi_hat3,
                                             &gamma_hat3);

      // check if frame has been detected
      if (detected) {
        printf("***** frame3 detected! tau-hat:%8.4f, dphi-hat:%8.4f, gamma:%8.2f dB\n",
               tau_hat3, dphi_hat3, 20*log10f(gamma_hat3));

        // push buffered samples through synchronizer
        // NOTE: this will set internal state appropriately
        //       to STATE_SEEKPN
        pushpn3();
      }
    }

    // push buffered p/n sequence through synchronizer
    void framesync::pushpn3()
    {
      // reset filterbanks
      firpfb_crcf_reset(mf3);
      firpfb_crcf_reset(dmf3);

      // read buffer
      std::complex<float> * rc;
      windowcf_read(buffer3, &rc);

      // compute delay and filterbank index
      //  tau_hat < 0 :   delay = 2*k*m-1, index = round(   tau_hat *npfb), flag = 0
      //  tau_hat > 0 :   delay = 2*k*m-2, index = round((1-tau_hat)*npfb), flag = 0
      assert(tau_hat3 < 0.5f && tau_hat3 > -0.5f);
      unsigned int delay = 2*k*m - 1; // samples to buffer before computing output
      pfb_soft3       = -tau_hat3*npfb3;
      pfb_index3      = (int) roundf(pfb_soft3);
      while (pfb_index3 < 0) {
        delay         -= 1;
        pfb_index3 += npfb3;
        pfb_soft3  += npfb3;
      }
      pfb_timer3 = 0;

      // set coarse carrier frequency offset
      nco_crcf_set_frequency(nco_coarse3, dphi_hat3);
    
      unsigned int buffer_len = (seq_len + m)*k;
      for (unsigned int i = 0; i < buffer_len; i++) {
        if (i < delay) {
          std::complex<float> y;
          nco_crcf_mix_down(nco_coarse3, rc[i]*0.5f/gamma_hat3, &y);
          nco_crcf_step(nco_coarse3);

          // push initial samples into filterbanks
          firpfb_crcf_push(mf3,  y);
          firpfb_crcf_push(dmf3, y);
        }
        else {
            // run remaining samples through p/n sequence recovery
          if(state == STATE_DETECTFRAME3)
          {
            state = STATE_RXPN3;
            execute_rxpn3(rc[i]);
          }
          else if(state == STATE_RXPN3)
          {
            execute_rxpn3(rc[i]);
          }
          else if(state == STATE_DETECTFRAME1)
          {
            return;
          }
        }
      }
    }

    void framesync::execute_rxpn3(std::complex<float> _x)
    {
      // validate input
      if(pn3_counter == seq_len)
      {
        fprintf(stderr, "Buff3 Full\n");
        return;
      }

      // mix signal down
      std::complex<float> y;
      nco_crcf_mix_down(nco_coarse3, _x*0.5f/gamma_hat3, &y);
      nco_crcf_step(nco_coarse3);

      // update symbol synchronizer
      std::complex<float> mf_out = 0.0f;
      int sample_available = update_symsync3(y, &mf_out);

      // compute output if timeout
      if (sample_available) {
        // save output in p/n symbols buffer
        preamble_rx3[pn3_counter] = mf_out;

        // update p/n counter
        pn3_counter++;

        if (pn3_counter == seq_len) {
          syncpn3();
          state = STATE_DETECTFRAME1;
          frame3_count++;
          if(__SAVE__)
            fwrite((void *)preamble_rx3, sizeof(std::complex<float>), seq_len, f_pn3);
          reset3();
        }
      }
    }

    int framesync::update_symsync3(std::complex<float>   _x,
                           std::complex<float> * _y)
    {
      // push sample into filterbanks
      firpfb_crcf_push(mf3,  _x);
      firpfb_crcf_push(dmf3, _x);

      //
      std::complex<float> mf_out  = 0.0f;    // matched-filter output
      std::complex<float> dmf_out = 0.0f;    // derivatived matched-filter output

      int sample_available = 0;

      // compute output if timeout
      if (pfb_timer3 <= 0) {
        sample_available = 1;

        // reset timer
        pfb_timer3 = k;  // k samples/symbol

        firpfb_crcf_execute(mf3, pfb_index3, &mf_out);
        firpfb_crcf_execute(dmf3, pfb_index3, &dmf_out);

        // update filtered timing error
        // hi  bandwidth parameters: {0.92, 1.20}, about 100 symbols settling time
        // med bandwidth parameters: {0.98, 0.20}, about 200 symbols settling time
        // lo  bandwidth parameters: {0.99, 0.05}, about 500 symbols settling time
        pfb_q3 = 0.99f*pfb_q3 + 0.05f*std::real(std::conj(mf_out)*dmf_out );

        // accumulate error into soft filterbank value
        pfb_soft3 += pfb_q3;

        // compute actual filterbank index
        pfb_index3 = roundf(pfb_soft3);

        // contrain index to be in [0, npfb-1]
        while (pfb_index3 < 0) {
          pfb_index3 += npfb3;
          pfb_soft3  += npfb3;

          // adjust pfb output timer
          pfb_timer3--;
        }
        while (pfb_index3 > npfb3 - 1) {
          pfb_index3 -= npfb3;
          pfb_soft3  -= npfb3;

          // adjust pfb output timer
          pfb_timer3++;
        }
      }

      // decrement symbol timer
      pfb_timer3--;

      // set output and return
      *_y = mf_out;
    
      return sample_available;
    }

    void framesync::syncpn3()
    {
      // estimate residual carrier frequency offset from p/n symbols
      std::complex<float> dphi_metric = 0.0f;
      std::complex<float> r0 = 0.0f;
      std::complex<float> r1 = 0.0f;
      for (unsigned int i = 0; i < seq_len; i++) {
        r0 = r1;
        r1 = preamble_rx3[i]*preamble_pn3[i];
        dphi_metric += r1 * std::conj(r0);
      }
      float dphi_hat = std::arg(dphi_metric);

      // estimate carrier phase offset from p/n symbols
      std::complex<float> theta_metric = 0.0f;
      for (unsigned int i = 0; i < seq_len; i++)
        theta_metric += preamble_rx3[i]*liquid::math::cexpjf(-dphi_hat*i)*preamble_pn3[i];
      float theta_hat = std::arg(theta_metric);
      // TODO: compute gain correction factor

      // initialize fine-tuned nco
      nco_crcf_set_frequency(nco_fine3, dphi_hat);
      nco_crcf_set_phase(nco_fine3, theta_hat);

      // correct for carrier offset, pushing through phase-locked loop
      for (unsigned int i = 0; i < seq_len; i++) {
        // mix signal down
        nco_crcf_mix_down(nco_fine3, preamble_rx3[i], &preamble_rx3[i]);
        
        // push through phase-locked loop
        float phase_error = std::imag(preamble_rx3[i]*preamble_pn3[i]);
        nco_crcf_pll_step(nco_fine3, phase_error);
        // update nco phase
        nco_crcf_step(nco_fine3);
      }
    }

    unsigned long int framesync::get_frame1_count() {
      return frame1_count;
    }

    unsigned long int framesync::get_frame2_count() {
      return frame2_count;
    }

    unsigned long int framesync::get_frame3_count() {
      return frame3_count;
    }
  }   // namespace mimo
}     // namespace liquid
