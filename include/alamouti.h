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

#ifndef ALAMOUTI_H
#define ALAMOUTI_H

#include <math.h>
#include <assert.h>
#include <complex>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <liquid/liquid.h>
#include <limits.h>

#include "liquid_math.h"

namespace liquid {
  namespace alamouti {

    class framegen
    {
      private:
        unsigned int training_seq_len;
        unsigned int training_symbol_counter;
        unsigned int payload_len;
        unsigned int payload_symbol_counter;
        unsigned int frame_len;

        std::complex<float> * sps1;
        std::complex<float> * sps2;
        std::complex<float> * pn1;
        std::complex<float> * pn2;
        std::complex<float> * payload;

        firinterp_crcf interp1;
        firinterp_crcf interp2;
        unsigned int k;
        unsigned int m;
        float beta;
        float gain1;
        float gain2;

        enum {
            STATE_TXPN1 = 0,
            STATE_TXPN2,
            STATE_TXPAYLOAD,
            STATE_PNZ,
        } state;

      public:
        framegen(unsigned int, unsigned int, float);
        ~framegen();
        void reset();
        void set_gains(float, float);
        unsigned int get_training_seq_len();
        unsigned int get_frame_len();
        unsigned int work(std::complex<float> **);
    };

    class framesync
    {
      private:
        //
        unsigned long int training_seq_n[2][2];
        unsigned int training_match_index[2][2];
        // synchronizer objects
        unsigned int training_seq_len;
        unsigned int payload_len;
        unsigned int frame_len;
        bool frame_detect_flag[2][2];
        detector_cccf frame_detector[2][2];   // pre-demod detector
        float tau_hat[2][2];                  // fractional timing offset estimate
        float dphi_hat[2][2];                 // carrier frequency offset estimate
        float gamma_hat[2][2];                // channel gain estimate
        windowcf buffer[2];                   // pre-demod buffered samples, size: k*(pn_len+m)
        float nco_coarse_freq;
        nco_crcf nco_coarse[2];            // coarse carrier frequency recovery
        nco_crcf nco_fine[2];              // fine carrier recovery (after demod)
        
        // timing recovery objects, states
        unsigned int k;                  // interp samples/symbol (fixed at 2)
        unsigned int m;                  // interp filter delay (symbols)
        float        beta;               // excess bandwidth factor

        unsigned int npfb;
        firpfb_crcf mf[2];                 // matched filter decimator
        firpfb_crcf dmf[2];                // derivative matched filter decimator
        float pfb_q[2];                    //
        float pfb_soft[2];                 // soft filterbank index
        int pfb_index[2];                  // hard filterbank index
        int pfb_timer[2];                  // filterbank output flag
        std::complex<float> symsync_out[2];// symbol synchronizer output
        
        // preamble
        std::complex<float> * preamble_pn[2];  // known 63-symbol p/n sequence
        std::complex<float> * preamble_rx[2][2];  // received p/n symbols
        
        // status variables
        enum {
            STATE_DETECTFRAME1=0,           // detect frame (seek p/n sequence)
            STATE_RXPN1,                    // receive p/n sequence
            STATE_DETECTFRAME2,             // detect frame (seek p/n sequence)
            STATE_RXPN2,                    // receive p/n sequence
            STATE_DETECTFRAME3,             // detect frame (seek p/n sequence)
            STATE_RXPN3,                    // receive p/n sequence
            STATE_RXPAYLOAD
        } state;

        void execute_seekpn1(std::complex<float> ** _x, unsigned int i);
        void execute_seekpn2(std::complex<float> ** _x, unsigned int i);

        // update symbol synchronizer internal state (filtered error, index, etc.)
        //  _x      :   input sample
        //  _y      :   output symbol
        int update_symsync1(std::complex<float>   _x,
                           std::complex<float> * _y);
        int update_symsync2(std::complex<float>   _x,
                           std::complex<float> * _y);
        int update_symsync3(std::complex<float>   _x,
                           std::complex<float> * _y);
        // push buffered p/n sequence through synchronizer
        void pushpn1();
        void pushpn2();

        // push samples through synchronizer, saving received p/n symbols
        void execute_rxpn1(std::complex<float> ** _x, unsigned int i);
        void execute_rxpn2(std::complex<float> ** _x, unsigned int i);

        // once p/n symbols are buffered, estimate residual carrier
        // frequency and phase offsets, push through fine-tuned NCO
        void syncpn1();
        void syncpn2();

        void execute_rxpayload(std::complex<float> ** _x, unsigned int i);
        void update_training_sequence_n();
        int find_index_of_corr();
        void training_match_index_diffs();

      public:
        framesync(unsigned int, unsigned int, float);
        ~framesync();
        void reset();
        void work(std::complex<float> **, unsigned int);
    };
  }
}

#endif
