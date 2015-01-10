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

#ifndef MIMO_H
#define MIMO_H

#include <math.h>
#include <assert.h>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <liquid/liquid.h>

#include "liquid_math.h"

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

    class framesync
    {
      private:
        // callback
        unsigned long int frame1_count;
        unsigned long int frame2_count;

        framesync_callback callback;     // user-defined callback function
        void * userdata;                 // user-defined data structure
        framesyncstats_s framestats1;    // frame statistic object
        framesyncstats_s framestats2;    // frame statistic object
        
        // synchronizer objects
        unsigned int seq_len_exp;
        unsigned int seq_len;
        detector_cccf frame1_detector;   // pre-demod detector
        detector_cccf frame2_detector;   // pre-demod detector
        float tau_hat1;                  // fractional timing offset estimate
        float dphi_hat1;                 // carrier frequency offset estimate
        float gamma_hat1;                // channel gain estimate
        float tau_hat2;                  // fractional timing offset estimate
        float dphi_hat2;                 // carrier frequency offset estimate
        float gamma_hat2;                // channel gain estimate
        windowcf buffer1;                // pre-demod buffered samples, size: k*(pn_len+m)
        windowcf buffer2;                // pre-demod buffered samples, size: k*(pn_len+m)
        nco_crcf nco_coarse1;            // coarse carrier frequency recovery
        nco_crcf nco_fine1;              // fine carrier recovery (after demod)
        nco_crcf nco_coarse2;            // coarse carrier frequency recovery
        nco_crcf nco_fine2;              // fine carrier recovery (after demod)
        
        // timing recovery objects, states
        unsigned int k;                  // interp samples/symbol (fixed at 2)
        unsigned int m;                  // interp filter delay (symbols)
        float        beta;               // excess bandwidth factor
        firpfb_crcf mf1;                 // matched filter decimator
        firpfb_crcf dmf1;                // derivative matched filter decimator
        int npfb1;                       // number of filters in symsync
        float pfb_q1;                    //
        float pfb_soft1;                 // soft filterbank index
        int pfb_index1;                  // hard filterbank index
        int pfb_timer1;                  // filterbank output flag
        std::complex<float> symsync_out1;// symbol synchronizer output
        firpfb_crcf mf2;                 // matched filter decimator
        firpfb_crcf dmf2;                // derivative matched filter decimator
        int npfb2;                       // number of filters in symsync
        float pfb_q2;                    //
        float pfb_soft2;                 // soft filterbank index
        int pfb_index2;                  // hard filterbank index
        int pfb_timer2;                  // filterbank output flag
        std::complex<float> symsync_out2;// symbol synchronizer output
        
        // preamble
        std::complex<float> preamble_pn1[64];  // known 64-symbol p/n sequence
        std::complex<float> preamble_rx1[64];  // received p/n symbols
        std::complex<float> preamble_pn2[64];  // known 64-symbol p/n sequence
        std::complex<float> preamble_rx2[64];  // received p/n symbols
        
        // status variables
        enum {
            STATE_DETECTFRAME1=0,           // detect frame (seek p/n sequence)
            STATE_RXPN1,                    // receive p/n sequence
            STATE_DETECTFRAME2,             // detect frame (seek p/n sequence)
            STATE_RXPN2,                    // receive p/n sequence
        } state;
        unsigned int pn1_counter;        // counter: num of p/n syms received
        unsigned int pn2_counter;        // counter: num of p/n syms received

        // push samples through detection stage 1
        void execute_seekpn1(std::complex<float> _x);
        void execute_seekpn2(std::complex<float> _x);

        // update symbol synchronizer internal state (filtered error, index, etc.)
        //  _x      :   input sample
        //  _y      :   output symbol
        int update_symsync1(std::complex<float>   _x,
                           std::complex<float> * _y);
        int update_symsync2(std::complex<float>   _x,
                           std::complex<float> * _y);
        // push buffered p/n sequence through synchronizer
        void pushpn1();
        void pushpn2();

        // push samples through synchronizer, saving received p/n symbols
        void execute_rxpn1(std::complex<float> _x);
        void execute_rxpn2(std::complex<float> _x);

        // once p/n symbols are buffered, estimate residual carrier
        // frequency and phase offsets, push through fine-tuned NCO
        void syncpn1();
        void syncpn2();

      public:
        framesync(framesync_callback, void *, unsigned int,
                  unsigned int, float);
        ~framesync();
        void reset();
        void work(std::complex<float> *, unsigned int);
        unsigned long int get_frame1_count();
        unsigned long int get_frame2_count();
    };
  }
}

#endif
