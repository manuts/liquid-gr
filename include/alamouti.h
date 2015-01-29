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

#define NUM_FILTS 32
#define RANGE 30

namespace liquid {
  namespace alamouti {

    class tx_buff_transformer
    {
      private:
        bool initialized;
        unsigned int usrp_buff_len;
        unsigned int packer_buff_len;
        unsigned int items_txed;
        unsigned int hist;
        std::complex<float> * usrp_buff[2];
        std::complex<float> * packer_buff[2];
        std::complex<float> * history[2];
      public:

        tx_buff_transformer();
        ~tx_buff_transformer();
        void set_usrp_buff_len(unsigned int);
        void set_packer_buff_len(unsigned int);
        void set_usrp_buff(std::complex<float> **);
        void set_packer_buff(std::complex<float> **);
        void copy_buff();
        bool fill_buff();
    };

    class rx_buff_transformer
    {
      private:
        bool initialized;
        unsigned int usrp_buff_len;
        unsigned int unpacker_buff_len;
        unsigned int items_rxed;
        unsigned int hist;
        std::complex<float> * usrp_buff[2];
        std::complex<float> * unpacker_buff[2];
        std::complex<float> * history[2];
      public:

        rx_buff_transformer();
        ~rx_buff_transformer();
        void set_usrp_buff_len(unsigned int);
        void set_unpacker_buff_len(unsigned int);
        void set_usrp_buff(std::complex<float> **);
        void set_unpacker_buff(std::complex<float> **);
        void copy_buff();
        bool dump_buff();
    };

    class delay
    {
      private:
        std::complex<float> * hist;
        unsigned int d;
      public:
        delay(unsigned int);
        ~delay();
        void set_delay(unsigned int);
        void work(std::complex<float> *, unsigned int);
    };

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

    class channel_estimator
    {
      private:
        std::complex<float> * training_seq;
        unsigned int training_seq_len;
        float threshold;
        float dphi_max;
        unsigned int m;
        float dphi_step;
        float dphi[NUM_FILTS];
        float rxy[NUM_FILTS];
        dotprod_cccf * dp;

        unsigned int loc_max_index, glo_max_index, glo_max_freq_index;
        float glo_max, glo_max_freq;

      public:
        channel_estimator(std::complex<float> * _training_seq,
                          unsigned int _training_seq_len,
                          float _threshold,
                          float _dphi_max);
        ~channel_estimator();
        int find_corr_index(std::complex<float> *, unsigned int);
        void reset();
    };

    class framesync
    {
      private:
        //
        channel_estimator * estimator;
        unsigned long int training_seq_n[2][2];
        unsigned long int num_frames_detected;
        unsigned long int errors[2];
        int training_match_index[2][2];
        unsigned int pn_count;
        unsigned int payload_count;
        unsigned int work_exec_cycle;
        // synchronizer objects
        unsigned int wait_for_pn2;
        unsigned int training_seq_len;
        unsigned int payload_len;
        unsigned int frame_len;
        bool frame_detect_flag[2][2];
        detector_cccf frame_detector[2][2];   // pre-demod detector
        float tau_hat[2][2];                  // fractional timing offset estimate
        float dphi_hat[2][2];                 // carrier frequency offset estimate
        float gamma_hat[2][2];                // channel gain estimate
        windowcf pn_window[2];                // pre-demod buffered samples, size: k*(pn_len+m)
        windowcf pn_window_unscaled[2];       // pre-demod buffered samples, size: k*(pn_len+m)
        windowcf payload_window[2];           // pre-demod buffered samples, size: k*(pn_len+m)
        dotprod_cccf pn_dotprods[2];
        std::complex<float> * rx_sig[2];
        std::complex<float> * last_rx_sig[2];
        std::complex<float> * curr_rx_sig[2];
        int rx_sig_index;
        float nco_coarse_freq;
        float nco_fine_freq[2];
        float nco_fine_phase[2];
        nco_crcf nco_coarse[2];            // coarse carrier frequency recovery
        nco_crcf nco_fine[2];              // fine carrier recovery (after demod)
        float pow[2];

        FILE * out_files[2];

        modem demod[2];
        unsigned char * rx_payload[2];
        unsigned char * expected_payload;
        
        // timing recovery objects, states
        unsigned int k;                  // interp samples/symbol (fixed at 2)
        unsigned int m;                  // interp filter delay (symbols)
        float        beta;               // excess bandwidth factor

        unsigned int npfb;
        firpfb_crcf mf[2];                 // matched filter decimator
        firpfb_crcf dmf[2];                // derivative matched filter decimator
        float pfb_q;                       //
        float pfb_soft;                    // soft filterbank index
        int pfb_index;                     // hard filterbank index
        int pfb_timer;                     // filterbank output flag
        std::complex<float> symsync_out[2];// symbol synchronizer output
        
        // preamble
        std::complex<float> * preamble_pn[2];  // known 63-symbol p/n sequence
        std::complex<float> * preamble_rx[2][2];  // received p/n symbols
        
        // status variables
        enum {
            STATE_DETECTFRAME1=0,           // detect frame (seek p/n sequence)
            STATE_DETECTFRAME2,             // detect frame (seek p/n sequence)
            STATE_RXPN,                    // receive p/n sequence
            STATE_RXPAYLOAD
        } state;

        void execute_seekpn1(std::complex<float> ** _x, unsigned int i);
        void execute_seekpn2(std::complex<float> ** _x, unsigned int i);

        // update symbol synchronizer internal state (filtered error, index, etc.)
        //  _x      :   input sample
        //  _y      :   output symbol
        int update_symsync(std::complex<float> *,
                           std::complex<float> *);
        // push buffered p/n sequence through synchronizer
        void pushpn();

        void execute_rxpayload(std::complex<float> ** _x, unsigned int i);
        void execute_rxpn(std::complex<float> ** _x, unsigned int i);
        void update_training_sequence_n();
        int find_index_of_corr();
        void training_match_index_diffs1();
        void training_match_index_diffs2();
        void find_channel_diffs();
      
        void sync_pn();
        void decode_payload();

      public:
        framesync(unsigned int, unsigned int, float);
        ~framesync();
        void reset();
        void work(std::complex<float> **, unsigned int);
        unsigned int get_frame_len();
        void print_frame_struct();
        unsigned long int get_num_frames_detected();
        unsigned long int get_num_errors1();
        unsigned long int get_num_errors2();
        unsigned long int get_num_bits_detected();
    };
  }
}

#endif
