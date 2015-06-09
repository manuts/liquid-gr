/*
 * Copyright (c) 2014, 2015 Manu T S
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

#ifndef OFDM_H
#define OFDM_H

#include <math.h>
#include <complex>
#include <liquid/liquid.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fftw3.h>
#include <assert.h>
#include <volk/volk.h>

//#define OFDMFRAME_VERSN   (104)
//#define OFDMFRAME_H_USR   (8)
//#define OFDMFRAME_H_LEN   (14)
//#define OFDMFRAME_H_CRC   (LIQUID_CRC_32)
//#define OFDMFRAME_H_EC1   (LIQUID_FEC_GOLAY2412)
//#define OFDMFRAME_H_EC2   (LIQUID_FEC_NONE)
//#define OFDMFRAME_H_ENC   (36)
//#define OFDMFRAME_H_MOD   (LIQUID_MODEM_BPSK)
//#define OFDMFRAME_H_BPS   (1)
//#define OFDMFRAME_H_SYM   (288)
//#define OFDMFRAME_P_LEN   (512)
//#define OFDMFRAME_P_CRC   (LIQUID_CRC_NONE)
//#define OFDMFRAME_P_EC1   (LIQUID_FEC_NONE)
//#define OFDMFRAME_P_EC2   (LIQUID_FEC_NONE)
//#define OFDMFRAME_P_ENC   (512)                // FIXME
//#define OFDMFRAME_P_MOD   (LIQUID_MODEM_QPSK)
//#define OFDMFRAME_P_BPS   (2)
//#define OFDMFRAME_P_SYM   (2048)

typedef enum {RX_STATE_HDR = 0,  // receive header
              RX_STATE_PLD       // receive payload
             } rx_states;
typedef enum {TX_STATE_S0A = 0,    // write S0 symbol (first)
              TX_STATE_S0B,        // write S0 symbol (second)
              TX_STATE_S1,         // write S1 symbol
              TX_STATE_HDR,        // write header symbols
              TX_STATE_PLD         // write payload symbols
             } tx_states;

typedef struct {
  unsigned int      VERSN;        // frame structure version
  unsigned int      H_USR;        // user modifiable header bytes
  unsigned int      H_LEN;        // total length of header
  fec_scheme        H_EC1;        // first error correction for header
  fec_scheme        H_EC2;        // second error correction for header
  crc_scheme        H_CRC;        // header crc
  modulation_scheme H_MOD;        // header modulation scheme
  unsigned int      P_LEN;        // length of payload
  fec_scheme        P_EC1;        // first error correction for paylaod
  fec_scheme        P_EC2;        // second error correction for payload
  crc_scheme        P_CRC;        // payload crc
  modulation_scheme P_MOD;        // payload modulation scheme
} OFDMFRAME_STRUCT;

namespace liquid {
  namespace ofdm {
    class modulator
    {
      private:
        unsigned int M;                   // number of subcarriers
        unsigned int cp_len;              // cyclic prefix length
        unsigned int taper_len;           // tapering window length
        unsigned char * p;                // subcarrier allocation
        
        unsigned int M_null;              // number of null subcarriers
        unsigned int M_pilot;             // number of pilot subcarriers
        unsigned int M_data;              // number of data subcarriers
        unsigned int M_S0;                // number of S0 subcarriers
        unsigned int M_S1;                // number of S1 subcarriers
    
        std::complex<float> * X;          // frequency-domain buffer
        std::complex<float> tx_gain;      // DSP gain

        OFDMFRAME_STRUCT frame_struct;
    
        modem mod_header;                         // header modulator
        packetizer p_header;                      // header packetizer
        unsigned int header_enc_len;              // header length encoded
        unsigned int header_mod_len;              // header mod length
        unsigned int header_bps;                  // header bits per symbol
        unsigned char * header_d;  // header data
        unsigned char * header_e;  // header encoded
        unsigned char * header_m;  // header symbols
    
        modem mod_payload;                        // payload modulator
        packetizer p_payload;                     // payload packetizer
        unsigned int payload_enc_len;             // payload length encoded
        unsigned int payload_mod_len;             // payload mod length
        unsigned int payload_bps;                 // payload bits per symbol
        unsigned char * payload_e; // encoded bytes
        unsigned char * payload_m; // encoded bytes
    

        unsigned int num_header_symbols;  // number of header OFDM symbols
        unsigned int num_payload_symbols; // number of payload OFDM symbols
        unsigned int frame_len;           // length of the frame(# OFDM Symbols)
        unsigned int symbol_number;       // symbol counter
        unsigned int header_symbol_index;
        unsigned int payload_symbol_index;
        int frame_assembled;  // frame assembled flag
        int frame_complete;   // frame completed flag
        tx_states state;

        ofdmframegen fg;
    
      public:
        // constructor - destructor
        modulator(unsigned int       _M,
                  unsigned int       _cp_len,
                  unsigned int       _taper_len,
                  unsigned char *    _p,
                  OFDMFRAME_STRUCT * _frame_struct);
        ~modulator();
    
        // get-set methods
        unsigned int get_M();
        void set_M(unsigned int _M);
        unsigned int get_cp_len();
        void set_cp_len(unsigned int _cp_len);
        unsigned int get_taper_len();
        void set_taper_len(unsigned int _taper_len);
        unsigned char * get_p();
        void set_p(unsigned char * _p);
        void set_tx_gain(std::complex<float> _tx_gain);
        std::complex<float> get_tx_gain();
        void print_p();
        unsigned int get_M_null();
        unsigned int get_M_pilot();
        unsigned int get_M_data();
        unsigned int get_M_S0();
        unsigned int get_M_S1();
        unsigned int get_num_header_symbols();
        unsigned int get_num_payload_symbols();
        unsigned int get_payload_dec_len();
        unsigned int get_payload_enc_len();
        unsigned int get_payload_mod_len();
        unsigned int get_h_usr_len();
        unsigned int get_header_dec_len();
        unsigned int get_header_enc_len();
        unsigned int get_header_mod_len();
        unsigned int get_frame_len();
        int get_frame_assembled();
        int get_frame_complete();
        unsigned int get_header_symbol_index();
        unsigned int get_payload_symbol_index();
    
        // work methods
        void assemble_frame(unsigned char * _header,
                            unsigned char * _payload);
        int write_symbol(std::complex<float> * _buffer);
        void encode_header();
        void modulate_header();
        // helper methods
        void reset();
        void print_config();
    
        void write_S0a(std::complex<float> * _buffer);
        void write_S0b(std::complex<float> * _buffer);
        void write_S1(std::complex<float> * _buffer);
        void write_header(std::complex<float> * _buffer);
        void write_payload(std::complex<float> * _buffer);
        void assemble_output_samples(std::complex<float> * _buffer);
    };

    class demodulator
    {
      private:
        unsigned int M;                   // number of subcarriers
        unsigned int cp_len;              // cyclic prefix length
        unsigned int taper_len;           // tapering window length
        unsigned char * p;                // subcarrier allocation
        
        unsigned int M_null;              // number of null subcarriers
        unsigned int M_pilot;             // number of pilot subcarriers
        unsigned int M_data;              // number of data subcarriers
        unsigned int M_S0;                // number of S0 subcarriers
        unsigned int M_S1;                // number of S1 subcarriers
    
        OFDMFRAME_STRUCT frame_struct;
    
        modem mod_header;                         // header modulator
        packetizer p_header;                      // header packetizer
        unsigned int header_enc_len;              // header length encoded
        unsigned int header_mod_len;              // header mod length
        unsigned int header_bps;                  // header bits per symbol
        unsigned char * header_d;  // header data
        unsigned char * header_e;  // header encoded
        unsigned char * header_m;  // header symbols
    
        modem mod_payload;                        // payload modulator
        packetizer p_payload;                     // payload packetizer
        unsigned int payload_enc_len;             // payload length encoded
        unsigned int payload_mod_len;             // payload mod length
        unsigned int payload_bps;                 // payload bits per symbol
        unsigned char * payload_d; // encoded bytes
        unsigned char * payload_e; // encoded bytes
        unsigned char * payload_m; // encoded bytes
    

        unsigned int num_header_symbols;  // number of header OFDM symbols
        unsigned int num_payload_symbols; // number of payload OFDM symbols
        unsigned int frame_len;           // length of the frame(# OFDM Symbols)
        unsigned int symbol_number;       // symbol counter
        unsigned int header_symbol_index;
        unsigned int payload_symbol_index;
        unsigned int payload_buffer_index;
        rx_states state;
        int header_valid;
        int payload_valid;

        ofdmframesync fs;                   // frame synchronizer object
        framesync_callback callback;    // callback
        framesyncstats_s framestats;        // frame statistic object
        void * userdata;
        float evm_hat;                      // average error vector magnitude

      public:
        demodulator(unsigned int            _M,
                    unsigned int            _cp_len,
                    unsigned int            _taper_len,
                    unsigned char *         _p,
                    framesync_callback      _callback,
                    void *                  _userdata,
                    OFDMFRAME_STRUCT *      _frame_struct);
    
        // destructor
        ~demodulator();

        // get-set methods
        unsigned int get_M();
        void set_M(unsigned int _M);
        unsigned int get_cp_len();
        void set_cp_len(unsigned int _cp_len);
        unsigned int get_taper_len();
        void set_taper_len(unsigned int _taper_len);
        unsigned char * get_p();
        void set_p(unsigned char * _p);
        void print_p();
        unsigned int get_M_null();
        unsigned int get_M_pilot();
        unsigned int get_M_data();
        unsigned int get_M_S0();
        unsigned int get_M_S1();
        unsigned int get_num_header_symbols();
        unsigned int get_num_payload_symbols();
        unsigned int get_payload_enc_len();
        unsigned int get_payload_mod_len();
        unsigned int get_header_dec_len();
        unsigned int get_header_enc_len();
        unsigned int get_header_mod_len();
        unsigned int get_frame_len();
        unsigned int get_header_symbol_index();
        unsigned int get_payload_symbol_index();
        void increment_symbol_number();
        rx_states get_state();

        void demodulate_samples(std::complex<float> * _buff,
                                unsigned int num_samples);
        void print_config();
        float get_rssi();
        float get_cfo();
        void reset();
        friend int internal_callback(std::complex<float> * _X,
                                     unsigned char *       _p,
                                     unsigned int          _M,
                                     void *                _userdata);
        void rxheader(std::complex<float> * _X);
        void rxpayload(std::complex<float> * _X);
        void decode_header();
    };
  }
}

#endif
