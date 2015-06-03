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

#define OFDMFRAME_VERSN   (104)
#define OFDMFRAME_H_USR   (8)
#define OFDMFRAME_H_LEN   (14)
#define OFDMFRAME_H_CRC   (LIQUID_CRC_32)
#define OFDMFRAME_H_EC1   (LIQUID_FEC_GOLAY2412)
#define OFDMFRAME_H_EC2   (LIQUID_FEC_NONE)
#define OFDMFRAME_H_ENC   (36)
#define OFDMFRAME_H_MOD   (LIQUID_MODEM_BPSK)
#define OFDMFRAME_H_BPS   (1)
#define OFDMFRAME_H_SYM   (288)
#define OFDMFRAME_P_LEN   (1024)
#define OFDMFRAME_P_CRC   (LIQUID_CRC_32)
#define OFDMFRAME_P_EC1   (LIQUID_FEC_NONE)
#define OFDMFRAME_P_EC2   (LIQUID_FEC_NONE)
#define OFDMFRAME_P_ENC   (1028)                // FIXME
#define OFDMFRAME_P_MOD   (LIQUID_MODEM_QPSK)
#define OFDMFRAME_P_BPS   (2)
#define OFDMFRAME_P_SYM   (4112)

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
        unsigned int num_header_symbols;  // number of header OFDM symbols
        unsigned int num_payload_symbols; // number of payload OFDM symbols
        unsigned int frame_len;           // length of the frame(# OFDM Symbols)
    
        modem mod_header;                         // header modulator
        packetizer p_header;                      // header packetizer
        unsigned int header_dec_len;              // header length
        unsigned int header_enc_len;              // header length encoded
        unsigned int header_mod_len;              // header mod length
        unsigned char header_d[OFDMFRAME_H_LEN];  // header data
        unsigned char header_e[OFDMFRAME_H_ENC];  // header encoded
        unsigned char header_m[OFDMFRAME_H_SYM];  // header symbols
    
        modem mod_payload;                        // payload modulator
        packetizer p_payload;                     // payload packetizer
        unsigned int payload_dec_len;             // payload length
        unsigned int payload_enc_len;             // payload length encoded
        unsigned int payload_mod_len;             // payload mod length
        unsigned char payload_e[OFDMFRAME_P_ENC]; // encoded bytes
        unsigned char payload_m[OFDMFRAME_P_SYM]; // encoded bytes
    
        ofdmframegen fg;
    
        unsigned int symbol_number;               // symbol counter
        enum {
          STATE_S0A = 0,    // write S0 symbol (first)
          STATE_S0B,        // write S0 symbol (second)
          STATE_S1,         // write S1 symbol
          STATE_HDR,        // write header symbols
          STATE_PLD         // write payload symbols
        } state;
    
        int frame_assembled;  // frame assembled flag
        int frame_complete;   // frame completed flag
        unsigned int header_symbol_index;
        unsigned int payload_symbol_index;
    
      public:
        // constructor - destructor
        modulator(unsigned int       _M,
                  unsigned int       _cp_len,
                  unsigned int       _taper_len,
                  unsigned char *    _p);
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
    
        std::complex<float> * X;          // frequency-domain buffer
        unsigned int num_header_symbols;  // number of header OFDM symbols
        unsigned int num_payload_symbols; // number of payload OFDM symbols
        unsigned int frame_len;           // length of the frame(# OFDM Symbols)
    
        modem mod_header;                         // header modulator
        packetizer p_header;                      // header packetizer
        unsigned int header_dec_len;              // header length
        unsigned int header_enc_len;              // header length encoded
        unsigned int header_mod_len;              // header mod length
        unsigned char header_d[OFDMFRAME_H_LEN];  // header data
        unsigned char header_e[OFDMFRAME_H_ENC];  // header encoded
        unsigned char header_m[OFDMFRAME_H_SYM];  // header symbols
    
        modem mod_payload;                        // payload modulator
        packetizer p_payload;                     // payload packetizer
        unsigned int payload_dec_len;             // payload length
        unsigned int payload_enc_len;             // payload length encoded
        unsigned int payload_mod_len;             // payload mod length
        unsigned char payload_e[OFDMFRAME_P_ENC]; // encoded bytes
        unsigned char payload_m[OFDMFRAME_P_SYM]; // encoded bytes

        ofdmframesync fs;             // frame synchronizer object
        framesync_callback callback;
        void * userdata;

      public:
        demodulator(unsigned int       _M,
                    unsigned int       _cp_len,
                    unsigned int       _taper_len,
                    unsigned char *    _p,
                    framesync_callback _callback,
                    void *             _userdata);
    
        // destructor
        ~demodulator();
    
    };
  }
}

#endif
