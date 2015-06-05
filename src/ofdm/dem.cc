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

#include "ofdm.h"

int internal_callback(std::complex<float> * _X,
                      unsigned char *       _p,
                      unsigned int          _M,
                      void *                _userdata)
{
  liquid::ofdm::demodulator * _dem = 
    (liquid::ofdm::demodulator *)_userdata;
  #if DEBUG_OFDMFLEXFRAMESYNC
    printf("******* ofdmflexframesync callback invoked!\n");
  #endif
  _dem->increment_symbol_number();

  #if DEBUG_OFDMFLEXFRAMESYNC
    printf("received symbol %u\n", symbol_number);
  #endif

  // extract symbols
  switch (_dem->get_state()) {
    case RX_STATE_HDR:
      _dem->rxheader(_X);
      break;
    case RX_STATE_PLD:
      _dem->rxpayload(_X);
      break;
    default:
      fprintf(stderr,"error: internal_callback(), unknown/unsupported internal state\n");
      exit(1);
  }
  return 0;
}

namespace liquid {
  namespace ofdm {
    demodulator::demodulator(unsigned int     _M,
                       unsigned int           _cp_len,
                       unsigned int           _taper_len,
                       unsigned char *        _p,
                       framesync_callback     _callback,
                       void *                 _userdata)
    {
      // validate input
      if (_M < 8) {
        fprintf(stderr,"error: modulator::modulator(), number of subcarriers must be at least 8\n");
        exit(1);
      } else if (_M % 2) {
        fprintf(stderr,"error: modulator::modulator(), number of subcarriers must be even 8\n");
        exit(1);
      } else if (_cp_len < 1) {
        fprintf(stderr,"error: modulator::modulator(), cyclic prefix length must be at least 1\n");
        exit(1);
      } else if (_taper_len > _cp_len) {
        fprintf(stderr,"error: modulator::modulator(), taper length cannot exceed cyclic prefix length\n");
        exit(1);
      }

      // set internal properties
      M               = _M;
      cp_len          = _cp_len;
      taper_len       = _taper_len;
      callback        = _callback;
      userdata        = _userdata;

      p = (unsigned char *) malloc (M*sizeof(unsigned char));
      if(_p == NULL)
        ofdmframe_init_default_sctype(M, p);
      else
        memmove(p, _p, M*sizeof(unsigned char));
      ofdmframe_validate_sctype(p, M, &M_null, &M_pilot, &M_data);

      fs = ofdmframesync_create(M, cp_len, taper_len, p,
                                internal_callback,
                                (void *)this);
      p_header = packetizer_create(OFDMFRAME_H_LEN,
                                   OFDMFRAME_H_CRC,
                                   OFDMFRAME_H_EC1,
                                   OFDMFRAME_H_EC2);
      assert(packetizer_get_enc_msg_len(p_header) == OFDMFRAME_H_ENC);
      mod_header = modem_create(OFDMFRAME_H_MOD);
      p_payload = packetizer_create(OFDMFRAME_P_LEN,
                                    OFDMFRAME_P_CRC,
                                    OFDMFRAME_P_EC1,
                                    OFDMFRAME_P_EC2);
      assert(packetizer_get_enc_msg_len(p_payload) == OFDMFRAME_P_ENC);
      mod_payload = modem_create(OFDMFRAME_P_MOD);
    
      div_t d = div(OFDMFRAME_H_SYM, M_data);
      num_header_symbols = d.quot + (d.rem ? 1 : 0);
      d = div(OFDMFRAME_P_SYM, M_data);
      num_payload_symbols = d.quot + (d.rem ? 1 : 0);
      frame_len = 2 + 1 + num_header_symbols + num_payload_symbols;
      payload_dec_len = OFDMFRAME_P_LEN;
      payload_enc_len = OFDMFRAME_P_ENC;
      payload_mod_len = OFDMFRAME_P_SYM;
      header_dec_len = OFDMFRAME_H_LEN;
      header_enc_len = OFDMFRAME_H_ENC;
      header_mod_len = OFDMFRAME_H_SYM;
      reset();
    }
    
    demodulator::~demodulator()
    {
      ofdmframesync_destroy(fs);
      packetizer_destroy(p_header);
      modem_destroy(mod_header);
      packetizer_destroy(p_payload);
      modem_destroy(mod_payload);

      // free internal buffers/arrays
      free(p);
    }

    void demodulator::reset()
    {
      // reset internal state
      state = RX_STATE_HDR;
  
      // reset internal counters
      symbol_number         = 0;
      header_symbol_index   = 0;
      payload_symbol_index  = 0;
      payload_buffer_index  = 0;
      
      // reset error vector magnitude estimate
      evm_hat = 1e-12f;   // slight offset to ensure no log(0)
  
      // reset framestats object
      framesyncstats_init_default(&framestats);
  
      // reset internal OFDM frame synchronizer object
      ofdmframesync_reset(fs);
    }
    
    unsigned int demodulator::get_M()
    {
      return M;
    }
    
    void demodulator::set_M(unsigned int _M)
    {
      M = _M;
      reset();
    }
    
    unsigned int demodulator::get_cp_len()
    {
      return cp_len;
    }
    
    void demodulator::set_cp_len(unsigned int _cp_len)
    {
      cp_len = _cp_len;
      reset();
    }
    
    unsigned int demodulator::get_taper_len()
    {
      return taper_len;
    }
    
    void demodulator::set_taper_len(unsigned int _taper_len)
    {
      taper_len = _taper_len;
      reset();
    }
    
    unsigned char * demodulator::get_p()
    {
      return p;
    }
    
    void demodulator::set_p(unsigned char * _p)
    {
      memmove(p, _p, M*sizeof(unsigned char));
      ofdmframe_validate_sctype(p, M, &M_null, &M_pilot, &M_data);
      reset();
    }
    
    unsigned int demodulator::get_header_dec_len()
    {
      return header_dec_len;
    }
    
    unsigned int demodulator::get_header_enc_len()
    {
      return header_enc_len;
    }
    
    unsigned int demodulator::get_header_mod_len()
    {
      return header_mod_len;
    }
    
    unsigned int demodulator::get_payload_dec_len()
    {
      return payload_dec_len;
    }
    
    unsigned int demodulator::get_payload_enc_len()
    {
      return payload_enc_len;
    }
    
    unsigned int demodulator::get_payload_mod_len()
    {
      return payload_mod_len;
    }
    
    unsigned int demodulator::get_M_null()
    {
      return M_null;
    }
    
    unsigned int demodulator::get_M_pilot()
    {
      return M_pilot;
    }
    
    unsigned int demodulator::get_M_data()
    {
      return M_data;
    }
    
    unsigned int demodulator::get_M_S0()
    {
      return M_S0;
    }
    
    unsigned int demodulator::get_M_S1()
    {
      return M_S1;
    }
    
    unsigned int demodulator::get_num_header_symbols()
    {
      return num_header_symbols;
    }
    
    unsigned int demodulator::get_num_payload_symbols()
    {
      return num_payload_symbols;
    }
    
    unsigned int demodulator::get_frame_len()
    {
      return frame_len;
    }
    
    unsigned int demodulator::get_header_symbol_index()
    {
      return header_symbol_index;
    }
    
    unsigned int demodulator::get_payload_symbol_index()
    {
      return payload_symbol_index;
    }

    void demodulator::demodulate_samples(std::complex<float> * _buff,
                                         unsigned int num_samples)
    {
      ofdmframesync_execute(fs, _buff, num_samples);
    }

    // received signal strength indication
    float demodulator::get_rssi()
    {
      return ofdmframesync_get_rssi(fs);
    }

    // received carrier frequency offset
    float demodulator::get_cfo()
    {
      return ofdmframesync_get_cfo(fs);
    }

    void demodulator::print_config()
    {
      printf("modulator:\n");
      printf("    num subcarriers     :   %-u\n", M);
      printf("      * NULL            :   %-u\n", M_null);
      printf("      * pilot           :   %-u\n", M_pilot);
      printf("      * data            :   %-u\n", M_data);
      printf("    cyclic prefix len   :   %-u\n", cp_len);
      printf("    taper len           :   %-u\n", taper_len);
      printf("    properties:\n");
      if (1) {
        printf("    payload:\n");
        printf("      * decoded bytes   :   %-u\n", payload_dec_len);
        printf("      * encoded bytes   :   %-u\n", payload_enc_len);
        printf("      * modulated syms  :   %-u\n", payload_mod_len);
        printf("    total OFDM symbols  :   %-u\n", frame_len);
        printf("      * S0 symbols      :   %-u @ %u\n", 2, M + cp_len);
        printf("      * S1 symbols      :   %-u @ %u\n", 1, M + cp_len);
        printf("      * header symbols  :   %-u @ %u\n", num_header_symbols, M + cp_len);
        printf("      * payload symbols :   %-u @ %u\n", num_payload_symbols, M + cp_len);
    
        // compute asymptotic spectral efficiency
        unsigned int num_bits = 8*payload_dec_len;
        unsigned int num_samples = (M + cp_len)*(3 + num_header_symbols + num_payload_symbols);
        printf("    spectral efficiency :   %-6.4f b/s/Hz\n", (float)num_bits / (float)num_samples);
      }
    }
    
    void demodulator::print_p()
    {
      unsigned int i;
      std::cout << "Subcarrier Allocation:\n";
      for(i = 0; i < M - 1; i++)
      {
        if(p[i] == OFDMFRAME_SCTYPE_NULL)
          std::cout << "N, ";
        else if(p[i] == OFDMFRAME_SCTYPE_PILOT)
          std::cout << "P, ";
        else if(p[i] == OFDMFRAME_SCTYPE_DATA)
          std::cout << "D, ";
      }
      if(p[i] == OFDMFRAME_SCTYPE_NULL)
        std::cout << "N\n";
      else if(p[i] == OFDMFRAME_SCTYPE_PILOT)
        std::cout << "P\n";
      else if(p[i] == OFDMFRAME_SCTYPE_DATA)
        std::cout << "D\n";
    }

    // receive header data
    void demodulator::rxheader(std::complex<float> * _X)
    {
      #if DEBUG_OFDMFLEXFRAMESYNC
        printf("  ofdmflexframesync extracting header...\n");
      #endif

      // demodulate header symbols
      unsigned int i;
      int sctype;
      for (i = 0; i < M; i++) {
        // subcarrier type (PILOT/NULL/DATA)
        sctype = p[i];

        // ignore pilot and null subcarriers
        if (sctype == OFDMFRAME_SCTYPE_DATA) {
          // unload header symbols
          // demodulate header symbol
          unsigned int sym;
          modem_demodulate(mod_header, _X[i], &sym);
          header_m[header_symbol_index] = sym;
          header_symbol_index++;
          //printf("  extracting symbol %3u / %3u (x = %8.5f + j%8.5f)\n", _q->header_symbol_index, OFDMFLEXFRAME_H_SYM, crealf(_X[i]), cimagf(_X[i]));

          // get demodulator error vector magnitude
          float evm = modem_get_demodulator_evm(mod_header);
          evm_hat += evm*evm;

          // header extracted
          if (header_symbol_index == OFDMFRAME_H_SYM) {
            // decode header
            decode_header();
          
            // compute error vector magnitude estimate
            framestats.evm = 10*log10f(evm_hat/OFDMFRAME_H_SYM);

            // invoke callback if header is invalid
            if(header_valid)
              state = RX_STATE_PLD;
            else {
              //printf("**** header invalid!\n");
              // set framestats internals
              framestats.rssi             = ofdmframesync_get_rssi(fs);
              framestats.cfo              = ofdmframesync_get_cfo(fs);
              framestats.framesyms        = NULL;
              framestats.num_framesyms    = 0;
              framestats.mod_scheme       = LIQUID_MODEM_UNKNOWN;
              framestats.mod_bps          = 0;
              framestats.check            = LIQUID_CRC_UNKNOWN;
              framestats.fec0             = LIQUID_FEC_UNKNOWN;
              framestats.fec1             = LIQUID_FEC_UNKNOWN;

              // invoke callback method
              callback(header_d,
                       header_valid,
                       NULL,
                       0,
                       0,
                       framestats,
                       userdata);

              reset();
            }
            break;
          }
        }
      }
    }

    // decode header
    void demodulator::decode_header()
    {
      // pack 1-bit header symbols into 8-bit bytes
      unsigned int bps = modulation_types[OFDMFRAME_H_MOD].bps;
      unsigned int num_written;
      liquid_repack_bytes(header_m,
                          bps,
                          OFDMFRAME_H_SYM,
                          header_e,
                          8,
                          OFDMFRAME_H_ENC,
                          &num_written);
      assert(num_written == OFDMFRAME_H_ENC);

      // unscramble header
      unscramble_data(header_e, OFDMFRAME_H_ENC);

      // run packet decoder
      header_valid = packetizer_decode(p_header,
                                       header_e,
                                       header_d);

      #if DEBUG_OFDMFLEXFRAMESYNC
        printf("****** header extracted [%s]\n", header_valid ? "valid" : "INVALID!");
      #endif
      if (!header_valid)
        return;

      unsigned int n = OFDMFRAME_H_USR;

      // first byte is for expansion/version validation
      if (header_d[n+0] != OFDMFRAME_VERSN) {
        fprintf(stderr,"warning: decode_header(), invalid framing version\n");
        header_valid = 0;
      }

      // strip off payload length
      unsigned int payload_len = 
        (header_d[n+1] << 8) | (header_d[n+2]);

      // strip off modulation scheme/depth
      unsigned int mod_scheme = header_d[n+3];
      if (mod_scheme != OFDMFRAME_P_MOD) {
        fprintf(stderr,"warning: decode_header(), invalid modulation scheme\n");
        header_valid = 0;
      }
      if (payload_len != OFDMFRAME_P_LEN) {
        fprintf(stderr,"warning: decode_header(), invalid payload length\n");
        header_valid = 0;
      }

      // strip off CRC, forward error-correction schemes
      //  CRC     : most-significant 3 bits of [n+4]
      //  fec0    : least-significant 5 bits of [n+4]
      //  fec1    : least-significant 5 bits of [n+5]
      unsigned int check = (header_d[n+4] >> 5 ) & 0x07;
      unsigned int fec0  = (header_d[n+4]      ) & 0x1f;
      unsigned int fec1  = (header_d[n+5]      ) & 0x1f;

      // validate properties
      if (check != OFDMFRAME_P_CRC) {
        fprintf(stderr,"warning: decode_header(), decoded CRC exceeds available\n");
        check = LIQUID_CRC_UNKNOWN;
        header_valid = 0;
      }
      if (fec0 != OFDMFRAME_P_EC1) {
        fprintf(stderr,"warning: decode_header(), decoded FEC (inner) exceeds available\n");
        fec0 = LIQUID_FEC_UNKNOWN;
        header_valid = 0;
      }
      if (fec1 != OFDMFRAME_P_EC2) {
        fprintf(stderr,"warning: decode_header(), decoded FEC (outer) exceeds available\n");
        fec1 = LIQUID_FEC_UNKNOWN;
        header_valid = 0;
      }

      // print results
      #if DEBUG_OFDMFLEXFRAMESYNC
        printf("    properties:\n");
        printf("      * mod scheme      :   %s\n", modulation_types[mod_scheme].fullname);
        printf("      * fec (inner)     :   %s\n", fec_scheme_str[fec0][1]);
        printf("      * fec (outer)     :   %s\n", fec_scheme_str[fec1][1]);
        printf("      * CRC scheme      :   %s\n", crc_scheme_str[check][1]);
        printf("      * payload length  :   %u bytes\n", payload_len);
      #endif
    }

// receive payload data
    void demodulator::rxpayload(std::complex<float> * _X)
    {
      // demodulate paylod symbols
      unsigned int i;
      int sctype;
      for (i = 0; i < M; i++) {
        // subcarrier type (PILOT/NULL/DATA)
        sctype = p[i];

        // ignore pilot and null subcarriers
        if (sctype == OFDMFRAME_SCTYPE_DATA) {
          // unload payload symbols
          unsigned int sym;
          modem_demodulate(mod_payload, _X[i], &sym);

          // pack decoded symbol into array
          liquid_pack_array(payload_e,
                            payload_enc_len,
                            payload_buffer_index,
                            OFDMFRAME_P_BPS,
                            sym);

          // increment...
          payload_buffer_index += OFDMFRAME_P_BPS;

          // increment symbol counter
          payload_symbol_index++;

          if (payload_symbol_index == payload_mod_len) {
            // payload extracted

            // decode payload
            payload_valid = packetizer_decode(p_payload,
                                              payload_e,
                                              payload_d);
            #if DEBUG_OFDMFLEXFRAMESYNC
              printf("****** payload extracted [%s]\n", payload_valid ? "valid" : "INVALID!");
            #endif

            // ignore callback if set to NULL
            if (callback == NULL) {
              reset();
              break;
            }

            // set framestats internals
            framestats.rssi             = ofdmframesync_get_rssi(fs);
            framestats.cfo              = ofdmframesync_get_cfo(fs);
            framestats.framesyms        = NULL;
            framestats.num_framesyms    = 0;
            framestats.mod_scheme       = OFDMFRAME_P_MOD;
            framestats.mod_bps          = OFDMFRAME_P_BPS;
            framestats.check            = OFDMFRAME_P_CRC;
            framestats.fec0             = OFDMFRAME_P_EC1;
            framestats.fec1             = OFDMFRAME_P_EC2;

            // invoke callback method
            callback(header_d,
                     header_valid,
                     payload_d,
                     payload_dec_len,
                     payload_valid,
                     framestats,
                     userdata);


            // reset object
            reset();
            break;
          }
        }
      }
    }

    void demodulator::increment_symbol_number()
    {
      symbol_number++;
    }

    rx_states demodulator::get_state()
    {
      return state;
    }
  }
}
