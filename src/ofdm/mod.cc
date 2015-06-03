/*
 * Copyright (c) 2014, 2015 Manu T S
 *
 * This is free software: you can redistribute it and/or modulatorify
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

namespace liquid {
  namespace ofdm {
    modulator::modulator(unsigned int       _M,
                         unsigned int       _cp_len,
                         unsigned int       _taper_len,
                         unsigned char *    _p)
    {
      // validate input
      if (_M < 8) {
        fprintf(stderr,"error: modulator::modulator(), number of subcarriers must be at least 8\n");
        throw 0;
      } else if (_cp_len < 1) {
        fprintf(stderr,"error: modulator::modulator(), cyclic prefix length must be at least 1\n");
        throw 0;
      } else if (_taper_len > _cp_len) {
        fprintf(stderr,"error: modulator::modulator(), taper length cannot exceed cyclic prefix length\n");
        throw 0;
      }
    
      // set internal properties
      M               = _M;
      cp_len          = _cp_len;
      taper_len       = _taper_len;
    
      // allocate memory
      X = (std::complex<float> *) malloc (M*sizeof(std::complex<float>));
      p = (unsigned char *) malloc (M*sizeof(unsigned char));
      if(_p == NULL)
        ofdmframe_init_default_sctype(M, p);
      else
        memmove(p, _p, M*sizeof(unsigned char));
      ofdmframe_validate_sctype(p, M, &M_null, &M_pilot, &M_data);
    
      fg = ofdmframegen_create(M, cp_len, taper_len, p);
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
    
    modulator::~modulator()
    {
      free(X);
      free(p);
      ofdmframegen_destroy(fg);
      packetizer_destroy(p_header);
      packetizer_destroy(p_payload);
      modem_destroy(mod_header);
      modem_destroy(mod_payload);
    }
    
    void modulator::reset()
    {
      symbol_number = 0;
      state = STATE_S0A;
      frame_assembled = 0;
      frame_complete = 0;
      header_symbol_index = 0;
      payload_symbol_index = 0;
      ofdmframegen_reset(fg);
    }
    
    unsigned int modulator::get_M()
    {
      return M;
    }
    
    void modulator::set_M(unsigned int _M)
    {
      M = _M;
      reset();
    }
    
    unsigned int modulator::get_cp_len()
    {
      return cp_len;
    }
    
    void modulator::set_cp_len(unsigned int _cp_len)
    {
      cp_len = _cp_len;
      reset();
    }
    
    unsigned int modulator::get_taper_len()
    {
      return taper_len;
    }
    
    void modulator::set_taper_len(unsigned int _taper_len)
    {
      taper_len = _taper_len;
      reset();
    }
    
    unsigned char * modulator::get_p()
    {
      return p;
    }
    
    void modulator::set_p(unsigned char * _p)
    {
      memmove(p, _p, M*sizeof(unsigned char));
      ofdmframe_validate_sctype(p, M, &M_null, &M_pilot, &M_data);
      reset();
    }
    
    unsigned int modulator::get_header_dec_len()
    {
      return header_dec_len;
    }
    
    unsigned int modulator::get_header_enc_len()
    {
      return header_enc_len;
    }
    
    unsigned int modulator::get_header_mod_len()
    {
      return header_mod_len;
    }
    
    unsigned int modulator::get_payload_dec_len()
    {
      return payload_dec_len;
    }
    
    unsigned int modulator::get_payload_enc_len()
    {
      return payload_enc_len;
    }
    
    unsigned int modulator::get_payload_mod_len()
    {
      return payload_mod_len;
    }
    
    unsigned int modulator::get_M_null()
    {
      return M_null;
    }
    
    unsigned int modulator::get_M_pilot()
    {
      return M_pilot;
    }
    
    unsigned int modulator::get_M_data()
    {
      return M_data;
    }
    
    unsigned int modulator::get_M_S0()
    {
      return M_S0;
    }
    
    unsigned int modulator::get_M_S1()
    {
      return M_S1;
    }
    
    unsigned int modulator::get_num_header_symbols()
    {
      return num_header_symbols;
    }
    
    unsigned int modulator::get_num_payload_symbols()
    {
      return num_payload_symbols;
    }
    
    unsigned int modulator::get_frame_len()
    {
      return frame_len;
    }
    
    int modulator::get_frame_assembled()
    {
      return frame_assembled;
    }
    
    int modulator::get_frame_complete()
    {
      return frame_complete;
    }
    
    unsigned int modulator::get_header_symbol_index()
    {
      return header_symbol_index;
    }
    
    unsigned int modulator::get_payload_symbol_index()
    {
      return payload_symbol_index;
    }

    void modulator::assemble_output_samples(
            std::complex<float> * _buffer)
    {
      unsigned int frame_index = 0;
      assert(frame_assembled);
      while(!write_symbol(_buffer + (M + cp_len)*frame_index))
        frame_index++;
    }

    int modulator::write_symbol(std::complex<float> * _buffer)
    {
      // check if frame is actually assembled
      if ( !frame_assembled ) {
        fprintf(stderr,"warning: write_symbol(), frame not assembled\n");
        return 1;
      }
        // increment symbol counter
      symbol_number++;
  
      switch (state) {
      case STATE_S0A:
        // write S0 symbol (first)
        write_S0a(_buffer);
        break;
    
        case STATE_S0B:
          // write S0 symbol (second)
          write_S0b(_buffer);
          break;
    
        case STATE_S1:
          // write S1 symbols
          write_S1(_buffer);
          break;
    
        case STATE_HDR:
          // write header symbols
          write_header(_buffer);
          break;
    
        case STATE_PLD:
          // write payload symbols
          write_payload(_buffer);
          break;
    
        default:
          fprintf(stderr,"error: write_symbol(), unknown/unsupported internal state\n");
          exit(1);
        }
        if (frame_complete) {
          reset();
          return 1;
        }
        return 0;
    }
    
    void modulator::assemble_frame(unsigned char * _header,
                                   unsigned char * _payload)
    {
      memmove(header_d, _header, OFDMFRAME_H_USR*sizeof(unsigned char));
      encode_header();
      modulate_header();
      packetizer_encode(p_payload, _payload, payload_e);
      memset(payload_m, 0x00, payload_mod_len);
      unsigned int bps = modulation_types[OFDMFRAME_P_MOD].bps;
      unsigned int num_written;
      liquid_repack_bytes(payload_e,  8,  payload_enc_len,
                          payload_m, bps, payload_mod_len,
                          &num_written);
      frame_assembled = 1;
      #if DEBUG_OFDMFLEXFRAMEGEN
        printf("wrote %u symbols (expected %u)\n", num_written, _q->payload_mod_len);
      #endif
    }
    
    void modulator::encode_header()
    {
        // first 'n' bytes user data
        unsigned int n = OFDMFRAME_H_USR;
    
        // first byte is for expansion/version validation
        header_d[n + 0] = OFDMFRAME_VERSN;
    
        // add payload length
        header_d[n + 1] = (payload_dec_len >> 8) & 0xff;
        header_d[n + 2] = (payload_dec_len     ) & 0xff;
    
        // add modulation scheme/depth (pack into single byte)
        header_d[n + 3]  = OFDMFRAME_P_MOD;
    
        // add CRC, forward error-correction schemes
        //  CRC     : most-significant 3 bits of [n+4]
        //  fec0    : least-significant 5 bits of [n+4]
        //  fec1    : least-significant 5 bits of [n+5]
        header_d[n + 4]  = (OFDMFRAME_P_CRC & 0x07) << 5;
        header_d[n + 4] |= (OFDMFRAME_P_EC1) & 0x1f;
        header_d[n + 5]  = (OFDMFRAME_P_EC2) & 0x1f;
    
        // run packet encoder
        packetizer_encode(p_header, header_d, header_e);
    
        // scramble header
        scramble_data(header_e, OFDMFRAME_H_ENC);
    }
    
    void modulator::modulate_header()
    {
      unsigned int bps = modulation_types[OFDMFRAME_H_MOD].bps;
      unsigned int num_written;
      liquid_repack_bytes(header_e, 8, OFDMFRAME_H_ENC,
                          header_m, bps, OFDMFRAME_H_SYM,
                          &num_written);
    }
    
    // write first S0 symbol
    void modulator::write_S0a(std::complex<float> * _buffer)
    {
    #if DEBUG_OFDMFLEXFRAMEGEN
        printf("writing S0[a] symbol\n");
    #endif
        // write S0 symbol into front of buffer
        ofdmframegen_write_S0a(fg, _buffer);
        // update state
        state = STATE_S0B;
    }
    
    // write second S0 symbol
    void modulator::write_S0b(std::complex<float> * _buffer)
    {
    #if DEBUG_OFDMFLEXFRAMEGEN
        printf("writing S0[b] symbol\n");
    #endif
    
        // write S0 symbol into front of buffer
        ofdmframegen_write_S0b(fg, _buffer);
    
        // update state
        state = STATE_S1;
    }
    
    // write S1 symbol
    void modulator::write_S1(std::complex<float> * _buffer)
    {
    #if DEBUG_OFDMFLEXFRAMEGEN
        printf("writing S1 symbol\n");
    #endif
    
        // write S1 symbol into end of buffer
        ofdmframegen_write_S1(fg, _buffer);
    
        // update state
        symbol_number = 0;
        state = STATE_HDR;
    }
    
    // write header symbol
    void modulator::write_header(std::complex<float> * _buffer)
    {
    #if DEBUG_OFDMFRAMEGEN
      printf("writing header symbol\n");
    #endif
      // load data onto data subcarriers
      unsigned int i;
      int sctype;
      for (i = 0; i < M; i++) {
        sctype = p[i];
        if (sctype == OFDMFRAME_SCTYPE_DATA) {
          if (header_symbol_index < OFDMFRAME_H_SYM) {
            modem_modulate(mod_header,
                           header_m[header_symbol_index++],
                           X + i);
          }
          else {
            unsigned int sym = modem_gen_rand_sym(mod_header);
            modem_modulate(mod_header, sym, X + i);
          }
        }
        else {
          // ignore subcarrier (ofdmframegen handles nulls and pilots)
          X[i] = 0.0f;
        }
      }
      // write symbol
      ofdmframegen_writesymbol(fg, X, _buffer);
      // check state
      if(symbol_number == num_header_symbols) {
        symbol_number = 0;
        state = STATE_PLD;
      }
    }
    
    // write payload symbol
    void modulator::write_payload(std::complex<float> * _buffer)
    {
    #if DEBUG_OFDMFRAMEGEN
      printf("writing payload symbol\n");
    #endif
    
      // load data onto data subcarriers
      unsigned int i;
      int sctype;
      for (i = 0; i < M; i++) {
        sctype = p[i];
        if (sctype == OFDMFRAME_SCTYPE_DATA) {
          if (payload_symbol_index < payload_mod_len) {
            // modulate payload symbol onto data subcarrier
            modem_modulate(mod_payload,
                           payload_m[payload_symbol_index++],
                           X + i);
          } else {
            //printf("  random payload symbol\n");
            // load random symbol
            unsigned int sym = modem_gen_rand_sym(mod_payload);
            modem_modulate(mod_payload, sym, X + i);
          }
        } else {
          // ignore subcarrier (ofdmframegen handles nulls and pilots)
          X[i] = 0.0f;
        }
      }
      // write symbol
      ofdmframegen_writesymbol(fg, X, _buffer);
    
      // check to see if this is the last symbol in the payload
      if (symbol_number == num_payload_symbols)
          frame_complete = 1;
    }
    
    void modulator::print_config()
    {
      printf("modulator:\n");
      printf("    num subcarriers     :   %-u\n", M);
      printf("      * NULL            :   %-u\n", M_null);
      printf("      * pilot           :   %-u\n", M_pilot);
      printf("      * data            :   %-u\n", M_data);
      printf("    cyclic prefix len   :   %-u\n", cp_len);
      printf("    taper len           :   %-u\n", taper_len);
      printf("    properties:\n");
      printf("    frame assembled     :   %s\n", frame_assembled ? "yes" : "no");
      if (frame_assembled) {
        printf("    payload:\n");
        printf("      * decoded bytes   :   %-u\n", payload_dec_len);
        printf("      * encoded bytes   :   %-u\n", payload_enc_len);
        printf("      * modulated syms  :   %-u\n", payload_mod_len);
        printf("    total OFDM symbols  :   %-u\n", get_frame_len());
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
    
    void modulator::print_p()
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
  }
}
