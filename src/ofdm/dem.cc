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

namespace liquid {
  namespace ofdm {
    demodulator::demodulator(unsigned int       _M,
                       unsigned int       _cp_len,
                       unsigned int       _taper_len,
                       unsigned char *    _p,
                       framesync_callback _callback,
                       void *             _userdata)
    {
      M         = _M;
      cp_len    = _cp_len;
      taper_len = _taper_len;
      p         = _p;
      callback  = _callback;
      userdata  = _userdata;
      fs = ofdmflexframesync_create(M,
                                    cp_len,
                                    taper_len,
                                    p,
                                    callback,
                                    userdata);
      reset();
    }
    
    demodulator::~demodulator()
    {
      ofdmflexframesync_destroy(fs);
    }
    
    void demodulator::reset()
    {
      ofdmflexframesync_reset(fs);
    }
    
    unsigned int demodulator::get_M()
    {
      return M;
    }
    
    void demodulator::set_payload_len(unsigned int _payload_len)
    {
      payload_len = _payload_len;
    }
  }
}
