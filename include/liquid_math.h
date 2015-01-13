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

#ifndef LIQUID_MATH_H
#define LIQUID_MATH_H

#include <math.h>
#include <complex>
#include <liquid/liquid.h>

namespace liquid {
  namespace math {
    const std::complex<float> I(0.0, 1.0);
    const std::complex<float> Z(0.0, 0.0);
    std::complex<float> cexpjf(float theta);
  }
}

#endif
