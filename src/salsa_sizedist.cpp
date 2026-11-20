/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), mo_salsa_sizedist.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Number per bin of lognormal modes, for initial aerosol.
 */

#include <algorithm>
#include <cmath>

#include "salsa_sizedist.h"
#include "salsa_constants.h"


namespace Salsa
{
    // Fortran: mo_salsa_sizedist.f90: size_distribution
    // The modes are integrated over 10 equal diameter intervals per bin.
    std::vector<TF> lognormal_distribution(
            const std::vector<TF>& number,
            const std::vector<TF>& diameter,
            const std::vector<TF>& sigma,
            const Layout& layout,
            const int bin_start,
            const int bin_end)
    {
        const int n_intervals = 10;
        const TF minimum = 1.e-9_tf;

        std::vector<TF> number_bin(layout.aer_2b_end, 0._tf);

        for (int b=bin_start; b<bin_end; ++b)
        {
            const TF d_low = std::pow(layout.volume_low[b]/pi6, 1._tf/3._tf);
            const TF d_high = std::pow(layout.volume_high[b]/pi6, 1._tf/3._tf);
            const TF delta_d = (d_high - d_low)/n_intervals;

            for (int i=1; i<=n_intervals; ++i)
            {
                const TF d1 = d_low + (static_cast<TF>(i) - 1._tf)*delta_d;
                const TF d2 = d1 + delta_d;
                const TF d_mid = (d1 + d2)/2._tf;
                const TF delta_log_d = std::log(d2/d1);

                TF sum = 0._tf;
                for (size_t m=0; m<number.size(); ++m)
                {
                    const TF log_sigma = std::log(std::max(sigma[m], minimum));
                    const TF log_ratio = std::log(d_mid/std::max(diameter[m], minimum));

                    sum += number[m]*delta_log_d/(std::sqrt(2._tf*pi)*log_sigma)
                         * std::exp(-(log_ratio*log_ratio)/(2._tf*(log_sigma*log_sigma)));
                }

                number_bin[b] += sum;
            }
        }

        return number_bin;
    }
}
