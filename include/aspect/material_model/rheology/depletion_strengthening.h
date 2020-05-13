/*
  Copyright (C) 2019 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_material_model_rheology_depletion_strengthening_h
#define _aspect_material_model_rheology_depletion_strengthening_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {
      template <int dim>
      class DepletionStrengthening : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          DepletionStrengthening();

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * Compute depletion strengthening
           */

          double
          melt_fraction (const double temperature,
                         const double pressure,
                         const std::vector<double> &composition,
                         const Point<dim> &position) const;

          double
          compute_depl_effect (const double pressure,
                               const double temperature,
                               const std::vector<double> composition) const;

        private:

          /**
           * Melt fraction and depletion strengthening parameters
           */
          double A1;
          double A2;
          double A3;
          double B1;
          double B2;
          double B3;
          double C1;
          double C2;
          double C3;
          double r1;
          double r2;
          double beta;
          double M_cpx;
          double alpha_depletion;
          double delta_eta_depletion_max;
      };
    }
  }
}
#endif
