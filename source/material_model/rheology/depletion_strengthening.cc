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


#include <aspect/material_model/rheology/depletion_strengthening.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      DepletionStrengthening<dim>::DepletionStrengthening ()
      {}



      template <int dim>
      double
      DepletionStrengthening<dim>::compute_depl_effect (const double pressure,
                                              const double temperature,
                                              const unsigned int composition) const
      {
        // INFO HERE
        
        if (this->introspection().compositional_name_exists("maximum_melt_fraction"))
           {    
             // extract depletion = maximum experienced melt fraction:
             const unsigned int melt_index = this->introspection().compositional_index_for_name("maximum_melt_fraction");
             const double depletion_visc = std::min(1.0, std::max(composition[melt_index],0.0));

             // calculate strengthening due to depletion:
             const double depletion_strengthening = std::min(exp(alpha_depletion*depletion_visc),delta_eta_depletion_max);
           }
        else
           {    
             const double depletion_strengthening = 1.0;
           }

        return depletion_strengthening;
      }


      template <int dim>
      void
      DepletionStrengthening<dim>::declare_parameters (ParameterHandler &prm)
      {
          // Mantle melting parameterization following notation of Katz et al. (2003)
          prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: $°C$.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: $°C$.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: $°C$.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Exponential depletion strengthening factor", "0.0",
                             Patterns::Double (0),
                             "$\\alpha_F$: exponential dependency of viscosity on the depletion "
                             "field $F$ (called peridotite). "
                             "Dimensionless factor. With a value of 0.0 (the default) the "
                             "viscosity does not depend on the depletion. The effective viscosity increase"
                             "due to depletion is defined as $exp( \\alpha_F * F)$. "
                             "Rationale: melting dehydrates the source rock by removing most of the volatiles,"
                             "and makes it stronger. Hirth and Kohlstedt (1996) report typical values around a "
                             "factor 100 to 1000 viscosity contrast between wet and dry rocks, although some "
                             "experimental studies report a smaller (factor 10) contrast (e.g. Fei et al., 2013).");
          prm.declare_entry ("Maximum Depletion viscosity change", "1.0e3",
                             Patterns::Double (0),
                             "$\\Delta \\eta_{F,max}$: maximum depletion strengthening of viscosity. "
                             "Rationale: melting dehydrates the source rock by removing most of the volatiles,"
                             "and makes it stronger. Hirth and Kohlstedt (1996) report typical values around a "
                             "factor 100 to 1000 viscosity contrast between wet and dry rocks, although some "
                             "experimental studies report a smaller (factor 10) contrast (e.g. Fei et al., 2013).");
      }

      template <int dim>
      void
      DiffusionCreep<dim>::parse_parameters (ParameterHandler &prm)
      {
          // Mantle melting parameterization following notation of Katz et al. (2003)
          A1                                = prm.get_double ("A1");
          A2                                = prm.get_double ("A2");
          A3                                = prm.get_double ("A3");
          B1                                = prm.get_double ("B1");
          B2                                = prm.get_double ("B2");
          B3                                = prm.get_double ("B3");
          C1                                = prm.get_double ("C1");
          C2                                = prm.get_double ("C2");
          C3                                = prm.get_double ("C3");
          r1                                = prm.get_double ("r1");
          r2                                = prm.get_double ("r2");
          beta                              = prm.get_double ("beta");
          M_cpx                             = prm.get_double ("Mass fraction cpx");
          alpha_depletion                   = prm.get_double ("Exponential depletion strengthening factor");
          delta_eta_depletion_max           = prm.get_double ("Maximum Depletion viscosity change");
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class DepletionStrengthening<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
