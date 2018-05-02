/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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



#include <aspect/heating_model/compositional_heating.h>
#include <aspect/utilities.h>
#include <aspect/global.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    CompositionalHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &/*material_model_outputs*/,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      for (unsigned int q = 0; q < heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // Compute compositional volume fractions
          const std::vector<double> volume_fractions = MaterialModel::compute_volume_fractions(material_model_inputs.composition[q],
                                                       fields_used_in_heat_production_averaging);

          // Calculate average compositional heat production
          double compositional_heat_production = 0.;
          for (unsigned int c=0; c < volume_fractions.size(); ++c)
            compositional_heat_production += volume_fractions[c] * heating_values[c];

          heating_model_outputs.heating_source_terms[q] = compositional_heat_production;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }


    template <int dim>
    void
    CompositionalHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Compositional heating");
        {
          prm.declare_entry("Compositional heating values","0",
                            Patterns::List (Patterns::Double(0)),
                            "List of heat production per unit volume values for "
                            "background and compositional fields, for a total of "
                            "N+1 values, where N is the number of compositional fields. "
                            "Units: $W/m^3$.");
          prm.declare_entry ("Use compositional field for heat production averaging", "1",
                             Patterns::List(Patterns::Integer(0,1)),
                             "List of integers, detailing for each compositional field if it should be included in the "
                             "averaging scheme when the heat production is computed (if 1) or not (if 0).");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    CompositionalHeating<dim>::parse_parameters (ParameterHandler &prm)
    {

      // increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Compositional heating");
        {
          std::vector<int> used_fields = Utilities::possibly_extend_from_1_to_N (
                                           Utilities::string_to_int(
                                             Utilities::split_string_list(
                                               prm.get("Use compositional field for heat production averaging"))),
                                           n_fields,
                                           "Use compositional field for heat production averaging");

          fields_used_in_heat_production_averaging.resize(used_fields.size());
          for (unsigned int i=0; i<used_fields.size(); ++i)
            fields_used_in_heat_production_averaging[i] = used_fields[i];

          heating_values = Utilities::possibly_extend_from_1_to_N (
                             Utilities::string_to_double(
                               Utilities::split_string_list(
                                 prm.get("Compositional heating values"))),
                             n_fields,
                             "Compositional heating values");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(CompositionalHeating,
                                  "compositional heating",
                                  "Implementation of a model in which magnitude of internal heat production "
                                  "is determined from fixed values assigned to each compositional "
                                  "field. These values are interpreted as having units $W/m^3$.")
  }
}


