/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__utilities_h
#define __aspect__utilities_h

#include <aspect/global.h>

#include <deal.II/base/std_cxx1x/array.h>
#include <deal.II/base/point.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>

#include <aspect/geometry_model/interface.h>


namespace aspect
{
  /**
   * A namespace for utility functions that might be used in many different
   * places to prevent code duplication.
   */
  namespace Utilities
  {
    using namespace dealii;

    /**
     * Returns spherical coordinates of a cartesian point. The returned array
     * is filled with radius, phi and theta (polar angle). If the dimension is
     * set to 2 theta is omitted. Phi is always normalized to [0,2*pi].
     *
     */
    template <int dim>
    std_cxx1x::array<double,dim>
    spherical_coordinates(const Point<dim> &position);

    /**
     * Return the cartesian point of a spherical position defined by radius,
     * phi and theta (polar angle). If the dimension is set to 2 theta is
     * omitted.
     */
    template <int dim>
    Point<dim>
    cartesian_coordinates(const std_cxx1x::array<double,dim> &scoord);


    template <int dim, int grid_dim>
    class AsciiDataLookup
    {
      public:
        AsciiDataLookup(const GeometryModel::Interface<dim> &geometry_model,
                        const unsigned int components,
                        const double scale_factor,
                        const types::boundary_id boundary_id = numbers::invalid_boundary_id);

        /**
         * Checks whether a file named filename exists.
         *
         * @param filename File to check existence
         */
        bool fexists(const std::string &filename);

        /**
         * Outputs the AsciiData module information at model start.
         */
        void screen_output(const ConditionalOStream &pcout) const;

        /**
         * Loads a data text file. Throws an exception if the file does not exist,
         * if the data file format is incorrect or if the file grid changes over model runtime.
         */
        void
        load_file(const std::string &filename);

        /**
         * Returns the computed data (velocity, temperature, etc. - according to the used plugin)
         * in cartesian coordinates.
         *
         * @param position The current position to compute the data (velocity, temperature, etc.)
         * @param time_weight A weighting between the two current timesteps n and n+1
         */
        double
        get_data(const Point<dim> &position,
                 const unsigned int component,
                 const double time_weight) const;

      private:
        /**
         * The number of data components read in (=columns in the data file).
         */
        const unsigned int components;

        /**
         * A reference to the geometry model. Is needed to convert the
         * position into spherical coordinates if necessary.
         */
        const GeometryModel::Interface<dim> &geometry_model;

        /**
         * Interpolation functions to access the data.
         */
        std::vector<Functions::InterpolatedUniformGridData<grid_dim> *> data;
        std::vector<Functions::InterpolatedUniformGridData<grid_dim> *> old_data;

        /**
         * Model size
         */
        std_cxx11::array<std::pair<double,double>,grid_dim> grid_extent;

        /**
         * Number of points in the data grid.
         */
        TableIndices<grid_dim> table_points;

        /**
         * Dimensions of the boundary plane
         */
        unsigned int boundary_dimensions[grid_dim];

        /**
         * Scales the data boundary condition by a scalar factor. Can be
         * used to transform the unit of the data.
         */
        const double scale_factor;


        /**
         * Gets the extents of the model in the relevant dimensions and returns
         * the according minimum and maximum value of the boundary in each dimension.
         * In case of a spherical shell the function returns the values in
         * spherical coordinates.
         */
        std_cxx11::array<std::pair<double,double>,dim>
        get_model_extent (const GeometryModel::Interface<dim> &geometry_model) const;


        /**
         * Determines which of the dimensions of the position is used to find
         * the data point in the data grid. E.g. the left boundary of a box model extents in
         * the y and z direction (position[1] and position[2]), therefore the function
         * would return [1,2] for dim==3 or [1] for dim==2. We are lucky that these indices are
         * identical for the box and the spherical shell (if we use spherical coordinates for the
         * spherical shell), therefore we do not need to distinguish between them. For the initial
         * condition this function is trivial, because the position in the data grid is the same as
         * the actual position (the function returns [0,1,2] or [0,1]), but for the boundary
         * conditions it matters.
         */
        std_cxx11::array<unsigned int,grid_dim>
        get_boundary_dimensions (const types::boundary_id boundary_id) const;


        /**
         * Computes the table indices of each entry in the input data file.
         * The index depends on dim, grid_dim and the number of components.
         */
        TableIndices<grid_dim>
        compute_table_indices(const unsigned int i) const;

    };


  }
}

#endif
