//===----------------------- stl_utilities.cpp ----------------------------===//
//
//                                 ALPACA
//
// Part of ALPACA, under the GNU General Public License as published by
// the Free Software Foundation version 3.
// SPDX-License-Identifier: GPL-3.0-only
//
// If using this code in an academic setting, please cite the following:
// @article{hoppe2022parallel,
//  title={A parallel modular computing environment for three-dimensional
//  multiresolution simulations of compressible flows},
//  author={Hoppe, Nils and Adami, Stefan and Adams, Nikolaus A},
//  journal={Computer Methods in Applied Mechanics and Engineering},
//  volume={391},
//  pages={114486},
//  year={2022},
//  publisher={Elsevier}
// }
//
//===----------------------------------------------------------------------===//
#include "input_output/utilities/stl_utilities.h"

#include "utilities/mathematical_functions.h"

namespace StlUtilities {
/**
 * @brief Reads out and returns a double array from the STL file denoting either
 * a triangle point or normal.
 * @param stl_file_stream The stl file stream.
 * @return The array of read STL triangle.
 */
std::array<double, 3> ReadStlArray(std::ifstream &stl_file_stream) {
  double const x = static_cast<double>(ReadStlValue<float>(stl_file_stream));
  double const y = static_cast<double>(ReadStlValue<float>(stl_file_stream));
  double const z = static_cast<double>(ReadStlValue<float>(stl_file_stream));
  std::array<double, 3> const return_array({x, y, z});
  return return_array;
}

/**
 * @brief Reads out an STL file and returns geometry information. See \cite STL
 * for specifications of STL files.
 * @param stl_filename The name of the STL file.
 * @return The vector of read STL triangles.
 */
std::vector<Triangle> ReadStl(std::string const &stl_filename) {
  std::ifstream stl_file_stream(stl_filename.c_str(),
                                std::ios::in | std::ios::binary);
  if (!stl_file_stream) {
    throw std::invalid_argument("Specfied STL file does not exist: " +
                                stl_filename);
  }

  // shift by 80 bytes which constitute the header by STL standard
  stl_file_stream.seekg(80, std::ios::beg);

  // read the number of triangles (following four bytes)
  unsigned int const number_of_triangles =
      ReadStlValue<unsigned int>(stl_file_stream);

  std::vector<Triangle> return_vector;
  return_vector.reserve(number_of_triangles);

  // go through file and store triangle information in vector
  for (unsigned int i = 0; i < number_of_triangles; i++) {
    std::array<double, 3> const normal = ReadStlArray(stl_file_stream);
    std::array<double, 3> const v1 = ReadStlArray(stl_file_stream);
    std::array<double, 3> const v2 = ReadStlArray(stl_file_stream);
    std::array<double, 3> const v3 = ReadStlArray(stl_file_stream);
    return_vector.emplace_back(Triangle(normal, v1, v2, v3));
    // move by two bytes which are irrelevant by STL standard
    stl_file_stream.seekg(2, std::ios::cur);
  }
  return return_vector;
}

/**
 * @brief Computes the signed distance from the triangle to the grid point in
 * case the projected point is outside of the triangle.
 * @param point The point for which the distance is computed.
 * @param virtual_point_1 Projection of the point into the triangle plane.
 * @param start_to_virtual_1 Vector from the start of the closest edge to the
 * virtual point.
 * @param end_to_virtual_1 Vector from the edge of the closest edge to the
 * virtual point.
 * @param edge Vector of the closest edge.
 * @param start Start of the closest edge.
 * @param end End of the closest edge.
 * @param p1_to_virtual_1 Vector from point to the virtual point.
 * @param signed_distance_inside Distance between point and virtual point.
 */
double DetermineOutsideDistance(std::array<double, 3> const &point,
                                std::array<double, 3> const &virtual_point_1,
                                std::array<double, 3> const &start_to_virtual_1,
                                std::array<double, 3> const &end_to_virtual_1,
                                std::array<double, 3> const &edge,
                                std::array<double, 3> const &start,
                                std::array<double, 3> const &end,
                                std::array<double, 3> const &p1_to_virtual_1,
                                double const &signed_distance_inside) {
  std::array<double, 3> const R =
      VU::VectorTripleProduct(start_to_virtual_1, end_to_virtual_1, edge);
  double const one_R_abs = 1.0 / VU::L2Norm(R);
  double const distance_virtual_points =
      VU::DotProduct(p1_to_virtual_1, R) * one_R_abs;
  // Find the absolute maximum value of the edge
  auto const index = std::max_element(std::cbegin(edge), std::cend(edge),
                                      [](double const a, double const b) {
                                        return std::abs(a) < std::abs(b);
                                      }) -
                     std::cbegin(edge);
  // Compute the location
  double const location =
      (virtual_point_1[index] - distance_virtual_points * R[index] * one_R_abs -
       start[index]) /
      edge[index];

  if (location >
      1.0) { // projected point behind line, thereofore max point towards end
    return VU::Distance(point, end);
  } else if (location < 0.0) {
    return VU::Distance(point, start);
  } else {
    return std::sqrt(distance_virtual_points * distance_virtual_points +
                     signed_distance_inside * signed_distance_inside);
  }
}

/**
 * @brief Computes the signed distance from the triangle to the grid point and
 * writes in the implicit return parameter. See \cite Jones1995.
 * @param triangle The triangle for which the distance is computed.
 * @param x/y/z The coordinates of the grid point.
 * @param signed_distance Implicit return parameter for the signed distance.
 */
void Voxelization(Triangle const &triangle, std::array<double, 3> const &point,
                  double &signed_distance) {
  // Gets all edges from the triangle and the inverse counter part of it
  std::array<std::array<double, 3>, 3> const edges = {
      VU::Difference(triangle.p1, triangle.p2),
      VU::Difference(triangle.p2, triangle.p3),
      VU::Difference(triangle.p3, triangle.p1)};
  std::array<std::array<double, 3>, 3> const inverse_edges = {
      VU::Multiply(edges[0], -1.0), VU::Multiply(edges[1], -1.0),
      VU::Multiply(edges[2], -1.0)};

  // Computes the inner distance between the point and the triangle corner point
  double const signed_distance_inside =
      VU::DotProduct(VU::Difference(triangle.p1, point), triangle.normal);
  std::array<double, 3> const virtual_point_1(VU::Difference(
      VU::Multiply(triangle.normal, signed_distance_inside), point));

  std::array<double, 3> const P1virtualP0(
      VU::Difference(triangle.p1, virtual_point_1));
  std::array<double, 3> const P2virtualP0(
      VU::Difference(triangle.p2, virtual_point_1));
  std::array<double, 3> const P3virtualP0(
      VU::Difference(triangle.p3, virtual_point_1));

  double const f1 =
      VU::ScalarTripleProduct(triangle.v[0], P1virtualP0, triangle.normal);
  double const f2 =
      VU::ScalarTripleProduct(triangle.v[1], P2virtualP0, triangle.normal);
  double const f3 =
      VU::ScalarTripleProduct(triangle.v[2], P3virtualP0, triangle.normal);

  bool inside_point = false;
  double new_signed_distance = signed_distance_inside;
  // anticlockwise if > 0, clockwise G < 0
  if (f2 <= 0 && f1 > 0) { // clockwise v2, anticlockwise v1
    inside_point = VU::ScalarTripleProduct(P1virtualP0, P2virtualP0,
                                           triangle.normal) >= 0.0
                       ? true
                       : false;
    if (!inside_point) {
      new_signed_distance = DetermineOutsideDistance(
          point, virtual_point_1, P2virtualP0, P1virtualP0, edges[0],
          triangle.p1, triangle.p2, P1virtualP0, signed_distance_inside);
    }
  } else if (f1 <= 0 && f3 > 0) { // clockwise v1, anticlockwise v3
    inside_point = VU::ScalarTripleProduct(P3virtualP0, P1virtualP0,
                                           triangle.normal) >= 0.0
                       ? true
                       : false;
    if (!inside_point) {
      new_signed_distance = DetermineOutsideDistance(
          point, virtual_point_1, P1virtualP0, P3virtualP0, inverse_edges[2],
          triangle.p1, triangle.p3, P1virtualP0, signed_distance_inside);
    }
  } else if (f3 <= 0 && f2 > 0) { // clockwise v3, anticlockwise v2
    inside_point = VU::ScalarTripleProduct(P2virtualP0, P3virtualP0,
                                           triangle.normal) >= 0.0
                       ? true
                       : false;
    if (!inside_point) {
      new_signed_distance = DetermineOutsideDistance(
          point, virtual_point_1, P2virtualP0, P3virtualP0, edges[1],
          triangle.p2, triangle.p3, P1virtualP0, signed_distance_inside);
    }
  } else {
    throw std::logic_error("Wrong triangle distance in STL voxelization.");
  }

  new_signed_distance =
      std::abs(new_signed_distance) * Sign(signed_distance_inside);

  if (std::abs(new_signed_distance) < std::abs(signed_distance) ||
      (std::abs(new_signed_distance) == std::abs(signed_distance) &&
       inside_point)) {
    signed_distance = new_signed_distance;
  }
}
} // namespace StlUtilities
