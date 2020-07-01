#!/usr/bin/env python3
import os
import numpy as np
import h5py
import sympy
from argparse import ArgumentParser as argparser

def main( h5file ) :

   with h5py.File( h5file, "r") as data:
      density = np.array(data["cell_data"]["density"][:,0,0])
      cell_vertices = np.array(data["mesh_topology"]["cell_vertex_IDs"])
      vertex_coordinates = np.array(data["mesh_topology"]["cell_vertex_coordinates"])

   # Number of cells per dimension:
   # Only meaningful for cubic setups, otherwise full symmetry cannot be expected
   nc, is_integer = sympy.integer_nthroot( density.shape[0], 3)
   if( is_integer == False ) :
       log = ['ERROR: Domain has no cubic size']
       return ( False, log )

   ordered_vertex_coordinates = vertex_coordinates[cell_vertices]
   coords = np.mean( ordered_vertex_coordinates, axis = 1 )

   first_trafo = coords[:,0].argsort(kind='stable')
   coords = coords[first_trafo]
   second_trafo = coords[:,1].argsort(kind='stable')
   coords = coords[second_trafo]
   third_trafo = coords[:,2].argsort(kind='stable')
   coords = coords[third_trafo]

   trafo = first_trafo[second_trafo[third_trafo]]

   # Corner points
   log = ['---Test Symmetry at 8 Corner Points---']

   density = density[trafo]
   density = density.reshape( nc, nc, nc )
   coners_are_symmetric = (     density[0][0][0] == density[-1][ 0][ 0]
                            and density[0][0][0] == density[ 0][-1][ 0]
                            and density[0][0][0] == density[-1][-1][ 0]
                            and density[0][0][0] == density[ 0][ 0][-1]
                            and density[0][0][0] == density[-1][ 0][-1]
                            and density[0][0][0] == density[ 0][-1][-1]
                            and density[0][0][0] == density[-1][-1][-1]
                           )

   log.append( "Symmetry at corner points is fully recovered: " + str( coners_are_symmetric ) )

   # Points for full symmetry
   log.append( '---Test Full Symmetry at 48 Arbitrary Points---' )

   point1 = 3
   point2 = 7
   point3 = 1
   point4 = -2
   point5 = -8
   point6 = -4

   d1  = density[point1][point2][point3]
   d2  = density[point1][point2][point4]
   d3  = density[point1][point5][point3]
   d4  = density[point1][point5][point4]
   d5  = density[point1][point3][point2]
   d6  = density[point1][point3][point5]
   d7  = density[point1][point4][point2]
   d8  = density[point1][point4][point5]

   d9  = density[point2][point1][point3]
   d10 = density[point2][point1][point4]
   d11 = density[point2][point6][point3]
   d12 = density[point2][point6][point4]
   d13 = density[point2][point3][point1]
   d14 = density[point2][point3][point6]
   d15 = density[point2][point4][point1]
   d16 = density[point2][point4][point6]

   d17 = density[point3][point1][point2]
   d18 = density[point3][point1][point5]
   d19 = density[point3][point6][point2]
   d20 = density[point3][point6][point5]
   d21 = density[point3][point2][point1]
   d22 = density[point3][point2][point6]
   d23 = density[point3][point5][point1]
   d24 = density[point3][point5][point6]

   d25 = density[point4][point1][point2]
   d26 = density[point4][point1][point5]
   d27 = density[point4][point6][point2]
   d28 = density[point4][point6][point5]
   d29 = density[point4][point2][point1]
   d30 = density[point4][point2][point6]
   d31 = density[point4][point5][point1]
   d32 = density[point4][point5][point6]

   d33 = density[point5][point1][point3]
   d34 = density[point5][point1][point4]
   d35 = density[point5][point6][point3]
   d36 = density[point5][point6][point4]
   d37 = density[point5][point3][point1]
   d38 = density[point5][point3][point6]
   d39 = density[point5][point4][point1]
   d40 = density[point5][point4][point6]

   d41 = density[point6][point2][point3]
   d42 = density[point6][point2][point4]
   d43 = density[point6][point5][point3]
   d44 = density[point6][point5][point4]
   d45 = density[point6][point3][point2]
   d46 = density[point6][point3][point5]
   d47 = density[point6][point4][point2]
   d48 = density[point6][point4][point5]

   full_symm = d1 == d2  and d1 == d3  and d1 == d4  and d1 == d5  and d1 == d6  and d1 == d7  and d1 == d8  and d1 == d9  and d1 == d10 and \
               d1 == d11 and d1 == d12 and d1 == d13 and d1 == d14 and d1 == d15 and d1 == d16 and d1 == d17 and d1 == d18 and d1 == d19 and \
               d1 == d20 and d1 == d21 and d1 == d22 and d1 == d23 and d1 == d24 and d1 == d25 and d1 == d26 and d1 == d27 and d1 == d28 and \
               d1 == d29 and d1 == d30 and d1 == d31 and d1 == d32 and d1 == d33 and d1 == d34 and d1 == d35 and d1 == d36 and d1 == d37 and \
               d1 == d38 and d1 == d39 and d1 == d40 and d1 == d41 and d1 == d42 and d1 == d43 and d1 == d44 and d1 == d45 and d1 == d46 and \
               d1 == d47 and d1 == d48

   log.append( 'Full Domain Symmetry is fully recovered: ' + str( full_symm ) )

   if( full_symm == False ) :
      log.append( 'The actual values are: ' )
      log.append( d1 )
      log.append( d2 )
      log.append( d3 )
      log.append( d4 )
      log.append( d5 )
      log.append( d6 )
      log.append( d7 )
      log.append( d8 )

      log.append( d9 )
      log.append( d10 )
      log.append( d11 )
      log.append( d12 )
      log.append( d13 )
      log.append( d14 )
      log.append( d15 )
      log.append( d16 )

      log.append( d17 )
      log.append( d18 )
      log.append( d19 )
      log.append( d20 )
      log.append( d21 )
      log.append( d22 )
      log.append( d23 )
      log.append( d24 )

      log.append( d25 )
      log.append( d26 )
      log.append( d27 )
      log.append( d28 )
      log.append( d29 )
      log.append( d30 )
      log.append( d31 )
      log.append( d32 )

      log.append( d33 )
      log.append( d34 )
      log.append( d35 )
      log.append( d36 )
      log.append( d37 )
      log.append( d38 )
      log.append( d39 )
      log.append( d40 )

      log.append( d41 )
      log.append( d42 )
      log.append( d43 )
      log.append( d44 )
      log.append( d45 )
      log.append( d46 )
      log.append( d47 )
      log.append( d48 )

   return ( coners_are_symmetric and full_symm, log )

def CheckSymmetry( h5file ) :
   return main( h5file )

def ParseArguments() :
   parser = argparser( description = "Automatic Check for Floating-Point Symmertry Preservation in 3D Simulations. Therefore a results file of simulating \
                                      the three dimensional implosion case described in Fleischmann et al., 2019, 'Numerical symmetry-preserving techniques for low-dissipation shock-capturing schemes' \
                                      is analysed." )
   parser.add_argument( "hdffile", help = "The HDF5 file to be checked for symmetry" )
   arguments = parser.parse_args()
   h5file = os.path.abspath( arguments.hdffile )
   return h5file

if __name__ == "__main__":
   result, log = main( ParseArguments() )
   for line in log :
      print( line )
