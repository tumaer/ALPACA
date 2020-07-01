#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
import os.path
import h5py
import sys

def ErrorNormsAbsolute( cell_volume : np.array , solution : np.array , exact : np.array ) :
   assert cell_volume.size == solution.size == exact.size, "Array dimensions do not match. Norms cannot be computed"
   l1 = np.sum(   np.abs( solution - exact )          * cell_volume )
   l2 = np.sum( ( np.abs( solution - exact ) ** 2.0 ) * cell_volume ) ** 0.5
   li = np.max(   np.abs( solution - exact )                        )
   return np.array( [l1, l2, li], dtype = np.float64 )

def ErrorNormsRelative( cell_volume, solution, exact ) :
   assert cell_volume.size == solution.size == exact.size, "Array dimensions do not match. Norms cannot be computed"
   l1 = np.sum(   np.abs( solution / exact - 1.0 )       * cell_volume )
   l2 = np.sum( ( np.abs( solution / exact - 1.0 ) **2 ) * cell_volume ) ** 0.5
   li = np.max(   np.abs( solution / exact - 1.0 )                     )
   return np.array( [l1, l2, li], dtype = np.float64 )

def ComputeShockMachNumber( pressure_ratio : np.float64, sound_speed_ratio : np.float64, gamma : np.float64, initial_guess_mach_number : np.float64,\
                            termination_error : np.float64 = np.float64( 1.0e-15 ) ) :

   def dx(f, x) :
      return np.abs( 0.0 - f( x ) )

   def NewtonRaphson( f, df, x0, e ) :
      delta = dx( f, x0 )
      while delta > e:
         x0 = x0 - f( x0 ) / df( x0 )
         delta = dx( f, x0 )

      return x0

   gamma_one = np.float64( gamma - 1.0 )
   gamma_two = np.float64( gamma + 1.0 )
   gamma_thr = np.float64( 1.0 / gamma_one )

   def y( shock_mach_number ) :
      smnsq_one = shock_mach_number**2.0 - 1.0
      numerator = ( 1.0 - gamma_one / gamma_two * sound_speed_ratio * smnsq_one / shock_mach_number )**( 2.0 * gamma * gamma_thr )
      denominator = 1.0 + 2.0 * gamma / gamma_two * smnsq_one
      return  numerator / denominator - pressure_ratio

   def dy( shock_mach_number ) :
      smnsq_one = shock_mach_number**2.0 - 1.0
      cr_gam = sound_speed_ratio * gamma_one
      first_numerator = 2 * gamma_thr * gamma * ( ( 2.0 * sound_speed_ratio * gamma_one ) / gamma_two - ( sound_speed_ratio * gamma_one * smnsq_one ) / \
                        ( shock_mach_number**2 * gamma_two ) ) * ( 1.0 - ( sound_speed_ratio * gamma_one * smnsq_one ) / \
                        ( shock_mach_number * gamma_two ) )**( 2.0 * gamma_thr * gamma - 1.0 )
      first_denominator = ( 2.0 * gamma * smnsq_one ) / gamma_two + 1.0
      first_term = - first_numerator / first_denominator

      second_numerator = 4.0 * shock_mach_number * gamma * \
                         ( 1.0 - ( sound_speed_ratio * gamma_one * smnsq_one ) / ( shock_mach_number * gamma_two ) )**( 2.0 * gamma_thr * gamma )
      second_denominator = gamma_two * ( ( 2.0 * gamma * smnsq_one ) / gamma_two + 1.0 )**2
      second_term = - second_numerator / second_denominator

      return first_term + second_term

   return NewtonRaphson( y, dy, initial_guess_mach_number, termination_error )

def Riemann( x_cell_center : np.array , time : np.float64, gamma : np.float64, rho_left : np.float64, rho_right : np.float64, velocity_left : np.float64,\
             velocity_right : np.float64, pressure_left : np.float64, pressure_right : np.float64 ) :
   # Assumed structure of exact solution
   #
   #    \         /      |con |       |s|
   #     \   f   /       |tact|       |h|
   # left \  a  /  state |disc| state |o| right
   # state \ n /    2    |cont|   3   |c| state
   #   1    \ /          |tinu|       |k|   4
   #         |           |ity |       | |
   #      x = 0.5

   # Speeds of sound
   c_left =  np.float64( ( gamma * pressure_left  / rho_left  )**( 0.5 ) )
   c_right = np.float64( ( gamma * pressure_right / rho_right )**( 0.5 ) )

   pressure_ratio = pressure_right / pressure_left
   sound_speed_ratio = c_right / c_left

   # Call Newton's method to solve shock number
   shock_mach_number = ComputeShockMachNumber( pressure_ratio, sound_speed_ratio, gamma, np.float64( 1.0 ) )

   p34 = 1.0 + 2.0 * gamma / ( gamma + 1.0 ) * ( shock_mach_number**2.0 - 1.0 )
   p3 = p34 * pressure_right
   alpha = ( gamma + 1.0 ) / ( gamma - 1.0 )
   rho3 = rho_right * ( 1.0 + alpha * p34 ) / ( alpha + p34 )
   rho2 = rho_left * ( p34 * pressure_right / pressure_left )**( 1.0 / gamma )
   velocity2 = velocity_left - velocity_right + ( 2.0 / ( gamma - 1.0 ) ) * c_left *\
        ( 1.0 - ( p34 * pressure_right / pressure_left )**( ( gamma - 1.0 ) / ( 2.0 * gamma ) ) )
   c2 = ( gamma * p3 / rho2 )**( 0.5 )

	# Shock position
   shock_position = 0.5 + time * c_right * ( ( gamma - 1.0 ) / ( 2.0 * gamma ) + ( gamma + 1.0 ) / ( 2.0 * gamma ) * p34 )**0.5 + time * velocity_right
   # Position of contact discontinuity
   contact_position = 0.5 + velocity2 * time + time * velocity_right
   # Start of expansion fan
   expansion_start_position = 0.5 + ( velocity_left - c_left ) * time
   # End of expansion fan
   expansion_end_position = 0.5 + ( velocity2 + velocity_right - c2 ) * time

   p_exact = np.zeros( x_cell_center.size )
   velocity_exact = np.zeros( x_cell_center.size )
   rho_exact = np.zeros( x_cell_center.size )

   # Where to apply which condition
   left_condition   =       np.array( x_cell_center <  expansion_start_position, dtype = bool )
   fan_condition    = np.logical_and( x_cell_center >= expansion_start_position, x_cell_center < expansion_end_position )
   state2_condition = np.logical_and( x_cell_center >= expansion_end_position, x_cell_center < contact_position )
   state3_condition = np.logical_and( x_cell_center >= contact_position, x_cell_center < shock_position )
   right_condition  =       np.array( x_cell_center >= shock_position, dtype = bool )

   #Fill the Pressure
   p_exact = np.where( left_condition  , pressure_left, p_exact )
   p_exact = np.where( fan_condition   , pressure_left * ( 1.0 + ( expansion_start_position - x_cell_center ) / ( c_left * alpha * time ) )**( 2.0 * gamma / ( gamma - 1.0 ) ), p_exact )
   p_exact = np.where( state2_condition, p3, p_exact )
   p_exact = np.where( state3_condition, p3, p_exact )
   p_exact = np.where( right_condition , pressure_right, p_exact )

   #Fill the Velocities
   velocity_exact = np.where( left_condition  , velocity_left, velocity_exact )
   velocity_exact = np.where( fan_condition   , velocity_left + ( 2.0 / ( gamma + 1.0 ) )*( x_cell_center - expansion_start_position ) / time, velocity_exact )
   velocity_exact = np.where( state2_condition, velocity2 + velocity_right, velocity_exact )
   velocity_exact = np.where( state3_condition, velocity2 + velocity_right, velocity_exact )
   velocity_exact = np.where( right_condition , velocity_right, velocity_exact )

   #Fill the Densities
   rho_exact = np.where( left_condition  , rho_left, rho_exact )
   rho_exact = np.where( fan_condition   , rho_left * ( 1.0 + ( expansion_start_position - x_cell_center ) / ( c_left * alpha * time ) )**( 2.0 / ( gamma - 1.0 ) ), rho_exact )
   rho_exact = np.where( state2_condition, rho2, rho_exact )
   rho_exact = np.where( state3_condition, rho3, rho_exact )
   rho_exact = np.where( right_condition , rho_right, rho_exact )

   return [rho_exact, velocity_exact, p_exact]

def SodAnalysis( x_cell_center, density, velocity, volume, time ) :
   [rho_exact, velocity_exact, p_exact] = Riemann( x_cell_center, time, 1.4, 1.0, 0.125, 0.0, 0.0, 1.0, 0.1 )
   relative_rho = ErrorNormsRelative( volume, density, rho_exact )
   absolute_velocity = ErrorNormsAbsolute( volume, velocity, velocity_exact )
   return np.concatenate( ( relative_rho, absolute_velocity ) )

def ReadSimulationDataFromHdf( hdffile ) :
   h5file_data = h5py.File( hdffile, "r+" )
   cell_vertices = h5file_data["mesh_topology"]["cell_vertex_IDs"][:,:]
   vertex_coordinates = h5file_data["mesh_topology"]["cell_vertex_coordinates"][:,:]
   ordered_vertex_coordinates = vertex_coordinates[cell_vertices]
   cell_centers = np.mean( ordered_vertex_coordinates, axis = 1 )
   longest_axis = np.argmax( np.argmax( cell_centers, axis = 0 ) )
   x_cell_center = cell_centers[:, longest_axis]
   min_cell_coordinates = np.min( ordered_vertex_coordinates, axis = 1 )
   max_cell_coordinates = np.max( ordered_vertex_coordinates, axis = 1 )
   delta_xyz = max_cell_coordinates - min_cell_coordinates
   volume = np.prod( delta_xyz, axis = 1 )
   density = h5file_data["cell_data"]["density"][:,0,0]
   velocity = h5file_data["cell_data"]["velocity"][:,0,0]
   # Currently not used: pressure = h5file_data["cell_data"]["pressure"][:,0,0]

   return x_cell_center, density, velocity, volume,

def SodErrorsFromFile( hdffile, time ) :
   [x_cell_center, density, velocity, volume] = ReadSimulationDataFromHdf( hdffile )
   return SodAnalysis( x_cell_center, density, velocity, volume, time )

def main( arguments ) :

   hdffile = os.path.abspath( os.path.realpath( arguments.h5file ) )
   time = arguments.time

   [x_cell_center, density, velocity, volume] = ReadSimulationDataFromHdf( hdffile )
   res = SodAnalysis( x_cell_center, density, velocity, volume, time )

   np.set_printoptions( precision = 16, linewidth = 140 )
   print( "Computed Sod Error norms of file " + hdffile + " at time " + str( time ) )
   print( "L1-Rel-Rho, L2-Rel-rho, Li-Rel-Rho, L1-Abs-Vel, L2-Abs-Vel, Linf-Abs-Vel" )
   print( res )

   allowed_error = np.append( np.array( arguments.max_density_errors ), np.array( arguments.max_velocity_errors ) )
   if( ( res > allowed_error).any() ) :
      return 1
   else :
      return 0


def ParseArguments() :
   parser = ArgumentParser( description = "Compute the exact solution of the 1D sod shock tube problem" )
   parser.add_argument( "h5file", help = "The HDF5-file holding the data" )
   parser.add_argument( "time", help = "The time of the solution", type = np.float64 )
   parser.add_argument( "--max-density-errors", nargs = 3, default = np.array( [1000,1000,1000] , dtype = np.float64 ), help = "The maximum allowed density error in the three norms l1, l2, and linf", type = np.float64 )
   parser.add_argument( "--max-velocity-errors", nargs = 3, default = np.array( [1000,1000,1000], dtype = np.float64 ), help = "The maximum allowed density error in the three norms l1, l2, and linf", type = np.float64 )
   arguments = parser.parse_args()
   return arguments

if __name__ == "__main__" :
   arguments = ParseArguments()
   exit_code = main( arguments )
   sys.exit( exit_code )
