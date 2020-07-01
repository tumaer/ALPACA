#!/usr/bin/env python3

import os
import h5py
import numpy as np
import xml.etree.ElementTree as et
from argparse import ArgumentParser


def SanityCheck( use_geometric_res, use_levelset_res ) :
   assert np.array( use_geometric_res + [use_levelset_res] ).any(), "Neither any geometric nor the levelset restriction is set to True, why bother running this script?"

def CylindricalRestriciton( plane, radius, location, coords ) :
   array_circle = np.ones( ( len( coords ), 2 ) ) * np.array( location )
   distance = np.linalg.norm( array_circle - coords[:, [int( i ) - 1 for i in plane]], axis = 1 )
   indices_to_keep_circle = np.where( distance <= radius )[0]
   return indices_to_keep_circle

def SphericalRestriction( radius, location, coords ) :
   array_sphere = np.ones( ( len( coords ), 3 ) ) * np.array( location )
   distance = np.linalg.norm( array_sphere - coords, axis = 1 )
   indices_to_keep_sphere = np.where( distance <= radius )[0]
   return indices_to_keep_sphere

def LineRestriction( axis, location, coords ):
   indices_to_keep_ax = []
   for ax in axis:
     indices_to_keep_ax.append( np.where( ( coords[:,int( ax ) - 1] > location[int( ax ) - 1][0]) & (coords[:, int( ax ) - 1] < location[int( ax ) - 1][1] ) )[0] )
     coords = coords[indices_to_keep_ax[-1]]
   indices_to_keep_line = TraceBackIndices( indices_to_keep_ax )
   return indices_to_keep_line

def TraceBackIndices( itk_list ):
   if( len( itk_list ) > 0 ) :
      indices_to_keep = itk_list[-1]
      for i in reversed( range( len( itk_list ) - 1 ) ):
         indices_to_keep = itk_list[i][indices_to_keep]
      return indices_to_keep
   else :
      return 0

def TrimFiles( path_data, path_trimmed, cut_off_for_levelset, use_levelset_res, use_geometric_res, cyl_res,
               sph_res, line_res, quantities_to_keep, precision, rearange_vertex_coordinates ) :

   if not os.path.exists( path_trimmed ):
      os.mkdir( path_trimmed )

   files_h5 = []
   files_xdmf = []
   for file in os.listdir( path_data ):
      if file.endswith( "h5" ):
         files_h5.append( file )
      if file.endswith( "xdmf" ):
         if "data" in file and not "time_series" in file:
            files_xdmf.append( file )

   print( files_h5, files_xdmf )

   # SORTING FILES
   sort_h5 = []
   sort_xdmf = []
   for h5, xdmf in zip( files_h5, files_xdmf ):
      sort_h5.append( float( h5[5:13] ) )
      sort_xdmf.append( float( xdmf[5:13] ) )
   sort_h5 = np.array( sort_h5 )
   sort_xdmf = np.array( sort_xdmf )
   indices_h5 = np.argsort( sort_h5 )
   indices_xdmf = np.argsort( sort_xdmf )
   files_h5 = list( np.array( files_h5 )[indices_h5] )
   files_xdmf = list( np.array( files_xdmf )[indices_xdmf] )

   trimmed_data = {}
   for filename_h5, filename_xdmf in zip( files_h5, files_xdmf ):

      print( "\nREADING DATA FROM ", os.path.join( path_data, filename_h5 ) )
      h5file_data = h5py.File( os.path.join( path_data, filename_h5 ), "r+" )
      cell_vertices = h5file_data["mesh_topology"]["cell_vertex_IDs"][:,:]
      vertex_coordinates = h5file_data["mesh_topology"]["cell_vertex_coordinates"][:,:]
      if np.array( use_geometric_res ).any():
         coords = np.mean( vertex_coordinates[cell_vertices], axis = 1 )

      print( "FINDING INDICES TO KEEP" )
      if np.array(use_geometric_res).any():
         itk_geometric = []
         # LINE RESTRICTION
         if use_geometric_res[0] :
            itk_line = LineRestriction( line_res[0], line_res[1], coords )
            coords = coords[itk_line]
            itk_geometric.append(itk_line)

         # SPHERICAL RESTRICTION
         if use_geometric_res[1] :
            itk_cylinder = CylindricalRestriciton( cyl_res[0], cyl_res[1], cyl_res[2], coords )
            coords = coords[itk_cylinder]
            itk_geometric.append( itk_cylinder )

         # CYLINDRICAL RESTRICTION
         if use_geometric_res[2] :
            itk_sphere = SphericalRestriction( sph_res[0], sph_res[1], coords )
            itk_geometric.append( itk_sphere )

         itk_geometric = TraceBackIndices(itk_geometric)

      # LEVELSET RESTRICTION
      if use_levelset_res :
         if np.array( use_geometric_res ).any() :
            levelset = h5file_data["cell_data"]["levelset"][:,0,0]
            ITK_levelset = np.where( ( levelset[itk_geometric] >= cut_off_for_levelset[0] ) & ( levelset[itk_geometric] < cut_off_for_levelset[1] ) )[0]
            indices_to_keep = itk_geometric[ITK_levelset]
         else :
            levelset = h5file_data["cell_data"]["levelset"][:,0,0]
            indices_to_keep = np.where( ( levelset >= cut_off_for_levelset[0] ) & ( levelset < cut_off_for_levelset[1] ) )[0]
      else :
         indices_to_keep = itk_geometric

      cell_vertices = cell_vertices[indices_to_keep]

      # REARANGING VERTEX COORDINATES
      if rearange_vertex_coordinates :
         remaining_vertices = np.unique( cell_vertices )
         vertex_coordinates = vertex_coordinates[remaining_vertices]
         remaining_vertices = dict( zip( remaining_vertices, range( len( remaining_vertices ) ) ) )
         for j in range( cell_vertices.shape[1] ) :
            for i in range( cell_vertices.shape[0] ) :
               cell_vertices[i,j] = remaining_vertices[cell_vertices[i,j]]
            print( "Finished rearranging vertices with index %d of cell_vertices" % j )

      # WRITING NEW HDF5 FILE
      print( "WRITING TRIMMED HDF5 FILE ", os.path.join(path_trimmed, filename_h5 ) )
      h5file_trimmed = h5py.File( os.path.join(path_trimmed, filename_h5), "w" )
      h5file_trimmed.create_group( "mesh_topology" )
      h5file_trimmed.create_group( "cell_data" )
      h5file_trimmed["mesh_topology"].create_dataset( name = "cell_vertex_IDs", data = cell_vertices, dtype = precision[0] )
      h5file_trimmed["mesh_topology"].create_dataset( name = "cell_vertex_coordinates", data = vertex_coordinates, dtype = precision[1] )
      for quantity in h5file_data["cell_data"]:
         if quantity in quantities_to_keep:
            h5file_trimmed["cell_data"].create_dataset( name = quantity, data = h5file_data["cell_data"][quantity][indices_to_keep,...],
                                              dtype = precision[2] )

      h5file_data.close()
      h5file_trimmed.close()

      # WRITING NEW XDMF FILE
      no_cells = len( cell_vertices )
      no_vertices = len( vertex_coordinates )

      tree = et.parse( os.path.join( path_data, filename_xdmf ) )
      root = tree.getroot()

      for topology in root.iter( "Topology" ) :
         topology.set( "NumberOfElements", str( no_cells ) )
         topology[0].set( "Dimensions", str( no_cells ) + " 8" )

      for topology in root.iter( "Geometry" ) :
         topology.set( "NumberOfElements", str( no_vertices ) )
         topology[0].set( "Dimensions", str( no_vertices ) + " 3" )
         topology[0].set( "Precision", str( precision[1] ) )

      for quantity in root[0][0][0].findall( "Attribute" ) :
         # get the dimensions of the original file and add the new number of cells
         dimensions = quantity[0].attrib["Dimensions"].split(' ')
         dimensions[0] = str( no_cells )
         dimensions = ' '.join( dimensions )
         # Remove quantity if not desired or change dimension to restricted data
         if quantity.attrib["Name"] not in quantities_to_keep :
            root[0][0][0].remove( quantity )
         else:
            quantity[0].set( "Dimensions", dimensions )
            quantity[0].set( "Precision", str( precision[2] ) )

      tree.write( os.path.join( path_trimmed, filename_xdmf ) )

def main( arguments ):

   SanityCheck( arguments.geometric_restrictions_active, arguments.levelset_restriction_active )

   TrimFiles( arguments.directory_path, arguments.path_trimmed, arguments.cut_off_for_levelset, arguments.levelset_restriction_active, arguments.geometric_restrictions_active, \
              arguments.cylindrical_restriction, arguments.spherical_restriction, arguments.line_restriction, arguments.quantities_to_keep, arguments.precision, arguments.rearrange_vertex_coordinates )


def ParseArguments() :
   parser = ArgumentParser( description = "This script allows you to reduce the size of your output via levelset and geometrical restrictions.\n \
                                           It requires h5py and numpy. The script allows to reduce the size of the putput via different restrictions. \
                                           A mix of theses or even all restrictions can be applied simultaneously." )
   parser.add_argument( "--path", help = "Path to the .h5 and .xdmf file pair to be trimmed", default = "", dest = "directory_path" )
   parser.add_argument( "--res-path", help = "Path were the resulting files will be placed", default = "", dest = "path_trimmed" )
   parser.add_argument( "--quantities", nargs = "+", help = "Quantities to keep (e.g. levelset pressure) separated by single whitespace. Name \
                         needs to exactly match the name in the HDF5 file", default = ["levelset"], dest = "quantities_to_keep" )
   parser.add_argument( "--ls-restriction", nargs = 2, help = "The lower and upper bound of the kept levelset values", default = [], dest = "cut_off_for_levelset", type = float )
   parser.add_argument( "--line-restriction", nargs = 7, help = "Drops all data that is NOT within the given range of one or multiple axis. This way you can \
                         cut out blocks. Arguments: active axes, x0, x1, y0, y1, z0, z1", default = [], dest = "line_restriction" )
   parser.add_argument( "--cylinder-restriction", nargs = 4, help = "Drops all data that is NOT within the provided cylinder. Arguments: two orthogonal axes (as one \
                         string), the radius, coordinates of a center-line point", default = [], dest = "cylindrical_restriction" )
   parser.add_argument( "--spherical-restriction", nargs = 4, help = "drops all data that is NOT within the provided sphere. Arguments: radius, x-, y-, z-location of  center point "\
                        , default = [], dest = "spherical_restriction", type = float )
   parser.add_argument( "--precision", nargs = 3, help = "Precision for the cell vertices (int), the cell coordinates (float) and the simulation data (float). \
                         Attention - the precision must suffice to save all vertices in the trimmed output.", default = ["u8", "f8", "f8"], dest = "precision" )
   parser.add_argument( "--rearrange-vertices", help = "Drops all unused vertices, thereby reduces the file size significantly. \
                         Takes a long time. Recommended if any geometric restriction is used.", dest = "rearrange_vertex_coordinates", action = "store_true" )
   arguments = parser.parse_args()
   arguments.levelset_restriction_active = False
   arguments.geometric_restrictions_active = [False, False, False]
   if( len( arguments.cut_off_for_levelset ) != 0 ) :
      arguments.levelset_restriction_active = True
   if( len( arguments.line_restriction ) != 0 ):
      arguments.geometric_restrictions_active[0] = True
      temp = arguments.line_restriction
      temp[0] = temp[0].replace( "x", "1" )
      temp[0] = temp[0].replace( "y", "2" )
      temp[0] = temp[0].replace( "z", "3" )
      temp[1] = float( temp[1] )
      temp[2] = float( temp[2] )
      temp[3] = float( temp[3] )
      temp[4] = float( temp[4] )
      temp[5] = float( temp[5] )
      temp[6] = float( temp[6] )
      arguments.line_restriction = line_restriction = [temp[0], [(temp[1], temp[2]), (temp[3], temp[3]), (temp[5], temp[6])]]
   if( len( arguments.cylindrical_restriction ) != 0 ) :
      arguments.geometric_restrictions_active[1] = True
      temp = arguments.cylindrical_restriction
      temp[0] = temp[0].replace( "x", "1" )
      temp[0] = temp[0].replace( "y", "2" )
      temp[0] = temp[0].replace( "z", "3" )
      temp[1] = float( temp[1] )
      temp[2] = float( temp[2] )
      temp[3] = float( temp[3] )
      arguments.cylindrical_restriction = [temp[0], temp[1], ( temp[2], temp[3] )]
   if( len( arguments.spherical_restriction ) != 0 ) :
      arguments.geometric_restrictions_active[2] = True
      temp = arguments.spherical_restriction
      arguments.spherical_restriction = [temp[0], ( temp[1], temp[2], temp[3] )]

   return arguments

if __name__ == "__main__":
   arguments = ParseArguments()
   main( arguments )
