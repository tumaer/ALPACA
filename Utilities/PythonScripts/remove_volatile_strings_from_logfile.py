#!/usr/bin/env python3

import os.path
import subprocess as sp
from argparse import ArgumentParser

def main( output_file_full_path ) :
   sp.run( ["sed", "-i", "/Load\ Balancing\ *(/d", str( output_file_full_path )] )
   sp.run( ["sed", "-i", "/Number\ of\ MPI\ ranks/d", str( output_file_full_path )] )
   sp.run( ["sed", "-i", "/Output\ Folder/d", str( output_file_full_path )] )
   sp.run( ["sed", "-i", "/Total\ Time\ Spent/d", str( output_file_full_path )] )
   sp.run( ["sed", "-i", "/Restart\ File/d", str( output_file_full_path )] )
   sp.run( ["sed", "-i", "/Wall\ clock/d", str( output_file_full_path )] )

def RemoveVolatileStringFromLogFile( file_name ) :
   main( os.path.abspath( os.path.realpath( file_name ) ) )

def ParseArguments() :
   parser = ArgumentParser( description = "Removes all lines which (may) change between runs of the same test case like runtime, number of ranks etc." )
   parser.add_argument( "file", help = "The file to be altered" )
   arguments = parser.parse_args()
   output_file_path = os.path.abspath( os.path.realpath( arguments.file ) )
   return output_file_path

if __name__ == "__main__":
   arguments = ParseArguments()
   main( arguments )
