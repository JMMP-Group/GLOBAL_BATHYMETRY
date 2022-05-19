#! /usr/bin/env python

'''
Read in closea masks and use them to copy the bathymetry 
for inland seas from one bathymetry file to another. 

@author: Dave Storkey
@date: Nov 2021
'''

import xarray as xr
import numpy as np

def closea_copy_bathy(bathy_source=None, bathy_target=None, bathy_out=None, 
                      closea_file=None, lakes_to_copy=None):

    with xr.open_dataset(closea_file) as closea_data:
        closea_mask = closea_data.closea_mask

    # nan_to_num converts NaNs to zeroes
    closea_indices = np.unique(np.nan_to_num(closea_mask.values)).astype(int)

    with xr.open_dataset(bathy_source) as source_data:
        bathy_source = source_data.Bathymetry

    with xr.open_dataset(bathy_target) as target_data:
        bathy_target = target_data.Bathymetry

    if lakes_to_copy is not None:
        indices_to_copy = lakes_to_copy
    else:
        # Don't fill in the zeroes of course.
        indices_to_copy = closea_indices[1:]

    for lake_index in indices_to_copy:
        print('Copying lake index : ',lake_index)
        if lake_index in closea_indices:
            print('Number of lake points : ',np.count_nonzero(closea_mask.astype(int) == lake_index) )
            bathy_target.values = np.where(closea_mask.astype(int) == lake_index, bathy_source.values, bathy_target.values)

    outdata = bathy_target.to_dataset()    
    outdata.to_netcdf(bathy_out)
                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-C", "--closea_file", action="store",dest="closea_file",
                         help="name of file with closea_mask in it")
    parser.add_argument("-S", "--bathy_source", action="store",dest="bathy_source",
                         help="source bathymetry to copy from")
    parser.add_argument("-T", "--bathy_target", action="store",dest="bathy_target",
                         help="source bathymetry to copy to")
    parser.add_argument("-o", "--outfile", action="store",dest="bathy_out",
                         help="name of output file")
    parser.add_argument("-L", "--lakes_to_copy", action="store",dest="lakes_to_copy",type=int,nargs="+",
                         help="List of integer labels for lakes to be copied. (Default copy all lakes).")

    args = parser.parse_args()

    closea_copy_bathy(closea_file=args.closea_file,bathy_source=args.bathy_source,bathy_target=args.bathy_target,
                      lakes_to_copy=args.lakes_to_copy, bathy_out=args.bathy_out)

