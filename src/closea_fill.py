#! /usr/bin/env python

'''
Read in closea masks and use to fill in inland seas in tmask field. 

@author: Dave Storkey
@date: Oct 2021
'''

import xarray as xr
import numpy as np

def closea_fill(mesh_file=None, closea_file=None, outfile=None, lakes_to_fill=None):

    with xr.open_dataset(closea_file) as closea_data:
        closea_mask = closea_data.closea_mask

    # nan_to_num converts NaNs to zeroes
    closea_indices = np.unique(np.nan_to_num(closea_mask.values)).astype(int)

    with xr.open_dataset(mesh_file) as mesh_data:
        glamt = mesh_data.glamt
        gphit = mesh_data.gphit
        glamf = mesh_data.glamf
        gphif = mesh_data.gphif
        tmask = mesh_data.tmask

    tmask_dummy,closea_mask3d = xr.broadcast(tmask,closea_mask)

    print('tmask.shape:',tmask.shape)
    print('closea_mask3d.shape:',closea_mask3d.shape)

    if lakes_to_fill is not None:
        indices_to_fill = lakes_to_fill
    else:
        # Don't fill in the zeroes of course.
        indices_to_fill = closea_indices[1:]

    for lake_index in indices_to_fill:
        print('Filling lake index : ',lake_index)
        if lake_index in closea_indices:
            print('Number of lake points : ',np.count_nonzero(closea_mask3d.astype(int) == lake_index) )
            tmask.values = np.where(closea_mask3d.astype(int) == lake_index, 0, tmask.values)

    outdata = tmask.to_dataset()    
    outdata.update({'glamt':glamt ,
                    'gphit':gphit ,
                    'glamf':glamf ,
                    'gphif':gphif  })

    outdata.to_netcdf(outfile)
                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-C", "--closea_file", action="store",dest="closea_file",
                         help="name of file with closea_mask in it")
    parser.add_argument("-M", "--mesh_file", action="store",dest="mesh_file",
                         help="name of mesh mask file")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                         help="name of output file")
    parser.add_argument("-L", "--lakes_to_fill", action="store",dest="lakes_to_fill",type=int,nargs="+",
                         help="List of integer labels for lakes to be filled. (Default fill all lakes).")

    args = parser.parse_args()

    closea_fill(closea_file=args.closea_file,mesh_file=args.mesh_file,
                lakes_to_fill=args.lakes_to_fill, outfile=args.outfile)

