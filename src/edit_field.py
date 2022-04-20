#! /usr/bin/env python

'''
Script to a 2D field at particular (i,j) points as specified
by a text file file. Generalisation of edit_bathy.py but this 
version allows you to specify ranges of i and j which will 
take the same value. 

@author: Dave Storkey
@date: Apr 2022
'''

import xarray as xr
import numpy as np
import re

def field_edit(file_in=None, file_out=None, varname=None, edits_file=None):

    with xr.open_dataset(file_in) as indata:
        try:
            nav_lat = indata.nav_lat
        except(AttributeError):
            nav_lat = None
        try:
            nav_lon = indata.nav_lon
        except(AttributeError):
            nav_lon = None
        field = getattr(indata,varname).squeeze()

    ii_list=[]
    jj_list=[]
    newvalue_list=[]
    with open(edits_file) as edits:
        for fline in edits:
            # comment lines begin with # or !
            if not re.search("^[#!]",fline):
                if "-" in fline.split(",")[0]: # specifying a range for i
                    istart=fline.split(",")[0].split("-")[0]
                    iend=fline.split(",")[0].split("-")[1]
                    ii=np.array([i for i in range(int(istart),int(iend)+1)])
                else:
                    ii=np.array([int(fline.split(",")[0])])
                if "-" in fline.split(",")[1]: # specifying a range for j
                    jstart=fline.split(",")[1].split("-")[0]
                    jend=fline.split(",")[1].split("-")[1]
                    jj=np.array([j for j in range(int(jstart),int(jend)+1)])
                else:
                    jj=np.array([int(fline.split(",")[1])])
                broadcast_multiplier=np.ones((len(ii),len(jj))).astype(int)
                ii=(ii*broadcast_multiplier.transpose()).transpose()
                jj=jj*broadcast_multiplier
                ii_list = ii_list + ii.flatten().tolist()
                jj_list = jj_list + jj.flatten().tolist()
                newvalue_list = newvalue_list + [fline.split(",")[2]]*len(ii.flatten().tolist())

    for ii,jj,newvalue in zip(ii_list,jj_list,newvalue_list):
        print("ii,jj,newvalue : ",ii,jj,newvalue)
        newvalue=float(newvalue)
        print("old value : ",field.values[jj,ii])
        field.values[jj,ii] = newvalue
        print("new value : ",field.values[jj,ii])

    outdata = field.to_dataset()    
    if nav_lat is not None and nav_lon is not None:
        outdata.update({'nav_lat':nav_lon ,
                        'nav_lat':nav_lon })

    outdata.to_netcdf(file_out)
                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--file_in", action="store",dest="file_in",
                         help="input file")
    parser.add_argument("-o", "--file_out", action="store",dest="file_out",
                         help="output file")
    parser.add_argument("-v", "--varname", action="store",dest="varname",
                         help="name of variable to edit")
    parser.add_argument("-e", "--edits_file", action="store",dest="edits_file",
                         help="file specifying editing to do")

    args = parser.parse_args()

    field_edit(file_in=args.file_in,file_out=args.file_out,varname=args.varname,
               edits_file=args.edits_file)

