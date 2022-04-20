#! /usr/bin/env python

'''
Horizontal remapping of bathymetry dataset to model grid. 

Assume dataset higher resolution than model grid, so use
binning and gridbox averaging (median by default).   

This version uses xarray.

30/7/2021 : This version uses the method of determining if
            a source data point is in the grid cell from 
            NEMOBAT, which assumes the grid cell is a linear 
            quadrilateral in lat/lon space rather than a parallelogram.
            This is slightly more accurate and simpler than
            my original parallelogram method. DS. 

23/8/2021 : Add a redundant time dimension because cdfsmooth (and NEMO?) requires this.
            Get the x,y dimensions of the output file the "right" way round.
            Add option to fill land points with missing data indicators.
            DS.

26/8/2021 : Deal with the discontinuity in the longitudes at 180E based on the
            NEMOBAT solution for this. DS. 

27/8/2021 : Include filling (by interpolation) of "bad" cells. DS. 

@author: Dave Storkey
@date: July 2021
'''

import xarray as xr
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib
import socket
if 'spice' in socket.gethostname():
    # Note this disables plt.show()
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

def angle(v1, v2):
    mag1 = ma.sqrt(v1[0]*v1[0]+v1[1]*v1[1])
    mag2 = ma.sqrt(v2[0]*v2[0]+v2[1]*v2[1])
    # avoid rounding errors taking the argument to arccos bigger than 1.0
    angle = ma.arccos( ma.minimum( (v1[0]*v2[0] + v1[1]*v2[1])/(mag1*mag2), 1.0 ) )
    return angle * 180.0 / np.pi

def determinant(v1, v2):
    # Calculate determinant (= cross product) of two 2D vectors.
    # The individual components of v1 and v2 can be arrays. 
    determinant = v1[0]*v2[1] - v1[1]*v2[0]
    return determinant

def fill_cells(bathy_in=None, coords_list=None, nx=None, ny=None):
    # Fill in bad cells (that GEBCO thinks are land points but we want to be sea points)
    # using interpolation from nearest neighbours
    bad_cells=[]
    bathy_out=bathy_in
    for ii,jj in coords_list:
        # Exclude EW and northfold halo points from the filling
        # Rely on a halo swap happening at a later stage
        # (Note this exclusion is necessary even if we are creating a limited area model).
        if ii != nx-1 and jj != ny-1:
            bathy_out[jj,ii] = (1.0/12.0) * ( 2.0 * bathy_in[jj  ,ii-1] + \
                                              2.0 * bathy_in[jj  ,ii+1] + \
                                              2.0 * bathy_in[jj+1,ii  ] + \
                                              2.0 * bathy_in[jj-1,ii  ] + \
                                                    bathy_in[jj-1,ii-1] + \
                                                    bathy_in[jj-1,ii+1] + \
                                                    bathy_in[jj+1,ii-1] + \
                                                    bathy_in[jj+1,ii+1] )
            if bathy_out[jj,ii] <= 0.0:
                bad_cells.append([ii,jj])
    # NB. NEED TO CHECK FOR CASES THAT DON'T WORK!!
    return(bathy_out,bad_cells)

def bathy_regrid_horiz(in_bathyfile=None, meshfile=None, out_bathyfile=None, mask_in=None,
                       vischeck=None, bathy_source=None, dlon=None, dlat=None, mean=None, 
                       mask_land=None, fill_bad_cells=None, visualise_bad_cells=None, min_depth=None):

    if bathy_source is None:
        bathy_source = "gebco"
    print('bathy_source : ',bathy_source)

    if min_depth is None:
        min_depth = 0.0

    if dlon is None:
        dlon = 0.5
    if dlat is None:
        dlat = 0.5

    bathy_dataset = xr.open_dataset(in_bathyfile)
    if bathy_source == "nemo":
        lon_bathy_in = bathy_dataset.nav_lon
        lat_bathy_in = bathy_dataset.nav_lat
#        in_bathy = bathy_dataset.Bathymetry 
    elif bathy_source == "gebco":
        lon_bathy_in = bathy_dataset.lon
        lat_bathy_in = bathy_dataset.lat
#        in_bathy = bathy_dataset.elevation * -1.0
    else:
        raise Exception("Error. Only possible bathy sources: gebco or nemo.")

    if len(lon_bathy_in.shape) == 1:
        lat_bathy, lon_bathy = xr.broadcast(lat_bathy_in,lon_bathy_in)
    else:
        lat_bathy, lon_bathy = lat_bathy_in, lon_bathy_in

    print('lon_bathy.shape : ',lon_bathy.shape)
    print('lat_bathy.shape : ',lat_bathy.shape)

    with xr.open_dataset(meshfile) as grid_data:
        nx = len(grid_data.x)
        ny = len(grid_data.y)
        lonT = grid_data.glamt.squeeze('t').values
        latT = grid_data.gphit.squeeze('t').values
        lonF = grid_data.glamf.squeeze('t').values
        latF = grid_data.gphif.squeeze('t').values
        tmask = grid_data.tmask.squeeze('t').isel(z=0).values
        if not mask_in:
            # always read in tmask but set values to 1 everywhere
            # if we don't want to impose a land-sea mask. 
            tmask.values[:] = 1

    bathy_out = ma.zeros((ny,nx))
    npoints = ma.zeros((ny,nx))
    if visualise_bad_cells:
        bad_cells = ma.zeros((ny,nx))
    count_bad_cells = 0
    list_bad_cells = []
    count_coastal_cells_to_be_fixed = 0
    list_coastal_cells_to_be_fixed = []
    count_unfixed_coastal_cells = 0
    list_unfixed_coastal_cells = []
    ii_keep=-1
    ijlist_near180 = []
    for iteration in [0,1]:
        # First iteration to do all points not near the 180E line and
        # second iteration to do remaining points near the 180E line if necessary.
        if iteration == 0:
            ijlist = [[ii,jj] for ii in range(nx-1) for jj in range(ny-1)]
        elif len(ijlist_near180) > 0:
            # shift all longitudes to the range [0,360] instead of [-180,180] to deal with the strip
            # near 180E. 
            print('=== Working on the points near 180E ===')
            lon_bathy_in = xr.where(lon_bathy_in < 0.0, lon_bathy_in+360.0, lon_bathy_in)
            # rebroadcasting here is more memory efficient than shifting the existing lon_bathy array:
            lat_bathy, lon_bathy = xr.broadcast(lat_bathy_in,lon_bathy_in)
            # lon_bathy = xr.where(lon_bathy < 0.0, lon_bathy+360.0, lon_bathy)
            lonT = xr.where(lonT < 0.0, lonT+360.0, lonT)
            lonF = xr.where(lonF < 0.0, lonF+360.0, lonT)
            ijlist = ijlist_near180
        else:
            break
        for ii,jj in ijlist:
            # only process sea points
            if tmask[jj+1,ii+1]:
                if iteration == 0 and ( lonT[jj+1,ii+1] < -175.0 or lonT[jj+1,ii+1] > 175.0 ):
                    ijlist_near180.append([ii,jj])
                    continue
                if ii != ii_keep:
                    print("Working on column : ",ii)
                    ii_keep = ii
                # Note that (ii,jj) are the coordinates of the F-point at the "bottom-left" corner of the cell of
                # of interest, so the corresponding T-point where the bathymetry is valid is at (ii+1,jj+1). (NEMO
                # always has a redundant row and column of T points at ii=1 and jj=1).

                # optimisation: only check local points to see if inside grid cell
                if bathy_source == "gebco":
                    indx = np.asarray( ( lon_bathy_in > lonT[jj+1,ii+1] - dlon ) & (lon_bathy_in < lonT[jj+1,ii+1] + dlon ) ).nonzero()
                    indy = np.asarray( ( lat_bathy_in > latT[jj+1,ii+1] - dlat ) & (lat_bathy_in < latT[jj+1,ii+1] + dlat ) ).nonzero()
                    nix, niy = len(indx[0]), len(indy[0])
                    indx_bcst = ( indx[0] * np.ones((niy,nix)).astype(int) ).flatten()
                    indy_bcst = ( indy[0] * np.ones((nix,niy)).astype(int) ).transpose().flatten()
                    lon_local = lon_bathy.values[(indy_bcst,indx_bcst)]
                    lat_local = lat_bathy.values[(indy_bcst,indx_bcst)]
                    bathy_local = bathy_dataset.elevation.values[(indy_bcst,indx_bcst)] * -1.0
                else:
                    lon_local = lon_bathy.values
                    lat_local = lat_bathy.values
                    bathy_local = bathy_dataset.elevation.values

                EdgeVector1 = [lonF[jj  ,ii+1]-lonF[jj  ,ii  ],latF[jj  ,ii+1]-latF[jj  ,ii  ]]
                EdgeVector2 = [lonF[jj+1,ii+1]-lonF[jj  ,ii+1],latF[jj+1,ii+1]-latF[jj  ,ii+1]]
                EdgeVector3 = [lonF[jj+1,ii  ]-lonF[jj+1,ii+1],latF[jj+1,ii  ]-latF[jj+1,ii+1]]
                EdgeVector4 = [lonF[jj  ,ii  ]-lonF[jj+1,ii  ],latF[jj  ,ii  ]-latF[jj+1,ii  ]]

                SourceVector1 = [lon_local-lonF[jj  ,ii  ],lat_local-latF[jj  ,ii  ]]
                SourceVector2 = [lon_local-lonF[jj  ,ii+1],lat_local-latF[jj  ,ii+1]]
                SourceVector3 = [lon_local-lonF[jj+1,ii+1],lat_local-latF[jj+1,ii+1]]
                SourceVector4 = [lon_local-lonF[jj+1,ii  ],lat_local-latF[jj+1,ii  ]]

                det1 = determinant( EdgeVector1, SourceVector1 )
                det2 = determinant( EdgeVector2, SourceVector2 )
                det3 = determinant( EdgeVector3, SourceVector3 )
                det4 = determinant( EdgeVector4, SourceVector4 )

                indices = ma.where( ( det1 >= 0 ) & ( det2 >= 0 ) & ( det3 >= 0 ) & ( det4 >= 0 ) )

                if len(indices[0]) > 0:
                    # found some points in the T-cell
                    bathy_binned = bathy_local[indices]
                    npoints[jj+1,ii+1] = len(indices[0])
                    if mean:
                        bathy_out[jj+1,ii+1] = ma.average(bathy_binned[:]) 
                    else: # median
                        bathy_out[jj+1,ii+1] = ma.median(bathy_binned[:]) 
                    if mask_in and bathy_out[jj+1,ii+1] <= 0.0:
                        # in this case we are imposing a land sea mask, so any points that come 
                        # out as land at a pre-determined sea point must be "tuned" to be ocean. 
                        count_coastal_cells_to_be_fixed += 1
                        list_coastal_cells_to_be_fixed.append([ii+1,jj+1])
                        indices_sea = ma.where(bathy_binned > 0.0)
                        if len(indices_sea[0]) > 0:
                            npoints[jj+1,ii+1] = len(indices_sea[0])
                            if mean:
                                bathy_out[jj+1,ii+1] = ma.average(bathy_binned[indices_sea]) 
                            else: # median
                                bathy_out[jj+1,ii+1] = ma.median(bathy_binned[indices_sea]) 
                        else:
                            count_unfixed_coastal_cells +=1
                            list_unfixed_coastal_cells.append([ii+1,jj+1])
                    if vischeck and jj == 2 and ii == 2:
                        loncell = [lonF[jj,ii],lonF[jj,ii+1],lonF[jj+1,ii+1],lonF[jj+1,ii],lonF[jj,ii]]
                        latcell = [latF[jj,ii],latF[jj,ii+1],latF[jj+1,ii+1],latF[jj+1,ii],latF[jj,ii]]
                        lons_to_plot = lon_local[indices]
                        lats_to_plot = lat_local[indices]
                else:
                    count_bad_cells += 1
                    list_bad_cells.append([ii+1,jj+1])

    print("Number of bad cells (no source data) : ",count_bad_cells)
    if count_bad_cells > 0:
        print("Bad cells (no source data) : ",list_bad_cells[:])
        if visualise_bad_cells:
            for (ii,jj) in list_bad_cells:
                bad_cells[jj,ii] = 8.0

    if mask_in:
        
        print("Number of coastal cells to be fixed originally : ",count_coastal_cells_to_be_fixed)
        count_unfixed_new = 0
        if count_coastal_cells_to_be_fixed > 0:
            print("Coastal cells to be fixed : ",list_coastal_cells_to_be_fixed[:])
        print("Number of unfixed coastal cells : ",count_unfixed_coastal_cells)
        # Use filling to fill the no-source-data bad cells and the unfixed coastal cells
        if count_bad_cells > 0 or count_unfixed_coastal_cells > 0:
            list_unfixed_coastal_cells.extend(list_bad_cells)
            print("Unfixed coastal cells and bad cells : ",list_unfixed_coastal_cells[:])
            if fill_bad_cells:
                count_unfixed_old = count_unfixed_coastal_cells + count_bad_cells + 1
                count_unfixed_new = count_unfixed_coastal_cells + count_bad_cells
                count_while = 0
                # carry on calling fill_cells until it doesn't do anything
                while count_unfixed_new < count_unfixed_old:
                    count_while += 1
                    bathy_out, unfixed_cells = fill_cells(bathy_in=bathy_out,coords_list=list_unfixed_coastal_cells,nx=nx,ny=ny)
                    count_unfixed_old = count_unfixed_new
                    count_unfixed_new = len(unfixed_cells)
                    print("old, new counts unfixed cells : ",count_unfixed_old, count_unfixed_new)
                print("Number of filling iterations : ",count_while)
                if len(unfixed_cells) > 0:
                    print("Still have the following unfixed cells after filling: ",unfixed_cells)
        if visualise_bad_cells:
            if count_coastal_cells_to_be_fixed > 0:
                for (ii,jj) in list_coastal_cells_to_be_fixed:
                    bad_cells[jj,ii] += 1.0
            if count_unfixed_coastal_cells > 0:
                for (ii,jj) in list_unfixed_coastal_cells:
                    bad_cells[jj,ii] += 1.0
            if count_unfixed_new > 0:
                for (ii,jj) in unfixed_cells:
                    bad_cells[jj,ii] += 2.0

    if min_depth:
        if mask_in:
            bathy_out = ma.where(tmask == 1, ma.maximum(bathy_out,min_depth), 0.0)
        else:
            bathy_out = ma.where(bathy_out > 0.0, ma.maximum(bathy_out,min_depth), 0.0)
 
    if mask_land:
         if mask_in:
             bathy_out = ma.masked_where(tmask == 0, bathy_out)
         else:
             bathy_out = ma.masked_where(bathy_out <= 0.0, bathy_out)

    if visualise_bad_cells:
         if mask_in:
             bad_cells = ma.masked_where(tmask == 0, bad_cells)
         else:
             bad_cells = ma.masked_where(bathy_out <= 0.0, bad_cells)
    
    if vischeck:
        # this check assumes we are using the little test case in the Fram Strait, 
        # (imin,imax,jmin,jmax) = (1100,1104,1114,1118) for eORCA025 grid
        mercator_arctic = ccrs.Mercator(central_longitude=0.0,
                                        min_latitude=78.0, max_latitude=79.0, globe=None, latitude_true_scale=78.5)
        ax = plt.axes(projection=mercator_arctic)
        ax.set_xlim(left=80000,right=180000)
        plt.plot(lonT,latT,"bo",transform=ccrs.Geodetic())
        plt.plot(lonF,latF,"rx",transform=ccrs.Geodetic())
        plt.plot(loncell,latcell,"k-",transform=ccrs.Geodetic())
        plt.plot(lons_to_plot,lats_to_plot,"k.",transform=ccrs.Geodetic())
        plt.savefig('vischeck.png')
                      
    bathy_out_with_time = ma.zeros((1,ny,nx))
    bathy_out_with_time[0] = bathy_out
    time_counter = ma.zeros((1))
    out_dataset = xr.Dataset(
        {
            "Bathymetry" : (('time_counter','y','x'),bathy_out_with_time,{'standard_name':'sea_floor_depth','units':'m','_FillValue':-1.0e+20}),
            "npoints"    : (('y','x'),npoints,{'_FillValue':-1.0e+20}),
        },
        coords={
            "nav_lon"      : (('y','x'),lonT,{'standard_name':'longitude','units':'degrees_east','_FillValue':False}),
            "nav_lat"      : (('y','x'),latT,{'standard_name':'latitude','units':'degrees_north','_FillValue':False}),
            "time_counter" : (('time_counter'),time_counter,{'standard_name':'time','_FillValue':False}),
        },
    )

    if visualise_bad_cells:
        out_dataset["bad_cells"]=(('y','x'),bad_cells[:],{'_FillValue':-1.0e+20})

    out_dataset.to_netcdf(out_bathyfile)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-B", "--bathy", action="store",dest="in_bathyfile",
                         help="name of file with input bathymetry dataset")
    parser.add_argument("-M", "--meshfile", action="store",dest="meshfile",
                         help="name of mesh mask file (or domcfg file) for target grid")
    parser.add_argument("-m", "--mask_in", action="store_true",dest="mask_in",
                         help="Use surface tmask field in meshmask file to predetermine coastline.")
    parser.add_argument("-o", "--outfile", action="store",dest="out_bathyfile",
                         help="name of output bathymetry file")
    parser.add_argument("-V", "--vischeck", action="store_true",dest="vischeck",
                         help="switch for visual check on binning algorithm")
    parser.add_argument("-S", "--source", action="store",dest="bathy_source",
                         help="bathymetry data source, currently gebco or nemo.")
    parser.add_argument("--dlon", action="store",dest="dlon",type=float,
                         help="delta-longitude to define local search area.")
    parser.add_argument("--dlat", action="store",dest="dlat",type=float,
                         help="delta-latitude to define local search area.")
    parser.add_argument("-d", "--min_depth", action="store",dest="min_depth",type=float,
                         help="Minimum depth in metres to impose at sea points.")
    parser.add_argument("--mean", action="store_true",dest="mean",
                         help="Use mean to average binned values. (Default is median).")
    parser.add_argument("-L", "--mask_land", action="store_true",dest="mask_land",
                         help="Mask land points (where bathy <=0) with missing data indicators.")
    parser.add_argument("-F", "--fill_bad_cells", action="store_true",dest="fill_bad_cells",
                         help="Fill in bad cells with interpolation.")
    parser.add_argument("--visualise_bad_cells", action="store_true",dest="visualise_bad_cells",
                         help="Output masked field showing where bad cells are.")

    args = parser.parse_args()

    bathy_regrid_horiz(in_bathyfile=args.in_bathyfile,meshfile=args.meshfile,
                       out_bathyfile=args.out_bathyfile,mask_in=args.mask_in,
                       vischeck=args.vischeck, bathy_source=args.bathy_source,
                       dlon=args.dlon, dlat=args.dlat, mean=args.mean, 
                       mask_land=args.mask_land, fill_bad_cells=args.fill_bad_cells,
                       min_depth=args.min_depth, visualise_bad_cells=args.visualise_bad_cells )
