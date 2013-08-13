import numpy as np
import tracmass
import netCDF4 as netCDF
import matplotlib.pyplot as plt
from datetime import datetime
import time

from traits.api import HasTraits, Date, Float, Dict, String, Array, Bool, Int

import inout
import plotting
import tools


class ModelRun(HasTraits):
    """ Represents an entire tracpy run
    """

    # Location of TXLA model output file and grid, on a thredds server.
    loc = String('http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc')

    # Number of days to run the drifters.
    ndays = Int(5)

    # Start date in date time formatting
    date = Date(datetime(2009,11, 25, 0))

    # Time between outputs
    tseas = Int(4*3600) # 4 hours between outputs, in seconds 

    # Number of interpolation steps between model outputs.
    nsteps = Int(5)

    # Use ff = 1 for forward in time and ff = -1 for backward in time.
    ff = Int(1)

    ah = Float(0)  # m^2/s
    av = Float(0)  # m^2/s

    # turbulence/diffusion flag
    doturb = Int(0)

    grid = Dict

    lon0 = Array
    lat0 = Array

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0

    ## Choose method for vertical placement of drifters
    z0 = String('s')  #'z' #'salt' #'s' (this will need to be Float for 3-D)
    zpar = Int(29)    #-10 #grid['km']-1 # 30 #grid['km']-1

    # simulation name, used for saving results into netcdf file
    name = String('temp')

    dostream = Int(0)

    def _grid_default(self):
        return inout.readgrid(loc)

    def _lon0_default(self):
        return np.linspace(-98.5, -87.5, 55)

    def _lat0_default(self):
        return np.linspace(22.5, 31, 49)

    def get_initial_points(self):
        # Eliminate points that are outside domain or in masked areas
        # decimate starting locations
        lon0, lat0 = np.meshgrid(self.lon0, self.lat0)
        self.lon0, self.lat0 = tools.check_points(lon0, lat0, self.grid)

    def _grid_default(self):
        # Read in grid parameters into dictionary, grid
        return inout.readgrid(self.loc, self.nc)

    def run(self):
        '''

        To re-compile tracmass fortran code, type "make clean" and "make f2py", which will give 
        a file tracmass.so, which is the module we import above. Then in ipython, "run run.py"
        xend,yend,zend are particle locations at next step
        some variables are not specifically because f2py is hiding them from me:
            imt, jmt, km, ntractot
        Look at tracmass.step to see what it is doing and making optional at the end.
        Do this by importing tracmass and then tracmass.step?

        I am assuming here that the velocity field at two times are being input into tracmass
        such that the output is the position for the drifters at the time corresponding to the
        second velocity time. Each drifter may take some number of steps in between, but those
        are not saved.

        loc         Path to directory of grid and output files
        nsteps      Number of steps to do between model outputs (iter in tracmass)
        ndays       number of days to track the particles from start date
        ff          ff=1 to go forward in time and ff=-1 for backward in time
        date        Start date in datetime object
        tseas       Time between outputs in seconds
        ah          Horizontal diffusion in m^2/s. 
            See project values of 350, 100, 0, 2000. For -turb,-diffusion
        av          Vertical diffusion in m^2/s.
        do3d        for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
        doturb      turbulence/diffusion flag. 
            doturb=0 means no turb/diffusion,
            doturb=1 means adding parameterized turbulence
            doturb=2 means adding diffusion on a circle
            doturb=3 means adding diffusion on an ellipse (anisodiffusion)
        lon0        Drifter starting locations in x/zonal direction.
        lat0        Drifter starting locations in y/meridional direction.
        z0/zpar     For 3D drifter movement, turn off twodim flag in makefile.
            Then z0 should be an array of initial drifter depths. 
            The array should be the same size as lon0 and be negative
            for under water. Currently drifter depths need to be above 
            the seabed for every x,y particle location for the script to run.
            To do 3D but start at surface, use z0=zeros(ia.shape) and have
                either zpar='fromMSL'
            choose fromMSL to have z0 starting depths be for that depth below the base 
            time-independent sea level (or mean sea level).
            choose 'fromZeta' to have z0 starting depths be for that depth below the
            time-dependent sea surface. Haven't quite finished the 'fromZeta' case.
            For 2D drifter movement, turn on twodim flag in makefile.
            Then: 
            set z0 to 's' for 2D along a terrain-following slice
                and zpar to be the index of s level you want to use (0 to km-1)
            set z0 to 'rho' for 2D along a density surface
                and zpar to be the density value you want to use
                Can do the same thing with salinity ('salt') or temperature ('temp')
                The model output doesn't currently have density though.
            set z0 to 'z' for 2D along a depth slice
                and zpar to be the constant (negative) depth value you want to use
            To simulate drifters at the surface, set z0 to 's' 
                and zpar = grid['km']-1 to put them in the upper s level
                z0='s' is currently not working correctly!!!
                In the meantime, do surface using the 3d set up option but with 2d flag set
        xp          x-locations in x,y coordinates for drifters
        yp          y-locations in x,y coordinates for drifters
        zp          z-locations (depths from mean sea level) for drifters
        t           time for drifter tracks
        name        Name of simulation to be used for netcdf file containing final tracks
        grid        (optional) Grid information, as read in by tracpy.inout.readgrid().

        The following inputs are for calculating Lagrangian stream functions
        dostream    Calculate streamfunctions (1) or not (0). Default is 0.
        U0, V0      (optional) Initial volume transports of drifters (m^3/s)
        U, V  (optional) Array aggregating volume transports as drifters move [imt-1,jmt], [imt,jmt-1]
        '''

        tic_start = time.time()
        tic_initial = time.time()

        # Units for time conversion with netCDF.num2date and .date2num
        units = 'seconds since 1970-01-01'

        # Number of model outputs to use
        # Adding one index so that all necessary indices are captured by this number.
        # Then the run loop uses only the indices determined by tout instead of needing
        # an extra one beyond
        tout = np.int((self.ndays * (24 * 3600)) / self.tseas + 1)

        # Convert date to number
        ncdate = netCDF.date2num(self.date, units)

        # Figure out what files will be used for this tracking
        nc, tinds = inout.setupROMSfiles(self.loc, ncdate, self.ff, tout)
        self.nc = nc
        self.tinds = tinds

        self.get_initial_points()

        # Interpolate to get starting positions in grid space
        xstart0, ystart0, _ = tools.interpolate2d(self.lon0, self.lat0, self.grid, 'd_ll2ij')
        self.xstart0, self.ystart0 = xstart0, ystart0
        # Do z a little lower down

        # Initialize seed locations 
        ia = np.ceil(xstart0) #[253]#,525]
        ja = np.ceil(ystart0) #[57]#,40]

        # don't use nan's
        ind2 = ~np.isnan(ia) * ~np.isnan(ja)
        ia = ia[ind2]
        ja = ja[ind2]
        xstart0 = xstart0[ind2]
        ystart0 = ystart0[ind2]

        dates = self.nc.variables['ocean_time'][:]   
        t0save = dates[tinds[0]] # time at start of drifter test from file in seconds since 1970-01-01, add this on at the end since it is big

        # Initialize drifter grid positions and indices
        xend = np.ones((ia.size,len(tinds)*self.nsteps))*np.nan
        yend = np.ones((ia.size,len(tinds)*self.nsteps))*np.nan
        zend = np.ones((ia.size,len(tinds)*self.nsteps))*np.nan
        zp = np.ones((ia.size,len(tinds)*self.nsteps))*np.nan
        iend = np.ones((ia.size,len(tinds)*self.nsteps))*np.nan
        jend = np.ones((ia.size,len(tinds)*self.nsteps))*np.nan
        kend = np.ones((ia.size,len(tinds)*self.nsteps))*np.nan
        ttend = np.ones((ia.size,len(tinds)*self.nsteps))*np.nan
        t = np.zeros((len(tinds)*self.nsteps+1))
        flag = np.zeros((ia.size),dtype=np.int) # initialize all exit flags for in the domain

        # Initialize vertical stuff and fluxes
        # Read initial field in - to 'new' variable since will be moved
        # at the beginning of the time loop ahead
        if isinstance(self.z0, basestring): # isoslice case
            ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[0], self.grid, nc, self.z0, self.zpar)
        else: # 3d case
            ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[0], self.grid, nc)

        ## Find zstart0 and ka
        # The k indices and z grid ratios should be on a wflux vertical grid,
        # which goes from 0 to km since the vertical velocities are defined
        # at the vertical cell edges. A drifter's grid cell is vertically bounded
        # above by the kth level and below by the (k-1)th level
        if isinstance(self.z0, basestring): # then doing a 2d isoslice
            # there is only one vertical grid cell, but with two vertically-
            # bounding edges, 0 and 1, so the initial ka value is 1 for all
            # isoslice drifters.
            ka = np.ones(ia.size) 

            # for s level isoslice, place drifters vertically at the center 
            # of the grid cell since that is where the u/v flux info is from.
            # For a rho/temp/density isoslice, we treat it the same way, such
            # that the u/v flux info taken at a specific rho/temp/density value
            # is treated as being at the center of the grid cells vertically.
            zstart0 = np.ones(ia.size)*0.5

        else:   # 3d case
            raise NotImplementedError

        # Find initial cell depths to concatenate to beginning of drifter tracks later
        zsave = tools.interpolate3d(xstart0, ystart0, zstart0, zwtnew)

        toc_initial = time.time()

        # break function here
        tinds = self.tinds
        nsteps = self.nsteps
        # j = 0 # index for number of saved steps for drifters
        tic_read = np.zeros(len(tinds))
        toc_read = np.zeros(len(tinds))
        tic_zinterp = np.zeros(len(tinds))
        toc_zinterp = np.zeros(len(tinds))
        tic_tracmass = np.zeros(len(tinds))
        toc_tracmass = np.zeros(len(tinds))
        # pdb.set_trace()
        xr3 = self.grid['xr'].reshape((self.grid['xr'].shape[0], self.grid['xr'].shape[1],1)).repeat(zwtnew.shape[2],axis=2)
        yr3 = self.grid['yr'].reshape((self.grid['yr'].shape[0], self.grid['yr'].shape[1],1)).repeat(zwtnew.shape[2],axis=2)
        # Loop through model outputs. tinds is in proper order for moving forward
        # or backward in time, I think.
        for j, tind in enumerate(tinds[:-1]):
            # pdb.set_trace()
            # Move previous new time step to old time step info
            ufold = ufnew
            vfold = vfnew
            dztold = dztnew
            zrtold = zrtnew
            zwtold = zwtnew

            tic_read[j] = time.time()
            # Read stuff in for next time loop
            if isinstance(self.z0, basestring): # isoslice case
                ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[j+1], self.grid, self.nc, self.z0, self.zpar)
            else: # 3d case
                ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[j+1], self.grid, self.nc)
            toc_read[j] = time.time()
            # print "readfields run time:",toc_read-tic_read

            print j
            #  flux fields at starting time for this step
            if j != 0:
                xstart = xend[:,j*self.nsteps-1]
                ystart = yend[:,j*self.nsteps-1]
                zstart = zend[:,j*self.nsteps-1]
                # mask out drifters that have exited the domain
                xstart = np.ma.masked_where(flag[:]==1,xstart)
                ystart = np.ma.masked_where(flag[:]==1,ystart)
                zstart = np.ma.masked_where(flag[:]==1,zstart)
                ind = (flag[:] == 0) # indices where the drifters are still inside the domain
            else: # first loop, j==0
                xstart = xstart0
                ystart = ystart0
                zstart = zstart0
                # TODO: Do a check to make sure all drifter starting locations are within domain
                ind = (flag[:] == 0) # indices where the drifters are inside the domain to start

            # Find drifter locations
            # only send unmasked values to step
            if not np.ma.compressed(xstart).any(): # exit if all of the drifters have exited the domain
                break
            else:
                # Combine times for arrays for input to tracmass
                # from [ixjxk] to [ixjxkxt]
                # Change ordering for these three arrays here instead of in readfields since
                # concatenate does not seem to preserve ordering
                uflux = np.asfortranarray(np.concatenate((ufold.reshape(np.append(ufold.shape,1)), \
                                        ufnew.reshape(np.append(ufnew.shape,1))), \
                                        axis=ufold.ndim))
                vflux = np.asfortranarray(np.concatenate((vfold.reshape(np.append(vfold.shape,1)), \
                                        vfnew.reshape(np.append(vfnew.shape,1))), \
                                        axis=vfold.ndim))
                dzt = np.asfortranarray(np.concatenate((dztold.reshape(np.append(dztold.shape,1)), \
                                        dztnew.reshape(np.append(dztnew.shape,1))), \
                                        axis=dztold.ndim))

                # Change the horizontal indices from python to fortran indexing 
                # (vertical are zero-based in tracmass)
                xstart, ystart = tools.convert_indices('py2f',xstart,ystart)

                # km that is sent to tracmass is determined from uflux (see tracmass?)
                # so it will be the correct value for whether we are doing the 3D
                # or isoslice case.
                # vec = np.arange(j*nsteps,j*nsteps+nsteps) # indices for storing new track locations
                tic_tracmass[j] = time.time()
                # pdb.set_trace()
                if self.dostream: # calculate Lagrangian stream functions
                    xend[ind,j*nsteps:j*nsteps+nsteps],\
                        yend[ind,j*nsteps:j*nsteps+nsteps],\
                        zend[ind,j*nsteps:j*nsteps+nsteps], \
                        iend[ind,j*nsteps:j*nsteps+nsteps],\
                        jend[ind,j*nsteps:j*nsteps+nsteps],\
                        kend[ind,j*nsteps:j*nsteps+nsteps],\
                        flag[ind],\
                        ttend[ind,j*nsteps:j*nsteps+nsteps], U, V = \
                            tracmass.step(np.ma.compressed(xstart),\
                                            np.ma.compressed(ystart),
                                            np.ma.compressed(zstart),
                                            self.tseas, uflux, vflux, self.ff, \
                                            self.grid['kmt'].astype(int), \
                                            dzt, self.grid['dxdy'], self.grid['dxv'], \
                                            self.grid['dyu'], self.grid['h'], nsteps, \
                                            self.ah, self.av, self.do3d, self.doturb, self.dostream, \
                                            t0=self.T0[ind],
                                            ut=self.U, vt=self.V)
                else: # don't calculate Lagrangian stream functions
                    xend[ind,j*nsteps:j*nsteps+nsteps],\
                        yend[ind,j*nsteps:j*nsteps+nsteps],\
                        zend[ind,j*nsteps:j*nsteps+nsteps], \
                        iend[ind,j*nsteps:j*nsteps+nsteps],\
                        jend[ind,j*nsteps:j*nsteps+nsteps],\
                        kend[ind,j*nsteps:j*nsteps+nsteps],\
                        flag[ind],\
                        ttend[ind,j*nsteps:j*nsteps+nsteps], _, _ = \
                            tracmass.step(np.ma.compressed(xstart),\
                                            np.ma.compressed(ystart),
                                            np.ma.compressed(zstart),
                                            self.tseas, uflux, vflux, self.ff, \
                                            self.grid['kmt'].astype(int), \
                                            dzt, self.grid['dxdy'], self.grid['dxv'], \
                                            self.grid['dyu'], self.grid['h'], nsteps, \
                                            self.ah, self.av, self.do3d, self.doturb, self.dostream)
                toc_tracmass[j] = time.time()
                # pdb.set_trace()

                # Change the horizontal indices from python to fortran indexing
                xend[ind,j*nsteps:j*nsteps+nsteps], \
                    yend[ind,j*nsteps:j*nsteps+nsteps] \
                                    = tools.convert_indices('f2py', \
                                        xend[ind,j*nsteps:j*nsteps+nsteps], \
                                        yend[ind,j*nsteps:j*nsteps+nsteps])

                # Calculate times for the output frequency
                if self.ff == 1:
                    t[j*nsteps+1:j*nsteps+nsteps+1] = t[j*nsteps] + np.linspace(self.tseas/nsteps,self.tseas,nsteps) # update time in seconds to match drifters
                else:
                    t[j*nsteps+1:j*nsteps+nsteps+1] = t[j*nsteps] - np.linspace(self.tseas/nsteps,self.tseas,nsteps) # update time in seconds to match drifters

                # Skip calculating real z position if we are doing surface-only drifters anyway
                if self.z0 != 's' and self.zpar != grid['km']-1:
                    tic_zinterp[j] = time.time()
                    # Calculate real z position
                    r = np.linspace(1./nsteps,1,nsteps) # linear time interpolation constant that is used in tracmass

                    for n in xrange(nsteps): # loop through time steps
                        # interpolate to a specific output time
                        # pdb.set_trace()
                        zwt = (1.-r[n])*zwtold + r[n]*zwtnew
                        zp[ind,j*nsteps:j*nsteps+nsteps], dt = tools.interpolate3d(xend[ind,j*nsteps:j*nsteps+nsteps], \
                                                                yend[ind,j*nsteps:j*nsteps+nsteps], \
                                                                zend[ind,j*nsteps:j*nsteps+nsteps], \
                                                                zwt)
                    toc_zinterp[j] = time.time()

        self.nc.close()
        t = t + t0save # add back in base time in seconds

        # pdb.set_trace()

        # Add on to front location for first time step
        xg=np.concatenate((xstart0.reshape(xstart0.size,1),xend),axis=1)
        yg=np.concatenate((ystart0.reshape(ystart0.size,1),yend),axis=1)
        # Concatenate zp with initial real space positions
        zp=np.concatenate((zsave[0].reshape(zstart0.size,1),zp),axis=1)

        # Delaunay interpolation
        # xp, yp, dt = tools.interpolate(xg,yg,grid,'d_ij2xy')
        # lonp, latp, dt = tools.interpolation(xg,yg,grid,'d_ij2ll')

        ## map coordinates interpolation
        # xp2, yp2, dt = tools.interpolate(xg,yg,grid,'m_ij2xy')
        # tic = time.time()
        lonp, latp, dt = tools.interpolate2d(xg,yg,self.grid,'m_ij2ll',mode='constant',cval=np.nan)
        # print '2d interp time=', time.time()-tic

        # pdb.set_trace()

        runtime = time.time()-tic_start

        print "run time:\t\t\t", runtime
        print "---------------------------------------------"
        print "Time spent on:"

        initialtime = toc_initial-tic_initial
        print "\tInitial stuff: \t\t%4.2f (%4.2f%%)" % (initialtime, (initialtime/runtime)*100)

        readtime = np.sum(toc_read-tic_read)
        print "\tReading in fields: \t%4.2f (%4.2f%%)" % (readtime, (readtime/runtime)*100)

        zinterptime = np.sum(toc_zinterp-tic_zinterp)
        print "\tZ interpolation: \t%4.2f (%4.2f%%)" % (zinterptime, (zinterptime/runtime)*100)

        tractime = np.sum(toc_tracmass-tic_tracmass)
        print "\tTracmass: \t\t%4.2f (%4.2f%%)" % (tractime, (tractime/runtime)*100)

        self.lonp, self.latp = lonp, latp
        self.zp = zp
        self.t = t

    def save_to_netcdf(self):
        """ Save results to netcdf file
        Not yet implemented in Traits model
        """
        raise NotImplementedError

        if dostream:
            inout.savetracks(lonp, latp, zp, t, name, nsteps, ff, tseas, ah, av, \
                                do3d, doturb, loc, T0, U, V)
            return lonp, latp, zp, t, grid, T0, U, V
        else:
            inout.savetracks(lonp, latp, zp, t, name, nsteps, ff, tseas, ah, av, \
                                do3d, doturb, loc)
            return lonp, latp, zp, t, grid

    def plot_tracks(self):
        lonp, latp, name, grid = self.lonp, self.latp, self.name, self.grid
        # Plot tracks
        plotting.tracks(lonp, latp, name, grid=grid)

        # Plot final location (by time index) histogram
        plotting.hist(lonp, latp, name, grid=grid, which='contour')
        plotting.hist(lonp, latp, name, grid=grid, which='pcolor')
        plt.show()

if __name__ == '__main__':
    model = ModelRun()
    model.run()
    model.plot_tracks()