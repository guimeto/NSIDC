# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:18:09 2019

@author: guillaume

"""
import os
import pandas as pd 
import numpy as np
from netCDF4 import Dataset, date2num
from datetime import datetime, timedelta, date
import calendar

cur_dir = os.path.realpath(os.path.dirname('__file__'))
    
# Get the latitude, longitude, and grid_area files.
ncols=304; nrows=448 #Number of columns and rows for the data matrices.

lonfile = open('psn25lons_v3.dat','rb')
latfile = open('psn25lats_v3.dat','rb')

lons = np.fromfile(lonfile,dtype='int32'); #lons[lons<=0]=lons[lons<=0]+360
lats = np.fromfile(latfile,dtype='int32')


# Numpy reads in the data into a 1-D array. Reshape into 2-D arrays.
lons = lons.reshape((ncols,nrows),order='F').astype('float')/1e5
lats = lats.reshape((ncols,nrows),order='F').astype('float')/1e5




for year in range(2017,2018):
    path =  str(year) + '/'
    for month in range(1,2):

        dmin, dmax =(1,calendar.monthrange(year,month)[1])        
        field = []
        for d in range(dmin,dmax+1):
            # Sample read-in for one daily file.
            icefile = open(path + 'nt_' + str(year) + '{:02d}'.format(month) + '{:02d}'.format(d) + '_f17_v1.1_n.bin','rb')  
            hdr1 = icefile.read(300)
            c = np.fromfile(icefile,dtype='uint8')
            
            # Numpy reads in the data into a 1-D array. Reshape into 2-D arrays.
            c = c.reshape((ncols,nrows),order='F').astype('float')
            icefile.close()
            
            # Values 0-250 represent fractional ice coverage scaled by 250.
            # 251: hole around pole due to orbital inclination
            # 253: coastlines
            # 254: land mask
            # 255: missing data
            c[c==251]=250; # treat hole as 100% concentration
            c[(c==253) | (c==254) | (c==255)]=np.nan; # treat coast, land, or missing as NaN
            c = c/250. # ice concentration (0 to 1)
            
            field.append(c)
        concatenated = np.dstack(field)
        
        start = date(year, month, 1)
        rng = pd.date_range(start, periods= (dmax - dmin +1), freq='D')
        tot = dmax - dmin +1
        
        ###### Écriture du fichier Netcdf en sortie
        C = Dataset('NOAA_NSIDC_SIC_'+str(year)+"_{:02d}".format(month)+'.nc', 'w', format="NETCDF4")
        C.description = 'NOAA/NSIDC Climate Data Record of Passive Microwave Daily Northern Hemisphere Sea Ice Concentration'
        C.conventions = 'CF-1.0'  
        C.model_id = 'NOAA/NSIDC'
        C.grid='latlon'
        C.institution = 'UQAM - ESCER Center, University of Quebec in Montreal'
        C.contact = 'Guillaume Dueymes'
        ########################################
        # Dimensions
        C.createDimension('x', nrows)
        C.createDimension('y', ncols)
        C.createDimension('time', tot)
        
        var=C.createVariable('seaice_conc_cdr', np.float32, ('y','x','time')) 
        var.long_name = 'Daily Northern Hemisphere Sea Ice Concentration'
        var.unit = '0-1'
        latitude=C.createVariable('latitude', np.float32, ('y','x'))
        longitude=C.createVariable('longitude', np.float32, ('y','x')) 
        latitude.units = 'degrees north'
        longitude.units = 'degrees east'
        
        
        time = C.createVariable('time', np.float64, ('time',))
        time.long_name = 'time'
        time.units = 'days since ' + str(year) + '-' + '{:02d}'.format(int(month))+'-01 00:00:00.0'
        time.calendar = "proleptic_gregorian"
        dates = [datetime(year,month,1)+n*timedelta(days=1) for n in range(concatenated.shape[2])]
        time[:] = date2num(dates,units=time.units,calendar=time.calendar)
        
        setattr(C.variables['latitude'],'Latitude','Lat')
        setattr(C.variables['longitude'],'Longitude','Lon')
         
        latitude[:,:] = lats
        longitude[:,:] = lons
        C.variables['seaice_conc_cdr'][:,:,:] = concatenated[::]
        C.close()




        
#












