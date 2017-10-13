# -------------------------------------------------------
# author		: Johannes Horak
# contact		: johannes.horak@uibk.ac.at
# -------------------------------------------------------
# this script generates two files:
# * high resolution idealized topography and
# * high resolution forcing file with idealized conditions.
# the specifics are supplied to the script with
# the corresponding parameters.
# VERY BASIC SCRIPT SO FAR!

import getopt
import xarray as xa
import pandas as pd
import sys
import glob as glob
import numpy as np
from math import radians, cos, sin, asin, sqrt
from datetime import datetime,timedelta

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2.0)**2.0 + cos(lat1) * cos(lat2) * sin(dlon/2.0)**2.0
    c = 2.0 * asin(sqrt(a)) 
    r = 6371.0 # Radius of earth in kilometers. Use 3956 for miles
    return c * r
    
def tpot_from_N(N,z):
	Tb=288.15
	g0=9.80665
	tpot = Tb*np.exp((N**2/g0)*z)
	return tpot

def tpot_from_t_and_p(T,p):
	theta = T*(10.0**5.0/p)**(2.0/7.0)
	return theta

def t_of_z(h):
	if 0<= h and h < 11000:
		Tb = 288.15
		Lb = -0.0065
		hb = 0
		b=0
	elif 11000 <= h and h < 20000:
		Tb = 216.65
		Lb = 0.0
		hb = 11000
		b=1
	elif 20000 <= h and h < 32000:
		Tb = 216.65
		Lb = 0.001
		hb = 20000
		b=2
	elif 32000 <= h and h < 47000:
		Tb = 228.65
		Lb = 0.0028
		hb = 32000
		b=3
	elif 47000 <= h and h < 51000:
		Tb = 270.65
		Lb = 0.0
		hb = 47000
		b=4
	elif 51000 <= h and h < 71000:
		Tb = 270.65
		Lb = -0.0028
		hb = 51000
		b=5
	else:
		Tb = 214.65
		Lb = -0.002
		hb = 70000
		b=6
	
	T = Tb+(h-hb)*Lb
	return T

def barometric_formula(h):
	M=0.0289644		# molar mass of air, kg/mol
	g0=9.80665		# gravitational acceleration m/s^2
	R=8.3144598		# universal gas constant, J/(mol*K)

	if 0<= h and h < 11000:
		Tb = 288.15
		Lb = -0.0065
		Pb = 101325.0
		hb = 0.0
		b=0
	elif 11000 <= h and h < 20000:
		Tb = 216.65
		Lb = 0.000001
		Pb = 22632.10
		hb = 11000.0
		b=1
	elif 20000 <= h and h < 32000:
		Tb = 216.65
		Lb = 0.001
		Pb = 5474.89
		hb = 20000.0
		b=2
	elif 32000 <= h and h < 47000:
		Tb = 228.65
		Lb = 0.0028
		Pb = 868.02
		hb = 32000.0
		b=3
	elif 47000 <= h and h < 51000:
		Tb = 270.65
		Lb = 0.000001
		Pb = 110.91
		hb = 47000.0
		b=4
	elif 51000 <= h and h < 71000:
		Tb = 270.65
		Lb = -0.0028
		Pb = 66.94
		hb = 51000.0
		b=5
	else:
		Tb = 214.65
		Lb = -0.002
		Pb = 3.96
		hb = 70000.0
		b=6

	P = Pb*(Tb/(Tb+Lb*(h-hb)))**(g0*M/(R*Lb))		# as seen in, e.g. https://en.wikipedia.org/wiki/Barometric_formula
	return P

def mwrite(string):
	sys.stdout.write(string)
	sys.stdout.flush()
	

dlon=None
dlat=None
dx=None
dy=None
Lx=None
Ly=None
Llon=None
Llat=None
a0=None
a1=None
topo=None
ws=None
ws_angle=None
ws_profile=None
Nz=None
ztop=None
Nbv=None
Nbvconst=False

# READ COMMAND LINE OPTIONS
try:
	opts, args = getopt.getopt(sys.argv[1:],"",["dlon=","dlat=","dx=","dy=","Lx=","Ly=","Llon=","Llat=","a0=","a1=","topo=","ws=","ws_angle=","Nz=","ztop=","Nbv="])
except getopt.GetoptError:
	print 'error occured, syntax:'
	print 'evaluator.py --mode=<mode> --scnconfig=<scenario_config> --reafile=<forcing file> --hrtopo=<high resolution topography>'
	sys.exit(2)
for opt, arg in opts:
	#print opt," ",arg
	if opt in ("-h","--help"):
		print 'evaluator.py --mode=<mode> --scnconfig=<scenario_config> --reafile=<forcing file> --hrtopo=<high resolution topography>'
		sys.exit()
	elif opt in ("--dlon"):
		dlon = float(arg)/60.0
	elif opt in ("--dlat"):
		dlat = float(arg)/60.0
	elif opt in ("--dx"):
		dx   = float(arg)
	elif opt in ("--dy"):
		dy   = float(arg)
	elif opt in ("--Lx"):
		Lx   = float(arg)
	elif opt in ("--Ly"):
		Ly   = float(arg)
	elif opt in ("--Llon"):
		Llon = float(arg)
	elif opt in ("--Llat"):
		Llat = float(arg)
	elif opt in ("--a0"):
		a0 = float(arg)
	elif opt in ("--a1"):
		a1 = float(arg)
	elif opt in ("--topo"):
		topo=arg
	elif opt in ("--ws"):
		ws=float(arg)
	elif opt in ("---ws_angle"):
		ws_angle=float(arg)
	elif opt in ("--Nz"):
		Nz = int(arg)
	elif opt in ("--ztop"):
		ztop = float(arg)
	elif opt in ("--Nbv"):
		Nbv = float(arg)
		Nbvconst = True

# dlon		... longitudinal resolution in degrees
# dlat		... latitudinal resolution in degrees
# dx		... longitudinal resolution in km
# dy		... latitudinal resolution in km
# Lx		... longitudinal width of domain in km
# Ly		... latitudinal width of domain in km
# Llon  	... longitudinal width of domain in degree
# Llat  	... latitudinal width of domain in degree
# ws		... windspeed in m/s
# ws_angle	... direction of flow in degree (0...north, 90deg ...east, 180deg ...south, 270deg ...west)
# ws_profile... z-profile to apply to windspeed. Standard: const (no change with z)
# Nz		... how many z-levels to use in forcing
# ztop		... how far above the topography the top level is to be placed in km
# Nbv		... Brunt Vaisala frequency, supply a constant Brunt Vaisala frequency and calculated potential temperature from there in s**-1

# topo	... which topography to use - standard: witch of agnesi
# a0	... parameter 0 for idealized topography
# a1	... parameter 1 for idealized topography


lonc = 0.0	# longitude around which domain is centered
latc = 0.0	# latitude around which domain is centered

if topo is None:	# if topo parameter is not supplied we set to a standard here
	topo='witch'
	
if topo == 'witch':
	if a0 is None or a1 is None:
		print "  witch of agnesi requires two characteristic parameters:"
		print ""
		print "  equation:"
		print "    h(x) = a0*(a1**2/(a1**2+x**2)"
		print ""
		print "  a0 ... height"
		print "  a1 ... width"
		print ""
		sys.exit(1)


print "* Checking supplied parameters"
print "-----------------------------------------------------------------"
if ws_profile is None:
	print "    {:20s}: not supplied, setting windspeed constant with height".format("ws-z-profile")
	ws_profile='const'
	
if ws is None:
	print "    {:20s}: no backgrond windspeed supplied, using ws=10m/s".format("ws")
	ws = 10.0
	
if ws_angle is None:
	print "    {:20s}: no angle supplied, using ws_angle=270.0 degree (westwind)".format("ws_angle")
	ws_angle = 270.0+90.0		# add 90 degrees so that sin and cos yield values that match the orientation (N..0,E..90,S..180,W..270)
else:
	ws_angle+=90.0				# add 90 degrees so that sin and cos yield values that match the orientation (N..0,E..90,S..180,W..270)

if Nz is None:
	print "    {:20s}: number of z-levels not supplied, using Nz=30".format("Nz")
	Nz = 30

if ztop is None:
	print "    {:20s}: top level height above topography not supplied, using ztop=40km".format("ztop")
	ztop = 40.0
	
if Nbvconst == True:
	print "    {:20s}: calculating potential temperature so that N =={:8.4f} is constant".format("N",Nbv)
else:
	print "    {:20s}: calculating potential temperature from T and p".format("N")	

if not(dlon is None) and not(dlat is None) and not(Llon is None) and not(Llat is None):
	dx = haversine(lonc-dlon/2.0,latc,lonc+dlon/2.0,latc)
	dy = haversine(lonc,latc-dlat/2.0,lonc,latc+dlat/2.0)
	Lx = haversine(lonc-Llon/2.0,latc,lonc+Llon/2.0,latc)
	Ly = haversine(lonc,latc-Llat/2.0,lonc,latc+Llat/2.0)
elif not(dx is None) and not(dy is None) and not (Lx is None) and not (Ly is None):
	dg_lon = haversine(lonc-0.5,latc,lonc+0.5,latc)			# yields km per degree of longitude (km/deg)
	dg_lat = haversine(lonc,latc-0.5,lonc,latc+0.5)			# yields km per degree of latitude  (km/deg)
	# find Llon and Llat
	

	# find longitudinal extension by calling haversine until
	# we find the angle that best approximates the distance
	test_min=99999.0
	angle = None
	for test_llon in np.arange(0.0,90.0,0.001):
		test_result1 = haversine(lonc-test_llon/2.0,latc,lonc+test_llon/2.0,latc)
		test_result2 = np.abs(test_result1 - Lx)
		if  test_result2 < test_min:
			test_min = test_result2
			angle = test_llon
	Llon=angle
	# find latitudinal extension by calling haversine until
	# we find the angle that best approximates the distance
	test_min=99999.0
	angle = None
	for test_llat in np.arange(0.0,90.0,0.001):
		test_result1 = haversine(lonc,latc-test_llat/2.0,lonc,latc+test_llat/2.0)
		test_result2 = np.abs(test_result1 - Ly)
		if  test_result2 < test_min:
			test_min = test_result2
			angle = test_llat
	Llat = angle

	dlon = 1.0/dg_lon
	dlat = 1.0/dg_lat

else:
	print " error: either (dlon,dlat) and (Llon,Llat) or (dx,dy) and (Lx,Ly) must be specified!"
	sys.exit(1)


dz	   = ztop*1000.0/float(Nz)							# thickness of layers in m
dg_lon = haversine(lonc-0.5,latc,lonc+0.5,latc)			# km per degree of longitude
dg_lat = haversine(lonc,latc-0.5,lonc,latc+0.5)			# km per degree of latitude

Nlon = int(np.ceil(Llon/dlon))							# grid cells along longitudinal axis
if Nlon % 2 == 0:										# go for uneven number of grid cells so that
	Nlon+=1												# N=0 is at lonc
Nlat = int(np.ceil(Llat/dlat))							# grid cells along latitudinal axis
if Nlat % 2 == 0:										# see above
	Nlat+=1

print "* Generating idealized topography and forcing for ICAR experiment"
print "-----------------------------------------------------------------"
print "    centering longitude: {:2.1f}".format(lonc)
print "    centering latitude : {:2.1f}".format(latc)
print "    1 deg. longitude   : {:2.0f}km".format(dg_lon)
print "    1 deg. latitude    : {:2.0f}km".format(dg_lat)
print "    *"
print "    dlon,dlat          : {:6.2f} '      {:6.2f}'".format(dlon*60.0,dlat*60.0)
print "    dx,dy              : {:6.2f} km     {:6.2f} km".format(dx,dy)
print "    domain extension   : {:6.2f} deg.   {:6.2f} deg.".format(Llon,Llat)
print "    domain extension   : {:6.2f} km     {:6.2f} km".format(Lx,Ly)
print "    Nz                 : {:n}".format(Nz)
print "    ztop               : {:6.2f} km".format(ztop)
print "    dz                 : {:6.2f} m".format(dz)
print "    *"
print "    grid cells         : {:6n}        {:6n}".format(Nlon,Nlat)
# generate the grid

Nlon_l = -(Nlon-1.0)/2.0				# lowest longitudinal N 
Nlon_h = -Nlon_l						# highest longitudinal N
Nlat_l = -(Nlat-1.0)/2.0				# lowest latitudinal N
Nlat_h = -Nlat_l						# highest latitudinal N

lonN=np.arange(Nlon_l,Nlon_h+1.0,1.0)	# array that contains number of each longitudinal cell
latN=np.arange(Nlat_l,Nlat_h+1.0,1.0)	# array that contains number of each latitudinal cell

lon=lonN*dlon							# array that contains all discrete longitudes
lat=latN*dlat							# array that contains all discrete latitudes

lon_gridded, lat_gridded = np.meshgrid(lon, lat)	# grids of longitude and latitude

Ndays		= 1 						# number of days for which forcing should be generated
Ntime		= Ndays*24

# GENERATE ICAR HIGH RESOLUTION TOPGRAPHY
i_topo  	= np.empty(Nlat*Nlon).reshape(Nlat,Nlon)	# array that will contain topography elevation in meters
i_landmask 	= np.ones(Nlat*Nlon).reshape(Nlat,Nlon)	# array that contains the landmask. if topography > 0 => 1, otherwise 0.
i_xlong_m 	= lon_gridded #np.empty(Nlat*Nlon).reshape(Nlat,Nlon)
i_xlat_m	= lat_gridded #np.empty(Nlat*Nlon).reshape(Nlat,Nlon)



for nlat,ny in enumerate(latN):
	y=ny*dy*1000.0
	for nlon,nx in enumerate(lonN):
		x=nx*dx*1000.0
		
		h=a0*(a1**2/(x**2+a1**2))				# topography
		
		i_topo[nlat,nlon]=h
		if h <= 0:
			i_landmask[nlat,nlon]=0				# landmask

icar_topo_ds	= xa.Dataset(
							data_vars={
								'HGT_M':(['south_north','west_east'],i_topo),
								'XLONG_M':(['south_north','west_east'],i_xlong_m),
								'XLAT_M':(['south_north','west_east'],i_xlat_m),
								'LANDMASK':(['south_north','west_east'],i_landmask)
							},
							coords={
								'south_north':range(0,Nlat),
								'west_east':range(0,Nlon)
							}
						)

icar_topo_ds.to_netcdf("./ideal_output/ideal_topo.nc",format='NETCDF4')	

dtime_base  	  = datetime(1900,1,1,0,0,0)
dtime_start 	  = datetime(2017,1,1,0,0,0)
dtime_end   	  = dtime_start + timedelta(days=Ndays)
dtime_diff_start  = (dtime_start-dtime_base)
dtime_diff_end	  = (dtime_end-dtime_base)
hours_since_start = dtime_diff_start.days*24.0	+ dtime_diff_start.seconds/3600.0
hours_since_end	  = dtime_diff_end.days*24.0	+ dtime_diff_end.seconds/3600.0
time_list	= np.arange(hours_since_start,hours_since_end,1.0)

i_Time		= time_list
i_hgt 		= np.empty(Ntime*Nlat*Nlon).reshape(Ntime,Nlat,Nlon)
i_xlong		= np.empty(Ntime*Nlat*Nlon).reshape(Ntime,Nlat,Nlon)
i_xlat		= np.empty(Ntime*Nlat*Nlon).reshape(Ntime,Nlat,Nlon)
i_u			= np.empty(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# u component of background flow
i_v			= np.empty(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# v component of background flow
i_p			= np.empty(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# pressured at every cell
i_tpot		= np.empty(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# potential temperatur at every cell
i_sp		= np.empty(Ntime*Nlat*Nlon).reshape(Ntime,Nlat,Nlon)		# surface pressure
i_ph		= np.empty(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# elevation of each cell
# unused variables
#i_phb		= np.zeros(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# height base (not used, zero everywhere)
i_phb		= np.zeros(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# phb (not used, zero everywhere)
i_pb		= np.zeros(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# pb (not used, zero everywhere)
i_qvapor	= np.zeros(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# liquid water content of cell (not used, zero everywhere)
i_qcloud	= np.zeros(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# liquid water content of cell (not used, zero everywhere)
i_qice		= np.zeros(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# ice content of cell (not used, zero everywhere)
i_tsk		= np.zeros(Ntime*Nz*Nlat*Nlon).reshape(Ntime,Nz,Nlat,Nlon)		# tsk (not used, zero everywhere)


print "    generating idealized forcing (this may take a while)..."
print "    *"

p0 			= 101325.0 		# pressure at h=0m, Pa

i_xlong[:]	= lon_gridded
i_xlat[:]	= lat_gridded

# GENERATE ICAR HIGH RESOLUTION FORCING
for ntime in range(0,1):
	mwrite("\r    working on timestep {:3n}/{:3n}".format(ntime+1,Ntime))
	for nlat,ny in enumerate(latN):
		y=ny*dy*1000.0
		for nlon,nx in enumerate(lonN):
			x=nx*dx*1000.0				
			h=a0*(a1**2/(x**2+a1**2))				# topography
			
			i_hgt[ntime,nlat,nlon] = h
			i_sp[ntime,nlat,nlon]  = barometric_formula(h)
			
			for nz in range(0,Nz):
				zp=nz*dz+dz*0.5							# height of cell center above topograhy
				z=zp+h									# absolute height
				
				u=ws*np.cos(ws_angle*np.pi/180.0)		# U wind component
				v=ws*np.sin(ws_angle*np.pi/180.0)		# V wind component
				p=barometric_formula(z)
				
				if Nbvconst == True:
					tpot = tpot_from_N(Nbv,z)
				else:
					tpot = tpot_from_t_and_p(t_of_z(z),p)
					
				i_u[ntime,nz,nlat,nlon] = u
				i_v[ntime,nz,nlat,nlon] = v
				i_p[ntime,nz,nlat,nlon] = p
				i_ph[ntime,nz,nlat,nlon] = z
				i_tpot[ntime,nz,nlat,nlon] = tpot

for ntime in range(1,Ntime):
	mwrite("\r    working on timestep {:3n}/{:3n}".format(ntime+1,Ntime))
	i_hgt[ntime,:] = i_hgt[0,:]#.copy()
	i_u[ntime,:] = i_u[0,:]#.copy()
	i_v[ntime,:] = i_v[0,:]#.copy()
	i_p[ntime,:] = i_p[0,:]#.copy()
	i_tpot[ntime,:] = i_tpot[0,:]#.copy()
	i_ph[ntime,:] = i_ph[0,:]#.copy()
	i_sp[ntime,:] = i_sp[0,:]#.copy()



print ""
print "    saving..."
icar_topo_ds	= xa.Dataset(
							data_vars={
								'HGT':(['Time','south_north','west_east'],i_hgt),
								'XLONG':(['Time','south_north','west_east'],i_xlong),
								'XLAT':(['Time','south_north','west_east'],i_xlat),
								'P':(['Time','bottom_top','south_north','west_east'],i_p),
								'T':(['Time','bottom_top','south_north','west_east'],i_tpot),
								'PH':(['Time','bottom_top','south_north','west_east'],i_ph),
								'U':(['Time','bottom_top','south_north','west_east'],i_u),
								'V':(['Time','bottom_top','south_north','west_east'],i_v),
								'sp':(['Time','south_north','west_east'],i_sp),
								'PB':(['Time','bottom_top','south_north','west_east'],i_pb),
								'PHB':(['Time','bottom_top','south_north','west_east'],i_phb),
								'QVAPOR':(['Time','bottom_top','south_north','west_east'],i_qvapor),
								'QCLOUD':(['Time','bottom_top','south_north','west_east'],i_qcloud),
								'QICE':(['Time','bottom_top','south_north','west_east'],i_qice),
								'TSK':(['Time','bottom_top','south_north','west_east'],i_tsk)
							},
							coords={
								'Time':i_Time,
								'bottom_top':np.arange(0,Nz,1.0),
								'south_north':np.arange(0,Nlat,1.0),
								'west_east':np.arange(0,Nlon,1.0)
							}
						)

icar_topo_ds.Time.attrs['units']='hours since 1900-01-01'
icar_topo_ds.Time.attrs['calendar']='gregorian'
#print icar_topo_ds.Time.attrs

icar_topo_ds.to_netcdf("./ideal_output/ideal_forcing.nc",format='NETCDF4')	
