if [ "$1" == "" ]; then
	printf " please supply number of vertical levels as parameter!\n"
	exit
fi

if [ "$2" == "" ]; then
	printf " please supply ICAR options file!\n"
	exit
fi

if [ "$3" == "" ]; then
	printf " please specify outdir!\n"
	exit
fi

nz="$1"
icaropt="$2"
outdir="$3"

printf " generating topography and forcing for $nz vertical levels\n"

# generate simulation with a 100% RH cloud between 0 and 3 km
python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.25 --dlatf 0.25 --topo witch --a0 1000 --a1 8000 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.0 --rh 100 --rhprofile cloud_z0_1km_z1_3km_below_rh50 --outdir "$outdir" --nzicar "$nz" --icaropt "$icaropt"

# with 100% moisture, this time we remove the forcing topography entirely (ffactor=0)
#python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.25 --dlatf 0.25 --topo witch --a0 1000 --a1 8000 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.0 --rh 100 --outdir "$outdir" --nzicar "$nz" --icaropt "$icaropt"

# with moisture
#python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.25 --dlatf 0.25 --topo witch --a0 1000 --a1 8000 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.4 --rh 100 --outdir "$outdir" --nzicar "$nz" --icaropt "$icaropt"
# without moisture
#python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.25 --dlatf 0.25 --topo witch --a0 1000 --a1 8000 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.4 --outdir "$outdir" --nzicar "$nz" --icaropt "$icaropt"

# without moisture, forcing resolution = high resolution topo
#python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.0089932 --dlatf 0.0089932 --topo witch --a0 1000 --a1 8000 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=1.0 --outdir "$outdir" --nzicar "$nz" --icaropt "$icaropt"

# without moisture, forcing resolution = high resolution topo and z-levels topo = z-levels ICAR
#python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.0089932 --dlatf 0.0089932 --topo witch --a0 1000 --a1 8000 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=1.0 --outdir "$outdir" --nzicar "$nz" --icaropt "$icaropt"
