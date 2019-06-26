nz="$2"
option="$1"

outpath="/glade/u/home/horak/scratch1/sims/idealized/dtf_6.0/rh100/a0_1800m_a1_13200m/advection_modification/rh100_from_0km_to_4km/${option}"
mkdir -p "$outpath"
outdir="${outpath}/${option}"

printf " generating topography and forcing for $nz vertical levels\n"

# generate simulation with a 100% RH cloud between 0 and 3 km
python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.25 --dlatf 0.25 --topo witch --a0 1800 --a1 13200 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.0 --rh 100 --rhprofile cloud_z0_0km_z1_4km --outdir "$outdir" --nzicar "$nz" --icaropt "./data/witch5_options/icar-095.${option}.options"
#printf "changing to ${outdir}_nz${nz}\n"
cd "${outdir}_nz${nz}" && qsub jobscript
