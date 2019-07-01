#a0="1000"
#a1="8000"
ws="10"
N="0.005"
option="$1"
nz="$2"
a0="$3"
a1="$4"
dx="$5"
dy="$6"
Lx="$7"
Ly="$8"

outpath="/home/c707/c7071062/scratch/icar/results/idealized/dtf_6.0/rh100/a0_${a0}m_a1_${a1}m/advection_modification/rh100_from_0km_to_4km/ws${ws}_N${N}/${option}"
outdir="${outpath}/${option}"

if [ "$option" = "" ]; then
	printf " usage: $0 option nz [a0] [a1] [dx] [dy] [Lx] [Ly]\n"
	printf "   all except the first two parameters are optional\n"
	exit 1
fi

if [ "$a0" = "" ]; then
	a0="1000"
	a1="8000"
fi

if [ "$dx" = "" ]; then
	dx="1"
	dy="1"
fi

if [ "$Lx" = "" ]; then
	Lx="200"
fi

if [ "$Ly" = "" ]; then
	Ly="100"
fi

printf "setting up idealized WET simulation...\n"
printf "  grid-spacing        : dx = ${dx}\n"
printf "                      : dy = ${dy}\n"
printf "  domain extensions   : Lx = ${Lx}km\n"
printf "                        Ly = ${Ly}km\n"
printf "  bell shaped mountain: a0 = ${a0}m\n"
printf "                        a1 = ${a1}m\n"
printf "  option chosen       : ${option}\n"
printf "  z-levels            : ${nz}\n"
printf "\n"
printf "  writing to ${outdir}_nz${nz}\n"
printf "\n"
printf "NOTE: windspeed ${ws}m/s!\n"
printf "      N = ${N}\n"
sleep 5

mkdir -p "$outpath"

printf " generating topography and forcing for $nz vertical levels\n"

# generate simulation with a 100% RH cloud between 0 and 3 km
python idealized_input.py --dx "${dx}" --dy "${dy}" --Lx "${Lx}" --Ly "${Ly}" --dlonf 0.25 --dlatf 0.25 --topo witch --a0 "${a0}" --a1 "${a1}" --Nz 30 --ztop 30 --Nbv "${N}" --ws "${ws}" --ws_angle 270 --ffactor=0.0 --rh 100 --rhprofile cloud_z0_0km_z1_4km --outdir "$outdir" --nzicar "$nz" --icaropt "./data/witch5_options/icar-095.${option}.options"
cd "${outdir}_nz${nz}" && qsub jobscript
