#a0="1800"
#a1="13200"
#a0="1000"
#a1="8000"
option="$1"
nz="$2"
a0="$3"
a1="$4"
dx="$5"
dy="$6"
Lx="$7"
Ly="$8"

outpath="/home/c707/c7071062/scratch/icar/results/idealized/dtf_6.0/rh0/a0_${a0}m_a1_${a1}m/"
path1="$option"
basepath="$( pwd )"
printf "\n"

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

printf "setting up idealized DRY simulation...\n"
printf "  grid-spacing        : dx = ${dx}\n"
printf "                      : dy = ${dy}\n"
printf "  domain extensions   : Lx = ${Lx}km\n"
printf "                        Ly = ${Ly}km\n"
printf "  bell shaped mountain: a0 = ${a0}m\n"
printf "                        a1 = ${a1}m\n"
printf "  option chosen       : ${option}\n"
printf "  z-levels            : ${nz}\n"
printf "\n"
printf "  writing to ${outpath}/${option}\n"
printf "\n"
sleep 5

mkdir -p "${outpath}/${option}"

python idealized_input.py --dx "${dx}" --dy "${dy}" --Lx "$Lx" --Ly "$Ly" --dlonf 0.25 --dlatf 0.25 --topo witch --a0 "${a0}" --a1 "${a1}" --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.0 --rh 0 --outdir "${outpath}/${option}/$path1" --nzicar "$nz" --icaropt ./data/witch5_options/icar-095."$option".options && cd "${outpath}/${option}/${path1}_nz${nz}" && qsub -l h_rt=6:00:00 jobscript
