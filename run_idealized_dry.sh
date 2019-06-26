a0="1800"
a1="13200"
#a0="1000"
#a1="8000"
outpath="/glade/u/home/horak/scratch1/sims/idealized/dtf_6.0/rh0/a0_${a0}m_a1_${a1}m/"
option="$1"
path1="$option"
basepath="/glade/u/home/horak/scratch2/coding/icar-idealized_input/"
nz="$2"

mkdir -p "${outpath}/${option}"

python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.25 --dlatf 0.25 --topo witch --a0 "${a0}" --a1 "${a1}" --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.0 --rh 0 --outdir "${outpath}/${option}/$path1" --nzicar "$nz" --icaropt ./data/witch5_options/icar-095."$option".options && cd "${outpath}/${option}/${path1}_nz${nz}" && qsub jobscript
