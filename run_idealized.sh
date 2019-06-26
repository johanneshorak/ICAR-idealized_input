nzs=(
#       10
#       15
	23
#	30
#	43
#	60
)

printf "chosen option: $1\n"
outpath="/glade/u/home/horak/scratch1/sims/idealized/dtf_6.0/rh100/"
basepath="/glade/u/home/horak/scratch2/coding/icar-idealized_input"
for i in "${nzs[@]}"; do
	nz="$i"
        opt="$1"
        optfile1="icar-095.witch_${opt}.options"
	path1="w${opt}_lbcloud0-4kmrh100"
	#path2="w5c1a2bct2mp-1_lbcloud0-4kmrh100"

	printf "$outpath/${path1}_nz${nz}\t\t$optfile1\n"

	python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.25 --dlatf 0.25 --topo witch --a0 1800 --a1 13200 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.0 --rh 100 --rhprofile cloud_z0_0km_z1_4km --outdir "$outpath/$path1" --nzicar "$nz" --icaropt ./data/witch5_options/"$optfile1" && cd "$outpath/${path1}_nz${nz}" && qsub jobscript
	cd "$basepath"

        #python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.25 --dlatf 0.25 --topo witch --a0 1800 --a1 13200 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.0 --rh 100 --rhprofile cloud_z0_0km_z1_4km --outdir "$outpath/$path2" --nzicar "$nz" --icaropt ./data/witch5_options/icar-095.witch_5c1a2bct2mp-1.options && cd "$outpath/${path2}_nz${nz}" && qsub jobscript
        #cd "$basepath"

	#python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 100 --dlonf 0.25 --dlatf 0.25 --topo witch --a0 1800 --a1 13200 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270 --ffactor=0.0 --rh 100 --rhprofile cloud_z0_0km_z1_6km --outdir "$outpath/$path2" --nzicar "$nz" --icaropt ./data/witch5_options/icar-095.witch_5c1a2.options && cd "$outpath/${path2}_nz${nz}" && qsub jobscript
	#cd "$basepath"
done
