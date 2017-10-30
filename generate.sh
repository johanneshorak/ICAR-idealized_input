# generates an idealized witch of agnes topography
# python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 200 --topo witch --a0 1000 --a1 8000 --Nz 30 --ztop 30 --Nbv 0.01 --ws 10 --ws_angle 270

# generates a small forcing and topographie for testing with idealized witch of agnes topography and 100% rh upstream
# python idealized_input.py --dx 1 --dy 1 --Lx 10 --Ly 10 --topo witch --a0 1000 --a1 8000 --Nz 30 --ztop 30 --Nbv 0.01 --rh 100 --ws 10 --ws_angle 270

# generates an idealized witch of agnes topography with 100% rh upstream
python idealized_input.py --dx 1 --dy 1 --Lx 200 --Ly 200 --topo witch --a0 1000 --a1 8000 --Nz 30 --ztop 30 --Nbv 0.01 --rh 100 --ws 10 --ws_angle 270

