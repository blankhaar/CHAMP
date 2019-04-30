#module load intel
#
#ifort ./main_inv.f90 ./ned.f90 ./dipole.f90 ./rateq_inv.f90 ./initiate_read.f90 ./stimulate.f90 ./pumping.f90 ./trans_el.f90 ./expon.f90 -lmkl_lapack -fopenmp -o watned_inv_90.x

gfortran ./main_inv.f90 ./ned.f90 ./dipole.f90 ./rateq_inv.f90 ./initiate_read.f90 ./stimulate.f90 ./pumping.f90 ./trans_el.f90 ./expon.f90 -llapack -o watned_inv_90.x



