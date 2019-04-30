#module load intel
#
#ifort ./main_inv.f90 ./ned.f90 ./dipole.f90 ./zeta.f90 ./rateq_inv.f90  ./initiate_read.f90 ./stimulate.f90 ./trans_el.f90 ./analysis.f90 ./expon.f90 -lmkl_lapack -o watned_inv_92.x

gfortran ./main_inv.f90 ./ned.f90 ./dipole.f90 ./zeta.f90 ./rateq_inv.f90  ./initiate_read.f90 ./stimulate.f90 ./trans_el.f90 ./analysis.f90 ./expon.f90  -llapack -o watned_inv_92.x

