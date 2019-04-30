#module load intel

#ifort ./main_para.f90 ./ned.f90 ./dipole.f90 ./integrals.f90 ./pumping.f90 ./zeta.f90 ./rateq_inv.f90 ./initiate.f90 ./stimulate.f90 ./trans_el.f90 ./analysis.f90 ./expon.f90 -lmkl_lapack -fopenmp  -o watned_inv_94.x


gfortran ./main_para.f90 ./ned.f90 ./dipole.f90 ./integrals.f90 ./pumping.f90 ./zeta.f90 ./rateq_inv.f90 ./initiate.f90 ./stimulate.f90 ./trans_el.f90 ./analysis.f90 ./expon.f90 -llapack -fopenmp -o watned_inv_94.x

