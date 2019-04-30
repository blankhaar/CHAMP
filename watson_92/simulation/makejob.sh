#!/bin/csh -f

foreach vth (1 2 3)

if (! -e ./v$vth ) then
  mkdir ./v$vth
endif

foreach deg (0 45 75)
if (! -e ./v$vth/deg$deg/ ) then
  mkdir ./v$vth/deg$deg
endif


if ( $vth == 1 ) then

cat << EOF > ./v$vth/deg$deg/in_watned.dat
input file SiO 1-0
3

5 4 -33375.709d0 -0.79d0 -1.34d0 1.798E-9 2E7 0.02d0 1.d0 1.d0
6 5 0.d0 3.71d0 4.12d0 1.806E-9 2E7 0.02d0 1.d0 1.d0
7 6 43017.581d0 6.51d0 7.24d0 1.860E-9 2E7 0.02d0 1.d0 1.d0

22.2350798E9
20.d0 $deg 
60 0.6d0 3.d0 18.01528d0
0.1d0 5E9 1499
1.d0

EOF

endif

if ( $vth == 2 ) then

cat << EOF > ./v$vth/deg$deg/in_watned.dat
input file SiO 1-0
3

5 4 -33375.709d0 -0.79d0 -1.34d0 1.798E-9 2E7 0.02d0 1.d0 1.d0
6 5 0.d0 3.71d0 4.12d0 1.806E-9 2E7 0.02d0 1.d0 1.d0
7 6 43017.581d0 6.51d0 7.24d0 1.860E-9 2E7 0.02d0 1.d0 1.d0

22.2350798E9
20.d0 $deg 
60 1.0d0 3.d0 18.01528d0
0.1d0 5E9 1499
1.d0

EOF

endif

if ( $vth == 3 ) then

cat << EOF > ./v$vth/deg$deg/in_watned.dat
input file SiO 1-0
3

5 4 -33375.709d0 -0.79d0 -1.34d0 1.798E-9 2E7 0.02d0 1.d0 1.d0
6 5 0.d0 3.71d0 4.12d0 1.806E-9 2E7 0.02d0 1.d0 1.d0
7 6 43017.581d0 6.51d0 7.24d0 1.860E-9 2E7 0.02d0 1.d0 1.d0

22.2350798E9
20.d0 $deg 
60 2.0d0 3.d0 18.01528d0
0.1d0 5E9 1499
1.d0

EOF

endif

cp watned_inv_92.x ./v$vth/deg$deg

cd ./v$vth/deg$deg
./watned_inv_92.x
\rm watned_inv_92.x
cd ../../

end
end


