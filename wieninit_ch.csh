#!/bin/csh
#ciftostruct.pl $1.cif
ciftostructsym.pl $1.cif
#cartostruct.pl $1.car
x nn -c
x sgroup 
cp $1.struct $1.struct.org
cp $1.struct_sgroup $1.struct
x sgroup -settol 0.000001
cp $1.struct $1.struct.sgrp2
cp $1.struct_sgroup $1.struct

echo "EDIT STRUCT FILE FOR THE ATOMS WITH CH"
cp $1.struct $1.struct_before_CH
vi $1.struct
#instgen.pl $1.struct
instgen
init_lapw

# extend the energy region
makehe_wien.pl $1 


# make parallel script 
makestartwien.pl $1
