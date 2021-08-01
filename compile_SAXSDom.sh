#!/bin/bash -e

echo " Start compile SAXSDom (will take ~10 min)"

cd /storage/hpc/data/jh7x3/SAXSDom//installation/Mocapy++-1.07

export LD_LIBRARY_PATH=/storage/hpc/data/jh7x3/SAXSDom//tools/boost_1_55_0/lib:$LD_LIBRARY_PATH

/storage/hpc/data/jh7x3/SAXSDom//tools/cmake-2.8.12.2/bin/cmake -DBOOST_ROOT='/storage/hpc/data/jh7x3/SAXSDom//tools/boost_1_55_0/' -DLAPACK_LIBRARY:FILEPATH='/storage/hpc/data/jh7x3/SAXSDom//tools/lapack-3.4.1/liblapack.a' .

make

cp /storage/hpc/data/jh7x3/SAXSDom//installation/Mocapy++-1.07/examples/SAXSDom  /storage/hpc/data/jh7x3/SAXSDom//bin

