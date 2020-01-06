#!/bin/bash -e

echo " Start compile SAXSDom (will take ~10 min)"

cd /storage/jhou4/Projects/SAXSDom/SAXSDom//installation/Mocapy++-1.07

export LD_LIBRARY_PATH=/storage/jhou4/Projects/SAXSDom/SAXSDom//tools/boost_1_55_0/lib:$LD_LIBRARY_PATH

/storage/jhou4/Projects/SAXSDom/SAXSDom//tools/cmake-2.8.12.2/bin/cmake -DBOOST_ROOT='/storage/jhou4/Projects/SAXSDom/SAXSDom//tools/boost_1_55_0/' -DLAPACK_LIBRARY:FILEPATH='/storage/jhou4/Projects/SAXSDom/SAXSDom//tools/lapack-3.4.1/liblapack.a' .

make

cp /storage/jhou4/Projects/SAXSDom/SAXSDom//installation/Mocapy++-1.07/examples/SAXSDom  /storage/jhou4/Projects/SAXSDom/SAXSDom//bin

