#!/bin/bash

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 the number of parameter ($#) is not correct!" >&2

  exit 1
fi

export LD_LIBRARY_PATH=/home/jh7x3/Mocapy_tool/IMP2.6/lib:/home/jh7x3/Mocapy_tool/boost_1_55_0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/jh7x3/lib32/:/usr/lib64/:$LD_LIBRARY_PATH

targetid=$1
workdir=$2
outputdir=$3
arguments=$4 # '  -t   -g test_assembly  -d 1 -x  3  --scoreWeight 0_800_0_0 --scoreWeightInitial 0_800_0_0  --scoreWeightFinal 1_100_1_1  --scoreWeightPenalyInitial 0.0001_0.0001_0.0001_0.0001  --scoreWeightPenalyFinal 10_20_10_10 --reg   --scoreFun KL'
epoch=$5

cd $workdir
mkdir -p $outputdir
mkdir metapsicov/ 
python2 /home/casp13/SAXS_modeling/scripts/init_cm.py  --fasta ${targetid}.fasta   > metapsicov/${targetid}_initial_domain.cm

#sed -i 's/\/home\/saxs\/jie\/DomainOrientation_project/\/group\/birchler-cheng\/Jie_files\/DomainOrientation_project/g' ./domain_list_withPath_reindex

echo "###### Generating $epoch decoys!!! \n\n"


/home/jh7x3/Mocapy_tool/Mocapy++-1.07/examples/UniCon3D_DomainOrient_bySAXS_regularized_fox_V10_witheva -k 1  -q 0.5 -s 500 ${targetid}_reindex.pdb




cp ${targetid}_reindex.pdb.dat ${targetid}_reindex.dat  


if [ -f "./${targetid}_reindex.dat" ]
then
	echo "./${targetid}_reindex.dat  found!\n";
else
	echo "./${targetid}_reindex.dat not exists, pass!\n";
  exit 1;
fi

if [ -f "./SCRATCH/${targetid}.ss8" ]
then
	echo "./SCRATCH/${targetid}.ss8  found!\n";
else
	echo "./SCRATCH/${targetid}.ss8 not exists, need generate!\n";
    /home/jh7x3/tools/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh ${targetid}.fasta SCRATCH/${targetid}
fi


for ((decoy=1;decoy <= $epoch;decoy++))
{
    decoymodel=$outputdir/Assembly_docoy$decoy/${targetid}_siminit_Wsaxs_regularize_000001.rebuilt.pdb;
    if [ -f "$decoymodel" ]
    then
    	echo "$decoymodel found! Pass\n";
    	continue
    else
    	echo "start decoy $decoy\n\n"
    fi

    echo "/home/jh7x3/Mocapy_tool/Mocapy++-1.07/examples/UniCon3D_DomainOrient_bySAXS_regularized_fox_V12_almost_final  -i ${targetid}_siminit_Wsaxs_regularize -f  ${targetid}.fasta  -s SCRATCH/${targetid}.ss8    -c metapsicov/${targetid}_initial_domain.cm   -l ./domain_list_withPath_reindex  -m /home/casp13/Confold2-Unicon3D/UniCon3D/lib/UniCon.iohmm    -n ./${targetid}_reindex.pdb     -e ./${targetid}_reindex.dat -o $outputdir/Assembly_docoy$decoy   $arguments\n\n" 
   /home/jh7x3/Mocapy_tool/Mocapy++-1.07/examples/UniCon3D_DomainOrient_bySAXS_regularized_fox_V12_almost_final  -i ${targetid}_siminit_Wsaxs_regularize -f  ${targetid}.fasta  -s SCRATCH/${targetid}.ss8    -c metapsicov/${targetid}_initial_domain.cm   -l ./domain_list_withPath_reindex  -m /home/casp13/Confold2-Unicon3D/UniCon3D/lib/UniCon.iohmm    -n ./${targetid}_reindex.pdb     -e ./${targetid}_reindex.dat -o $outputdir/Assembly_docoy$decoy   $arguments
   
   rm $outputdir/Assembly_docoy$decoy/sample*
   rm $outputdir/Assembly_docoy$decoy/GlobalFoldon*pdb
   rm $outputdir/Assembly_docoy$decoy/*initial*pdb
}
