#!/bin/bash

if [ "$#" -ne 8 ]; then
  echo "Usage: $0 the number of parameter ($#) is not correct!" >&2

  exit 1
fi

GLOBAL_PATH=/storage/htc/bdm/jh7x3/DomainOrientation_project/SAXSDom/;
export LD_LIBRARY_PATH=$GLOBAL_PATH/tools/IMP2.6/lib:$GLOBAL_PATH/tools/boost_1_55_0/lib:$LD_LIBRARY_PATH

targetid=$1
seqfile=$2
saxsfile=$3
domainfile=$4
outputdir=$5
arguments=$6 # '  -t   -g test_assembly  -d 1 -x  3  --scoreWeight 0_800_0_0 --scoreWeightInitial 0_800_0_0  --scoreWeightFinal 1_100_1_1  --scoreWeightPenalyInitial 0.0001_0.0001_0.0001_0.0001  --scoreWeightPenalyFinal 10_20_10_10 --reg   --scoreFun KL'
epoch=$7
nativefile=$8

mkdir -p $outputdir

cd $outputdir

mkdir -p $outputdir/SCRATCH/
mkdir -p $outputdir/metapsicov/ 
python2 $GLOBAL_PATH/scripts/init_cm.py  --fasta ${seqfile}   > $outputdir/metapsicov/${targetid}_initial_domain.cm


echo "###### Generating $epoch decoys!!! \n\n"

if [ -f "${saxsfile}" ]
then
	echo "${saxsfile} found!\n";
else
	echo "${saxsfile} not exists, pass!\n";
  exit 1;
fi

if [ -f "./SCRATCH/${targetid}.ss8" ]
then
	echo "./SCRATCH/${targetid}.ss8  found!\n";
else
	echo "./SCRATCH/${targetid}.ss8 not exists, need generate!\n";
    $GLOBAL_PATH/tools/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh ${seqfile} $outputdir/SCRATCH/${targetid}
fi

`cp ${seqfile} ${targetid}.fasta`;

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

    echo "$GLOBAL_PATH/bin/SAXSDom  -i ${targetid}_siminit_Wsaxs_regularize -f  ${targetid}.fasta  -s SCRATCH/${targetid}.ss8    -c metapsicov/${targetid}_initial_domain.cm   -l $domainfile  -m $GLOBAL_PATH/lib/UniCon.iohmm    -n $nativefile  -o $outputdir/Assembly_docoy$decoy   $arguments\n\n" 
   $GLOBAL_PATH/bin/SAXSDom  -i ${targetid}_siminit_Wsaxs_regularize -f  ${targetid}.fasta  -s SCRATCH/${targetid}.ss8    -c metapsicov/${targetid}_initial_domain.cm   -l $domainfile  -m $GLOBAL_PATH/lib/UniCon.iohmm    -n $nativefile  -o $outputdir/Assembly_docoy$decoy   $arguments
   
   #rm $outputdir/Assembly_docoy$decoy/sample*
   rm $outputdir/Assembly_docoy$decoy/GlobalFoldon*pdb
   rm $outputdir/Assembly_docoy$decoy/*initial*pdb
}
