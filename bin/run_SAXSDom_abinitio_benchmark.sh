#!/bin/bash

if [ "$#" -ne 6 ]; then
  echo "Usage: $0 the number of parameter ($#) is not correct!" >&2

  exit 1
fi

GLOBAL_PATH=/data/jh7x3/SAXSDom/;
export LD_LIBRARY_PATH=$GLOBAL_PATH/tools/IMP2.6/lib:$GLOBAL_PATH/tools/boost_1_55_0/lib:$LD_LIBRARY_PATH

targetid=$1
seqfile=$2
domainfile=$3
outputdir=$4
epoch=50
ncpu=$5
nativefile=$6


mkdir -p $outputdir

cp $seqfile $outputdir/seq.fasta
cp $domainfile $outputdir/domainlist
cp $nativefile $outputdir/native.pdb

perl $GLOBAL_PATH/scripts/run_SAXSdom_abinitio_benchmark_parallel.pl $targetid $outputdir/seq.fasta  $outputdir/domainlist  $outputdir $epoch $outputdir/native.pdb $ncpu  2>&1 | tee $outputdir/run.log

printf "\nFinished.."


if [[ ! -f "$outputdir/${targetid}_SAXSDom_top1.pdb" ]];then 
	printf "!!!!! Failed to run SAXSDom on $targetid>\n\n"
else
	printf "\nJob successfully completed!"
	printf "\nResults: $outputdir/${targetid}_SAXSDom_top1.pdb\n\n"
fi

