#!/bin/bash

if [ "$#" -ne 6 ]; then
  echo "Usage: $0 the number of parameter ($#) is not correct!" >&2

  exit 1
fi

GLOBAL_PATH=/data/jh7x3/SAXSDom/;
export LD_LIBRARY_PATH=$GLOBAL_PATH/tools/IMP2.6/lib:$GLOBAL_PATH/tools/boost_1_55_0/lib:$LD_LIBRARY_PATH

targetid=$1
seqfile=$2
saxsfile=$3
domainfile=$4
outputdir=$5
epoch=$6


mkdir -p $outputdir
cd $outputdir


cp $seqfile $outputdir/seq.fasta
cp $saxsfile $outputdir/seq.dat
cp $domainfile $outputdir/domainlist

perl $GLOBAL_PATH/scripts/run_SAXSdom_parallel.pl $targetid $outputdir/seq.fasta $outputdir/seq.dat  $outputdir/domainlist  $outputdir $epoch  2>&1 | tee $outputdir/run.log

printf "\nFinished.."


if [[ ! -f "$outputdir/${targetid}_SAXSDom_top1.pdb" ]];then 
	printf "!!!!! Failed to run SAXSDom on $targetid>\n\n"
else
	printf "\nJob successfully completed!"
	printf "\nResults: $outputdir/${targetid}_SAXSDom_top1.pdb\n\n"
fi

