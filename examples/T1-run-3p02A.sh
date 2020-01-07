#!/bin/bash
#SBATCH -J  3p02A
#SBATCH -o 3p02A-%j.out
#SBATCH --partition Lewis,hpc5,hpc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=2G
#SBATCH --time 2-00:00

dtime=$(date +%m%d%y)

GLOBAL_PATH=/data/jh7x3/SAXSDom/;
export LD_LIBRARY_PATH=$GLOBAL_PATH/tools/IMP2.6/lib:$GLOBAL_PATH/tools/boost_1_55_0/lib:$LD_LIBRARY_PATH


mkdir -p /data/jh7x3/SAXSDom/test_out/3p02A_test/
cd /data/jh7x3/SAXSDom/test_out/3p02A_test/


echo /data/jh7x3/SAXSDom/examples/3p02A/3p02A_D1_prediction.pdb > /data/jh7x3/SAXSDom/test_out/3p02A_test/domain_pdb_list
echo /data/jh7x3/SAXSDom/examples/3p02A/3p02A_D2_prediction.pdb >> /data/jh7x3/SAXSDom/test_out/3p02A_test/domain_pdb_list

sh /data/jh7x3/SAXSDom/bin/run_SAXSDom_benchmark.sh 3p02A /data/jh7x3/SAXSDom/examples/3p02A/3p02A.fasta  /data/jh7x3/SAXSDom/examples/3p02A/saxs_profile.dat /data/jh7x3/SAXSDom/test_out/3p02A_test/domain_pdb_list /data/jh7x3/SAXSDom/test_out/3p02A_test/ 5 /data/jh7x3/SAXSDom/examples/3p02A/native.pdb  2>&1 | tee  /data/jh7x3/SAXSDom/test_out/3p02A_test.log

printf "\nFinished.."
printf "\nCheck log file </data/jh7x3/SAXSDom/test_out/3p02A_test.log>\n\n"



rm /data/jh7x3/SAXSDom/test_out/3p02A_test.running
