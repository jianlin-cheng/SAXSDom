
$num = @ARGV;
if($num !=3)
{
  die "wrong parameters!\n";
}

$outputdir = $ARGV[0];
$sbatchdir = $ARGV[1];
$epoch_all = $ARGV[2];


-d $outputdir || `mkdir -p $outputdir`;
-d $sbatchdir || `mkdir -p $sbatchdir`;
for($epoch =1; $epoch<=$epoch_all;$epoch++)
{

  $outfile = "$sbatchdir/sbatch_$epoch.sh";
  open(OUT,">$outfile") || die "Failed to open $outfile\n";
  
  print OUT "#!/bin/bash -l\n";
  print OUT "#SBATCH -J epo$epoch\n";
  print OUT "#SBATCH -o epo$epoch.out\n";
  print OUT "#SBATCH --partition hpc4,hpc5\n";
  print OUT "#SBATCH --nodes=1\n";
  print OUT "#SBATCH --ntasks=1         # leave at '1' unless using a MPI code\n";
  print OUT "#SBATCH --cpus-per-task=1  # cores per task\n";
  print OUT "#SBATCH --mem-per-cpu=20G  # memory per core (default is 1GB/core)\n";
  print OUT "#SBATCH --time 2-00:00     # days-hours:minutes\n";
  
  -d "$outputdir/epoch$epoch" || `mkdir $outputdir/epoch$epoch`;
  
  `cp -ar /storage/htc/bdm/jh7x3/DomainOrientation_project/SAXSDom/benchmark/RcPutA_modeling/SCRATCH/ $outputdir/epoch$epoch`;
  print OUT "/storage/htc/bdm/jh7x3/DomainOrientation_project/SAXSDom/bin/run_SAXSDom_withNative.sh RcPutA /storage/htc/bdm/jh7x3/DomainOrientation_project/SAXSDom/benchmark/RcPutA/multicom_prediction/RcPutA_full.fasta  /storage/htc/bdm/jh7x3/DomainOrientation_project/SAXSDom/benchmark/RcPutA/multicom_prediction/RcH7m23.dat /storage/htc/bdm/jh7x3/DomainOrientation_project/SAXSDom/benchmark/RcPutA/multicom_prediction/domains_examples/link12/domain_pdb_list  $outputdir/epoch$epoch  ' -t   -g test_assembly  -d 1 -x  1  --scoreWeight 30_700_700_700 --scoreWeightInitial 30_700_700_700  --scoreCombine' 1 /storage/htc/bdm/jh7x3/DomainOrientation_project/SAXSDom/benchmark/RcPutA/multicom_prediction/meta_prc2.pdb\n";
  
  close OUT;
  
  chdir($sbatchdir);
  `sbatch sbatch_$epoch.sh`;

}

