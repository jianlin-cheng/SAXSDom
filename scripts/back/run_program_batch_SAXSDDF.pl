
$num = @ARGV;
if($num !=4)
{
  die "wrong parameters!\n";
}

$outputdir = $ARGV[0];
$sbatchdir = $ARGV[1];
$epoch_start = $ARGV[2];
$epoch_end = $ARGV[3];


-d $outputdir || `mkdir -p $outputdir`;
-d $sbatchdir || `mkdir -p $sbatchdir`;
for($epoch =$epoch_start; $epoch<=$epoch_end;$epoch++)
{

  $outfile = "$sbatchdir/sbatch_$epoch.sh";
  open(OUT,">$outfile") || die "Failed to open $outfile\n";
  
  print OUT "#!/bin/bash -l\n";
  print OUT "#SBATCH -J epo$epoch\n";
  print OUT "#SBATCH -o epo$epoch.out\n";
  print OUT "#SBATCH --partition Lewis,hpc4,hpc5\n";
  print OUT "#SBATCH --nodes=1\n";
  print OUT "#SBATCH --ntasks=1         # leave at '1' unless using a MPI code\n";
  print OUT "#SBATCH --cpus-per-task=1  # cores per task\n";
  print OUT "#SBATCH --mem-per-cpu=10G  # memory per core (default is 1GB/core)\n";
  print OUT "#SBATCH --time 2-00:00     # days-hours:minutes\n";
  
  -d "$outputdir/epoch$epoch" || `mkdir $outputdir/epoch$epoch`;
  
  `cp -ar /data/jh7x3/SAXS_fold_project/SASDDF9/SCRATCH/ $outputdir/epoch$epoch`;
  print OUT "/storage/htc/bdm/jh7x3/DomainOrientation_project/SAXSDom/bin/run_SAXSDom_withNative.sh SASDDF9 /data/jh7x3/SAXS_fold_project/SASDDF9/SASDDF9.fasta  /data/jh7x3/SAXS_fold_project/SASDDF9/SASDDF9-A-new.dat /data/jh7x3/SAXS_fold_project/SASDDF9/domain_models/domain_pdb_list  $outputdir/epoch$epoch  ' -t   -g test_assembly  -d 1 -x  2  --scoreWeight 30_0_0_0 --scoreWeightInitial 30_0_0_0  --scoreCombine' 1 /data/jh7x3/SAXS_fold_project/SASDDF9/SASDDF9_native.pdb\n";
  
  close OUT;
  
  chdir($sbatchdir);
  sleep(2);
  `sh sbatch_$epoch.sh &> sbatch_$epoch.log &`;

}

