#!/usr/bin/perl -w

use Cwd 'abs_path';
use File::Basename;

$GLOBAL_PATH="SOFTWARE_PATH";

$numArgs = @ARGV;
if($numArgs < 7)
{   
	print "the number of parameters is not correct!\n";
	exit(1);
}

$proc_num = 1;

$targetid="$ARGV[0]";
$seqfile="$ARGV[1]";
$saxsfile="$ARGV[2]";
$domainfile="$ARGV[3]";
$outputdir="$ARGV[4]";
$epoch="$ARGV[5]";
$nativefile="$ARGV[6]";
$proc_num="$ARGV[7]";

if($proc_num>10)
{
	$proc_num = 10;
}

$script_dir = abs_path(dirname($0));

if(!(-e $seqfile))
{
  die "Failed to find $seqfile\n";
}

if(!(-e $saxsfile))
{
  die "Failed to find $saxsfile\n";
}

if(!(-e "$GLOBAL_PATH/tools/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh"))
{
  die "Failed to find $GLOBAL_PATH/tools/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh\n";
}


`mkdir -p $outputdir`;

chdir($outputdir);

`mkdir -p $outputdir/SCRATCH/`;
`mkdir -p $outputdir/metapsicov/`;
`python2 $GLOBAL_PATH/scripts/init_cm.py  --fasta ${seqfile}   > $outputdir/metapsicov/${targetid}_initial_domain.cm`;

if(!(-e "./SCRATCH/${targetid}.ss8"))
{
	print "./SCRATCH/${targetid}.ss8  found!\n";
}else
{
	print "./SCRATCH/${targetid}.ss8 not exists, need generate!\n";
    `$GLOBAL_PATH/tools/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh ${seqfile} $outputdir/SCRATCH/${targetid}`;
}

`cp ${seqfile} ${targetid}.fasta`;



opendir(DIR,"$pdb_dir") || die "failed to open directory $pdb_dir\n";

@files = readdir(DIR);

closedir(DIR);

%pdb2dfire=();

$shell_dir = "$outputdir/run_src";
if(-d $shell_dir)
{
	`rm -rf $shell_dir`;
	`mkdir $shell_dir`;
}else{
	`mkdir $shell_dir`;
}

$shell_indx = 0;
for($decoy=1;$decoy <= $epoch;$decoy++)
{

  
  $decoymodel="$outputdir/Assembly_docoy$decoy/${targetid}_saxsdom_000001.rebuilt.pdb";

  #### run qprob on pdb file 
  if(!(-e "$decoymodel"))
  {
  	print "start generate $decoymodel found!\n";

	$shell_indx++;
	open(RUNFILE,">$shell_dir/job_$shell_indx.sh") || die "Failed to write $shell_dir/job_$shell_indx.sh\n\n";
	`touch $shell_dir/job_$shell_indx.queued`;
	print RUNFILE "#!/bin/bash\n\n";
	print RUNFILE "mv $shell_dir/job_$shell_indx.queued $shell_dir/job_$shell_indx.running\n\n";
	print RUNFILE "mkdir $outputdir/Assembly_docoy$decoy\n";
	print RUNFILE "printf \"$GLOBAL_PATH/bin/SAXSDom  -i ${targetid}_saxsdom -f  ${targetid}.fasta  -s SCRATCH/${targetid}.ss8    -c metapsicov/${targetid}_initial_domain.cm   -l $domainfile  -m $GLOBAL_PATH/lib/UniCon.iohmm        -e $saxsfile -o $outputdir/Assembly_docoy$decoy -t   -g test_assembly  -d 1 -x  1  --scoreWeight 10_700_700_700 --scoreWeightInitial 10_700_700_700  --scoreCombine  -n $nativefile$GLOBAL_PATH/bin/SAXSDom  -i ${targetid}_saxsdom -f  ${targetid}.fasta  -s SCRATCH/${targetid}.ss8    -c metapsicov/${targetid}_initial_domain.cm   -l $domainfile  -m $GLOBAL_PATH/lib/UniCon.iohmm        -e $saxsfile -o $outputdir/Assembly_docoy$decoy -t   -g test_assembly  -d 1 -x  1  --scoreWeight 10_700_700_700 --scoreWeightInitial 10_700_700_700  --scoreCombine  -n $nativefile &> $outputdir/Assembly_docoy$decoy/run.log\\n\\n\"\n";
	print RUNFILE "$GLOBAL_PATH/bin/SAXSDom  -i ${targetid}_saxsdom -f  ${targetid}.fasta  -s SCRATCH/${targetid}.ss8    -c metapsicov/${targetid}_initial_domain.cm   -l $domainfile  -m $GLOBAL_PATH/lib/UniCon.iohmm        -e $saxsfile -o $outputdir/Assembly_docoy$decoy -t   -g test_assembly  -d 1 -x  1  --scoreWeight 10_700_700_700 --scoreWeightInitial 10_700_700_700  --scoreCombine  -n $nativefile$GLOBAL_PATH/bin/SAXSDom  -i ${targetid}_saxsdom -f  ${targetid}.fasta  -s SCRATCH/${targetid}.ss8    -c metapsicov/${targetid}_initial_domain.cm   -l $domainfile  -m $GLOBAL_PATH/lib/UniCon.iohmm        -e $saxsfile -o $outputdir/Assembly_docoy$decoy -t   -g test_assembly  -d 1 -x  1  --scoreWeight 10_700_700_700 --scoreWeightInitial 10_700_700_700  --scoreCombine  -n $nativefile &> $outputdir/Assembly_docoy$decoy/run.log\n\n";
	print RUNFILE "rm $outputdir/Assembly_docoy$decoy/sample*";
	print RUNFILE "rm $outputdir/Assembly_docoy$decoy/GlobalFoldon*pdb";
	print RUNFILE "rm $outputdir/Assembly_docoy$decoy/*initial*pdb";
	print RUNFILE "mv $shell_dir/job_$shell_indx.running $shell_dir/job_$shell_indx.done";
	close RUNFILE;
  }else{
	print "$decoymodel found! Pass\n";
  }
}


##########################  Submiting jobs in parallel

chdir($shell_dir);
opendir(DIR,"$shell_dir") || die "Failed to open directory $shell_dir\n";
@input_files = readdir(DIR);
closedir(DIR);

@running_files = ();
foreach $file (sort @input_files)
{
	if($file eq '.' or $file eq '..' or substr($file,length($file)-3) ne '.sh')
	{
		next;
	}
	$file_path = "$shell_dir/$file";
	push @running_files,$file_path;
}
	
foreach $file_path (sort @running_files)
{
	## check the running jobs
	$min_elaps=0;
	while(1)
	{
		opendir(DIR,"$shell_dir") || die "Failed to open directory $shell_dir\n";
		@out_files = readdir(DIR);
		closedir(DIR);
		
		$running_num = 0;
		foreach $check_file (sort @out_files)
		{
			if($check_file eq '.' or $check_file eq '..' or substr($check_file,length($check_file)-8) ne '.running')
			{
				next;
			}
			$running_num++;
		}
		if($running_num<$proc_num)
		{
			last;
		}
		sleep(60);
		$min_elaps++;
		if($min_elaps > 60)
		{
			last; # move to next;
		}
	}
	
	if(!(-e substr($file_path,0,length($file_path)-3).".done"))
	{
		print "run test $file_path\n";
		system("sh $file_path &> $file_path.log &");
	}else{
		print "$file_path has been done\n";
		$queue_file = substr($file_path,0,length($file_path)-3).".queued";
		if(-e $queue_file)
		{
			`rm $queue_file`;
		}
	}
	
	$running_jobs=0;
	$processed_jobs=0;
	opendir(DIR,"$shell_dir") || die "Failed to open directory $shell_dir\n";
	@out_files = readdir(DIR);
	closedir(DIR);
	foreach $check_file (sort @out_files)
	{
		if($check_file eq '.' or $check_file eq '..')
		{
			next;
		}
		if(substr($check_file,length($check_file)-5) eq '.done')
		{
			$processed_jobs++;
		}
		if(substr($check_file,length($check_file)-8) eq '.running')
		{
			$running_jobs++;
		}
	}
	$remain_jobs = @running_files-$processed_jobs-$running_jobs;
	print "Current running jobs ($running_num), processed jobs ($processed_jobs), unprocessed jobs ($remain_jobs)\n\n";
	sleep(5);
}

#### check if all files have finished
print "#### check if all files have finished\n";

while(1)
{

	opendir(DIR,"$shell_dir") || die "Failed to open directory $shell_dir\n";
	@out_files = readdir(DIR);
	closedir(DIR);

  $running_num = 0;
  foreach $check_file (sort @out_files)
  {
  	if($check_file eq '.' or $check_file eq '..' or substr($check_file,length($check_file)-3) eq '.sh')
  	{
  		next;
  	}
   
    if(substr($check_file,length($check_file)-8) eq '.running' or substr($check_file,length($check_file)-7) eq '.queued')
    {
  	  $running_num++;
    }
  }
  
  if($running_num>0)
  {
    print "$running_num jobs are still running, please wait\n";
  }else{
    print "All training jobs are done\n\n";
    last;
  }
  
  sleep(60*5);
  
}


##### summarize results

### collect all models for evalution
`perl $GLOBAL_PATH/scripts/collect_models.pl $targetid  $outputdir  $outputdir/all_models`;

### run qprob to rank the model
#`$GLOBAL_PATH//tools/DeepQA/tools/qprob_package/bin/Qprob.sh $outputdir/${targetid}.fasta  $outputdir/all_models/ $outputdir/all_models_qprob`;



=pod
foreach $file (sort @files)
{

  if($file eq '.' or $file eq '..' or index($file,'.pdb')<0 or index($file,'scwrl')>0 or index($file,'rebuilt')>0)
  {
    next;
  }  
  
  $pdbfile = "$pdb_dir/$file";
  
  if(index($pdbfile,'/')>=0)
  {
    @tmp = split(/\//,$pdbfile);
    $idname = pop @tmp;
    $filepath = join('/',@tmp);
  }
  if(index($idname,'.pdb')>0)
  {
    $idname = substr($idname,0,index($idname,'.pdb'));
  }
  

  if(!(-e "$pdb_dir/${idname}_qprob/$idname.Qprob_score"))
  {
    die "The $pdb_dir/${idname}_qprob/$idname.Qprob_score failed to be genearted\n";
  }
  
  $qprob_score=10000;  #initialize
  open(RWPLUS_CHECK, "$pdb_dir/${idname}_qprob/$idname.Qprob_score") || print "Can't open qprob output file.\n";
  while(<RWPLUS_CHECK>)
  {
      $line = $_;
      $line =~ s/\n//;
    @tem_split=split(/\s+/,$line);
    $qprob_score=$tem_split[1];
  }
  close RWPLUS_CHECK;
  
  print "qprob score of $idname: $qprob_score\n";
  $pdb2dfire{$file} =  $qprob_score;
  `rm $filepath/$idname.rebuilt.pdb`;
  `rm $filepath/$idname.rebuilt.scwrl.pdb`;
}

=cut


