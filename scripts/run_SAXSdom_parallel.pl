#!/usr/bin/perl -w

use Cwd 'abs_path';
use File::Basename;


$numArgs = @ARGV;
if($numArgs != 4)
{   
	print "the number of parameters is not correct!\n";
	exit(1);
}

$pdb_dir	= abs_path($ARGV[0]);  # 
$tool_dir	= "$ARGV[1]"; #
$scoreout	= "$ARGV[2]"; #
$proc_num 	= "$ARGV[3]"; #



$script_dir = abs_path(dirname($0));
$pulchar_program = "$tool_dir/pulchra_306/pulchra";
$scwrl4_program = "$tool_dir/scwrl4/Scwrl4";


if(!(-e $pulchar_program))
{
  die "Failed to find $pulchar_program\n";
}

if(!(-e $scwrl4_program))
{
  die "Failed to find $scwrl4_program\n";
}

if(!(-e "$tool_dir/qprob_package/bin/Qprob.sh"))
{
  die "Failed to find $tool_dir/qprob_package/bin/Qprob.sh\n";
}


opendir(DIR,"$pdb_dir") || die "failed to open directory $pdb_dir\n";

@files = readdir(DIR);

closedir(DIR);

%pdb2dfire=();

$shell_dir = "$pdb_dir/run_src";
if(-d $shell_dir)
{
	`rm -rf $shell_dir`;
	`mkdir $shell_dir`;
}else{
	`mkdir $shell_dir`;
}

$shell_indx = 0;
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
  
  #### run pulchar on pdb file 
  #print "$pulchar_program $pdbfile\n";
  `$pulchar_program $pdbfile`;
  if(!(-e "$filepath/$idname.rebuilt.pdb"))
  {
    die "The $filepath/$idname.rebuilt.pdb failed to be genearted\n";
  }
  
  #### run scwrl on pdb file 
  
  #print "$scwrl4_program -i $filepath/$idname.rebuilt.pdb -o $filepath/$idname.rebuilt.scwrl.pdb\n";
  `$scwrl4_program -i $filepath/$idname.rebuilt.pdb -o $filepath/$idname.rebuilt.scwrl.pdb`;
  if(!(-e "$filepath/$idname.rebuilt.scwrl.pdb"))
  {
    die "The $filepath/$idname.rebuilt.scwrl.pdb failed to be genearted\n";
  }
  
  #### run qprob on pdb file 
  if(!(-e "$pdb_dir/${idname}_qprob/$idname.Qprob_score"))
  {
		$shell_indx++;
		open(RUNFILE,">$shell_dir/job_$shell_indx.sh") || die "Failed to write $shell_dir/job_$shell_indx.sh\n\n";
		`touch $shell_dir/job_$shell_indx.queued`;
		print RUNFILE "#!/bin/bash\n\n";
		print RUNFILE "mv $shell_dir/job_$shell_indx.queued $shell_dir/job_$shell_indx.running\n\n";
		print RUNFILE "mkdir $pdb_dir/${idname}_qprob\n";
		print RUNFILE "cp $filepath/$idname.rebuilt.scwrl.pdb $filepath/${idname}_qprob/${idname}_scwrl.pdb\n"; 
		print RUNFILE "cd $pdb_dir/${idname}_qprob\n";
		print RUNFILE "perl $script_dir/pdb2fasta.pl $pdb_dir/${idname}_qprob/${idname}_scwrl.pdb $pdb_dir/${idname}_qprob/$idname $idname.fasta\n"; 

		print RUNFILE "mkdir models\n"; 
		print RUNFILE "cp $pdb_dir/${idname}_qprob/${idname}_scwrl.pdb models\n"; 
		print RUNFILE "printf \"$tool_dir/qprob_package/bin/Qprob.sh $pdb_dir/${idname}_qprob/$idname.fasta   $pdb_dir/${idname}_qprob/models  $pdb_dir/${idname}_qprob/ &> $pdb_dir/${idname}_qprob/run.log\\n\\n\"\n";
		print RUNFILE "$tool_dir/qprob_package/bin/Qprob.sh $pdb_dir/${idname}_qprob/$idname.fasta   $pdb_dir/${idname}_qprob/models  $pdb_dir/${idname}_qprob/ &> $pdb_dir/${idname}_qprob/run.log\n\n";
        #print RUNFILE "/storage/htc/bdm/tools/qprob_package/bin/Qprob.sh $pdb_dir/${idname}_qprob/$idname.fasta   $pdb_dir/${idname}_qprob/models  $pdb_dir/${idname}_qprob/ &> $pdb_dir/${idname}_qprob/run.log\n\n";
		print RUNFILE "mv $shell_dir/job_$shell_indx.running $shell_dir/job_$shell_indx.done";
		close RUNFILE;
  }else{
	print "$pdb_dir/${idname}_qprob/$idname.Qprob_score already generated\n";
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


open(OUT,">$scoreout") || die "Failed to open file $scoreout\n";
foreach $model (sort {$pdb2dfire{$b} <=> $pdb2dfire{$a}} keys %pdb2dfire)
{
  print OUT "$model\t".$pdb2dfire{$model}."\n";
}
close OUT;
