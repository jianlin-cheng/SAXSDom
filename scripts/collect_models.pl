
$num = @ARGV;
if($num !=2)
{
  die "wrong parameters!\n";
}

$work_dir = $ARGV[0];
$outputdir = $ARGV[1];


opendir(DIR,"$work_dir") || die "Failed to open directory $work_dir\n";
@files = readdir(DIR);
closedir(DIR);

foreach $subdir (@files)
{
	if($subdir eq '.' or $subdir eq '..')
	{
		next;
	}	
	if(-e "$work_dir/$subdir/Assembly_docoy1/RcPutA_siminit_Wsaxs_regularize_000001.rebuilt.pdb")	
	{
		`cp $work_dir/$subdir/Assembly_docoy1/RcPutA_siminit_Wsaxs_regularize_000001.rebuilt.pdb $outputdir/SAXS_$subdir.pdb`;
	}
}

