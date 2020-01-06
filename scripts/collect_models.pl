
$num = @ARGV;
if($num !=3)
{
  die "wrong parameters!\n";
}
$targetid = $ARGV[0];
$work_dir = $ARGV[1];
$outputdir = $ARGV[2];


opendir(DIR,"$work_dir") || die "Failed to open directory $work_dir\n";
@files = readdir(DIR);
closedir(DIR);

if(!(-d $outputdir))
{
	`mkdir -p $outputdir`;
}
$index=0;
foreach $subdir (@files)
{
	if($subdir eq '.' or $subdir eq '..')
	{
		next;
	}	
	if(-e "$work_dir/$subdir/${targetid}_saxsdom_000001.rebuilt.pdb")	
	{
		$index++;
		`cp $work_dir/$subdir/${targetid}_saxsdom_000001.rebuilt.pdb $outputdir/SAXS_model$index.pdb`;
	}
}

