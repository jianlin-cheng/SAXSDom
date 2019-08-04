#!/usr/bin/perl -w
 use FileHandle; # use FileHandles instead of open(),close()
 use Cwd;
 use Cwd 'abs_path';

######################## !!! customize settings here !!! ############################
#																					#
# Set directory of SAXSDom databases and tools								        #

$SAXSDom_db_tools_dir = "/storage/htc/bdm/tools/SAXSDom_database_tools/";							        

######################## !!! End of customize settings !!! ##########################

######################## !!! Don't Change the code below##############


$install_dir = getcwd;
$install_dir=abs_path($install_dir);


if(!-d $install_dir)
{
	die "The SAXSDom directory ($install_dir) is not existing, please revise the customize settings part inside the configure.pl, set the path as  your unzipped SAXSDom directory\n";
}

if(!-d $SAXSDom_db_tools_dir)
{
	die "The SAXSDom databases/tools folder ($SAXSDom_db_tools_dir) is not existing\n";
}

if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
        $install_dir .= "/";
}

if ( substr($SAXSDom_db_tools_dir, length($SAXSDom_db_tools_dir) - 1, 1) ne "/" )
{
        $SAXSDom_db_tools_dir .= "/";
}



print "checking whether the configuration file run in the installation folder ...";
$cur_dir = `pwd`;
chomp $cur_dir;
$configure_file = "$cur_dir/configure.pl";
if (! -f $configure_file || $install_dir ne "$cur_dir/")
{
        die "\nPlease check the installation directory setting and run the configure program under the main directory of SAXSDom.\n";
}
print " OK!\n";



if (! -d $install_dir)
{
	die "can't find installation directory.\n";
}
if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
	$install_dir .= "/"; 
}


######### check the SAXSDom database and tools

$database_dir = "$SAXSDom_db_tools_dir/databases";
$tools_dir = "$SAXSDom_db_tools_dir/tools";

if(!(-d $database_dir) or !(-d $tools_dir))
{
	die "Failed to find databases and tools under $SAXSDom_db_tools_dir/\n";
}

if($SAXSDom_db_tools_dir eq "$cur_dir/")
{
	die "Same directory as SAXSDom main folder. Differnt path for original databases/tools folder $SAXSDom_db_tools_dir is recommended.\n";
}
#create link for databases and tools
`rm ${install_dir}databases`; 
`rm ${install_dir}tools`; 
`ln -s $database_dir ${install_dir}databases`;
`ln -s $tools_dir ${install_dir}tools`;


if (prompt_yn("SAXSDom will be installed into <$install_dir> ")){

}else{
	die "The installation is cancelled!\n";
}
print "Start install SAXSDom into <$install_dir>\n"; 




print "#########  (1) Configuring option files\n";

$option_list = "$install_dir/installation/configure_list";

if (! -f $option_list)
{
        die "\nOption file $option_list not exists.\n";
}
configure_file2($option_list,'bin');
configure_file2($option_list,'installation');
print "#########  Configuring option files, done\n\n\n";


### compress benchmark dataset

$benchmark_dir = "$install_dir/installation";
chdir($benchmark_dir);

print "Updating benchmark dataset in $benchmark_dir\n\n";
if(-e "benchmark.tar.gz")
{
	`rm benchmark.tar.gz`;
}
`wget http://sysbio.rnet.missouri.edu/bdm_download/test/benchmark.tar.gz`;
`tar -zxf benchmark.tar.gz`;
`rm benchmark.tar.gz`;
if(! -d "benchmark")
{
	die "Failed to download benchmark.tar.gz from http://sysbio.rnet.missouri.edu/bdm_download/test/benchmark.tar.gz, please contact chengji\@missouri.edu\n";
}


system("mv $install_dir/installation/SAXSDom_test_codes/T0_run_SAXSDom*.sh $install_dir/examples");
system("chmod +x $install_dir/examples/*.sh");


=pod
print "#########  (3) Configuring scripts\n";

$option_list = "$install_dir/installation/SAXSDom_configure_files/SAXSDom_scripts_list";

if (! -f $option_list)
{
        die "\nOption file $option_list not exists.\n";
}
configure_file($option_list,'src');
print "#########  Configuring scripts, done\n\n\n";



print "#########  (4) Configuring examples\n";

$option_list = "$install_dir/installation/SAXSDom_configure_files/SAXSDom_examples_list";

if (! -f $option_list)
{
        die "\nOption file $option_list not exists.\n";
}
system("rm $install_dir/installation/SAXSDom_test_codes/*.sh");
configure_file2($option_list,'installation');
print "#########  Configuring examples, done\n\n\n";

system("chmod +x $install_dir/installation/SAXSDom_test_codes/*sh");



system("cp $install_dir/src/run_SAXSDom.sh $install_dir/bin/run_SAXSDom.sh");
system("chmod +x $install_dir/bin/run_SAXSDom.sh");



print "#########  (5) Configuring database update scripts\n";

$option_list = "$install_dir/installation/SAXSDom_configure_files/SAXSDom_db_list";

if (! -f $option_list)
{
        die "\nOption file $option_list not exists.\n";
}
configure_file2($option_list,'src');
print "#########  Configuring database update scripts, done\n\n\n";






##### generate scripts for methods, saved in bin
#installation/SAXSDom_programs/.P1_run_hhsearch.sh

print "#########  (7) Configuring SAXSDom programs\n";
$method_file = "$install_dir/method.list";
$option_list = "$install_dir/installation/SAXSDom_configure_files/SAXSDom_programs_list";

`rm $install_dir/installation/SAXSDom_programs/*sh`;
`rm $install_dir/bin/*sh`;

$python_env = 0;
$boost_enable = 0;
open(OUT,">$install_dir/Method_tutorial.txt") || die "Failed to open file $install_dir/Method_tutorial.txt\n";
if(!(-e $method_file) or !(-e $option_list))
{
	print "\nFailed to find method file ($method_file and $option_list), please contact us!\n\n";
}else{
	open(IN,$option_list) || die "Failed to open file $option_list\n";
	@contents = <IN>;
	close IN;
	%method_programs=();
	foreach $line (@contents)
	{
		chomp $line;
		if(substr($line,0,1) eq '#')
		{
			next;
		}
		$line =~ s/^\s+|\s+$//g;
		if($line eq '')
		{
			next;
		}
		@tmp = split(':',$line);
		$method_programs{$tmp[0]} = $tmp[1];
	}
	
	open(IN,$method_file) || die "Failed to open file $method_file\n";
	@contents = <IN>;
	open(TMP,">$install_dir/installation/SAXSDom_configure_files/option.tmp");
	foreach $method (@contents)
	{
		chomp $method;
		if(substr($method,0,1) eq '#')
		{
			next;
		}
		$method =~ s/^\s+|\s+$//g;
		if($method eq '')
		{
			next;
		}
		if(exists($method_programs{"${method}"}))
		{
			$file = $method_programs{"${method}"};
			print TMP "$file\n";
		}
		if(exists($method_programs{"${method}_easy"}))
		{
			$file = $method_programs{"${method}_easy"};
			print TMP "$file\n";
		}
		if(exists($method_programs{"${method}_hard"}))
		{
			$file = $method_programs{"${method}_hard"};
			print TMP "$file\n";
		}
	}
	close TMP;
	configure_file2("$install_dir/installation/SAXSDom_configure_files/option.tmp",'installation');
	`rm $install_dir/installation/SAXSDom_configure_files/option.tmp`;
	`cp $install_dir/installation/SAXSDom_programs/*sh $install_dir/bin/`;
	
	print "#########  Configuring examples, done\n\n\n";
	
	$method_indx = 0;
	foreach $method (@contents)
	{
		chomp $method;
		if(substr($method,0,1) eq '#')
		{
			next;
		}
		$method =~ s/^\s+|\s+$//g;
		if($method eq '')
		{
			next;
		}
		$method_indx++;
		
		print  "\n################################################################# Method $method_indx: $method  #################################################################\n\n";
		print  OUT "\n################################################################# Method $method_indx: $method  #################################################################\n\n";
		if(exists($method_programs{"${method}"}))
		{
			$file = $method_programs{"${method}"};
			@tmp = split(/\//,$file);
			$program_file = pop @tmp;
			if(-e "$install_dir/bin/$program_file")
			{
				print "Usage: $install_dir/bin/$program_file <target id> <fasta> <output-directory>\n\n";
				print OUT "Usage: $install_dir/bin/$program_file <target id> <fasta> <output-directory>\n\n";
				print "\t** Example: $install_dir/bin/$program_file T1006 $install_dir/examples/T1006.fasta $install_dir/test_out/T1006_$method\n\n";
				print OUT "\t** Example: $install_dir/bin/$program_file T1006 $install_dir/examples/T1006.fasta $install_dir/test_out/T1006_$method\n\n";
			}
			
		}
		if(exists($method_programs{"${method}_easy"}))
		{
			$file = $method_programs{"${method}_easy"};
			@tmp = split(/\//,$file);
			$program_file = pop @tmp;
			if(-e "$install_dir/bin/$program_file")
			{
				print "Usage: $install_dir/bin/$program_file <target id> <fasta> <output-directory>\n\n";
				print OUT "Usage: $install_dir/bin/$program_file <target id> <fasta> <output-directory>\n\n";
				print "\t** Example: $install_dir/bin/$program_file T1006 $install_dir/examples/T1006.fasta $install_dir/test_out/T1006_$method\n\n";
				print OUT "\t** Example: $install_dir/bin/$program_file T1006 $install_dir/examples/T1006.fasta $install_dir/test_out/T1006_$method\n\n";
			}
		}
		if(exists($method_programs{"${method}_hard"}))
		{
			$file = $method_programs{"${method}_hard"};
			@tmp = split(/\//,$file);
			$program_file = pop @tmp;
			if(-e "$install_dir/bin/$program_file")
			{
				print "Usage: $install_dir/bin/$program_file <target id> <fasta> <output-directory>\n\n";
				print OUT "Usage: $install_dir/bin/$program_file <target id> <fasta> <output-directory>\n\n";
				print "\t** Example: $install_dir/bin/$program_file T1006 $install_dir/examples/T1006.fasta $install_dir/test_out/T1006_${method}_hard\n\n";
				print OUT "\t** Example: $install_dir/bin/$program_file T1006 $install_dir/examples/T1006.fasta $install_dir/test_out/T1006_${method}_hard\n\n";
			}
		}
		
		if($method eq 'dncon2')
		{
			$python_env = 1;
			$boost_enable = 1;
		}
	}
	

}
close OUT;
system("chmod +x $install_dir/installation/SAXSDom_test_codes/*sh");

system("cp $install_dir/src/run_SAXSDom.sh $install_dir/bin/run_SAXSDom.sh");
system("chmod +x $install_dir/bin/*.sh");


system("chmod +x $install_dir/src/visualize_SAXSDom_cluster/*.sh");
=cut



sub prompt_yn {
  my ($query) = @_;
  my $answer = prompt("$query (Y/N): ");
  return lc($answer) eq 'y';
}
sub prompt {
  my ($query) = @_; # take a prompt string as argument
  local $| = 1; # activate autoflush to immediately show the prompt
  print $query;
  chomp(my $answer = <STDIN>);
  return $answer;
}


sub configure_file{
	my ($option_list,$prefix) = @_;
	open(IN,$option_list) || die "Failed to open file $option_list\n";
	$file_indx=0;
	while(<IN>)
	{
		$file = $_;
		chomp $file;
		if ($file =~ /^$prefix/)
		{
			$option_default = $install_dir.$file.'.default';
			$option_new = $install_dir.$file;
			$file_indx++;
			print "$file_indx: Configuring $option_new\n";
			if (! -f $option_default)
			{
					die "\nOption file $option_default not exists.\n";
			}	
			
			open(IN1,$option_default) || die "Failed to open file $option_default\n";
			open(OUT1,">$option_new") || die "Failed to open file $option_new\n";
			while(<IN1>)
			{
				$line = $_;
				chomp $line;

				if(index($line,'SOFTWARE_PATH')>=0)
				{
					$line =~ s/SOFTWARE_PATH/$install_dir/g;
					$line =~ s/\/\//\//g;
					print OUT1 $line."\n";
				}else{
					print OUT1 $line."\n";
				}
			}
			close IN1;
			close OUT1;
		}
	}
	close IN;
}


sub configure_tools{
	my ($option_list,$prefix,$DBtool_path) = @_;
	open(IN,$option_list) || die "Failed to open file $option_list\n";
	$file_indx=0;
	while(<IN>)
	{
		$file = $_;
		chomp $file;
		if ($file =~ /^$prefix/)
		{
			$option_default = $DBtool_path.$file.'.default';
			$option_new = $DBtool_path.$file;
			$file_indx++;
			print "$file_indx: Configuring $option_new\n";
			if (! -f $option_default)
			{
					next;
					#die "\nOption file $option_default not exists.\n";
			}	
			
			open(IN1,$option_default) || die "Failed to open file $option_default\n";
			open(OUT1,">$option_new") || die "Failed to open file $option_new\n";
			while(<IN1>)
			{
				$line = $_;
				chomp $line;

				if(index($line,'SOFTWARE_PATH')>=0)
				{
					$line =~ s/SOFTWARE_PATH/$DBtool_path/g;
					$line =~ s/\/\//\//g;
					print OUT1 $line."\n";
				}else{
					print OUT1 $line."\n";
				}
			}
			close IN1;
			close OUT1;
		}
	}
	close IN;
}



sub configure_file2{
	my ($option_list,$prefix) = @_;
	open(IN,$option_list) || die "Failed to open file $option_list\n";
	$file_indx=0;
	while(<IN>)
	{
		$file = $_;
		chomp $file;
		if ($file =~ /^$prefix/)
		{
			@tmparr = split('/',$file);
			$filename = pop @tmparr;
			chomp $filename;
			$filepath = join('/',@tmparr);
			$option_default = $install_dir.$filepath.'/.'.$filename.'.default';
			$option_new = $install_dir.$file;
			$file_indx++;
			print "$file_indx: Configuring $option_new\n";
			if (! -f $option_default)
			{
					die "\nOption file $option_default not exists.\n";
			}	
			
			open(IN1,$option_default) || die "Failed to open file $option_default\n";
			open(OUT1,">$option_new") || die "Failed to open file $option_new\n";
			while(<IN1>)
			{
				$line = $_;
				chomp $line;

				if(index($line,'SOFTWARE_PATH')>=0)
				{
					$line =~ s/SOFTWARE_PATH/$install_dir/g;
					$line =~ s/\/\//\//g;
					print OUT1 $line."\n";
				}else{
					print OUT1 $line."\n";
				}
			}
			close IN1;
			close OUT1;
		}
	}
	close IN;
}


