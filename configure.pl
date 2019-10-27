#!/usr/bin/perl -w
 use FileHandle; # use FileHandles instead of open(),close()
 use Cwd;
 use Cwd 'abs_path';

######################## !!! customize settings here !!! ############################
#																					#
# Set directory of SAXSDom databases and tools								        #

$SAXSDom_db_tools_dir = "/data/commons/SAXSDom_db_tools/";						        

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

$tools_dir = "$SAXSDom_db_tools_dir/tools";

if(!(-d $tools_dir))
{
	die "Failed to find tools under $SAXSDom_db_tools_dir/\n";
}

if($SAXSDom_db_tools_dir eq "$cur_dir/")
{
	die "Same directory as SAXSDom main folder. Differnt path for original databases/tools folder $SAXSDom_db_tools_dir is recommended.\n";
}
#create link for databases and tools
if(-d "${install_dir}tools")
{
  `rm ${install_dir}tools`; 
}
`ln -s $tools_dir ${install_dir}tools`;
=pod

if (prompt_yn("SAXSDom will be installed into <$install_dir> ")){

}else{
	die "The installation is cancelled!\n";
}
print "Start install SAXSDom into <$install_dir>\n"; 
=cut



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

if(-e "Mocapy++-1.07.tar.gz")
{
  `tar -zxf Mocapy++-1.07.tar.gz`;
  `cp -ar  Mocapy_src/* Mocapy++-1.07/`;
  `cp -ar  $install_dir/src/SAXSDom.cpp Mocapy++-1.07/examples`;
  
  open(OUT,">$install_dir/compile_SAXSDom.sh") || die "Failed to open file $install_dir/compile_SAXSDom\n";
  print OUT "#!/bin/bash -e\n\n";
  print OUT "echo \" Start compile SAXSDom (will take ~10 min)\"\n\n";
  print OUT "cd $install_dir/installation/Mocapy++-1.07\n\n";
  print OUT "export LD_LIBRARY_PATH=$install_dir/tools/boost_1_55_0/lib:\$LD_LIBRARY_PATH\n\n";
  print OUT "$install_dir/tools/cmake-3.5.2/bin/cmake -DBOOST_ROOT='$install_dir/tools/boost_1_55_0/' -DLAPACK_LIBRARY:FILEPATH='$install_dir/tools/lapack-3.4.1/liblapack.a' .\n\n";
  print OUT "make\n\n";
  close OUT;  
}else{
  die "Failed to find Mocapy++-1.07.tar.gz, contact us\n";
}

chdir($install_dir);
=pod

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
=cut

system("mv $install_dir/installation/SAXSDom_test_codes/T0_run_SAXSDom*.sh $install_dir/examples");
system("chmod +x $install_dir/examples/*.sh");



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


