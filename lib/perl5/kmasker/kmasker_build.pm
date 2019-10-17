package kmasker::kmasker_build;
use Exporter qw(import);
use File::Basename;
use warnings;
use POSIX; 

my $loctime = localtime;
my $loctime2= localtime;
$loctime = strftime('%Y%m%d',localtime); ## outputs 120817


our @ISA = qw(Exporter);
our @EXPORT = qw(
setup_index
build_kindex_jelly
make_config
remove_repository_entry
);
our @EXPORT_OK = qw(build_kindex_jelly remove_kindex set_kindex_external set_private_path set_external_path show_path_infos clean_repository_directory read_config);

## VERSION
my $version_PM_build 	= "1.1.0 rc190712";


sub build_kindex_jelly{	
	my $href_info		= $_[0];
	my $build_config 	= $_[1];
	my $href_repos		= $_[2];
	my $href_path		= $_[3];
	
	my %HASH_info 		= %{$href_info};
	my %HASH_repo		= %{$href_repos};
	my %HASH_path		= %{$href_path};
	my $path_fastqstats	= $HASH_path{"fastq-stats"};
	my $path_jellyfish  = $HASH_path{"jellyfish"};
		
	#LOAD info
	if(defined $build_config){
		#READ config for build
		my $href_this = &read_config($build_config, \%HASH_info, $href_repos, "build");	
		%HASH_info = %{$href_this};	
	}
	
	#INFO
	print "\n .. building kindex for ".$HASH_info{"kindex name"}."\n" if(exists $HASH_info{"kindex name"});
	
	
	#CHECK minimal requirement
	my $go_interactiv;
	$go_interactiv = 1 if(!exists $HASH_info{"genome size"}||!exists $HASH_info{"kindex name"}||!exists $HASH_info{"k-mer"}||!exists $HASH_info{"common name"});
	if(defined $go_interactiv){
		my $href_info_update = &make_config($HASH_info{"PATH_kindex_private"}, $HASH_info{"PATH_kindex_external"}, \%HASH_info);
		%HASH_info = %{$href_info_update};
		$build_config = "repository_".$HASH_info{"kindex name"}.".info";
	}
	else{
		&make_minimal_info_config(\%HASH_info);
		$build_config = "repository_".$HASH_info{"kindex name"}.".info";
	}
	
	#READ info
	my $user_name 		= $HASH_info{"user_name"};
	my $seq				= $HASH_info{"seq"};
	my $k				= $HASH_info{"k-mer"};
	my $short_tag		= $HASH_info{"kindex name"};
	my $path			= $HASH_info{"path_bin"};
	my $threads			= $HASH_info{"threads"};
	my $size			= $HASH_info{"size"};
	my $seq_type		= $HASH_info{"sequence type"};
	
	#LOAD expert setting for build
	my $setting = "-s ".$size."G -t ".$threads;	#server setting
	if(exists $HASH_info{"user setting jelly"}){
		$parameter_extern = $HASH_info{"user setting jelly"};
		my @ARRAY_config = split(";", $parameter_extern);
		for(my $s=0;$s<scalar(@ARRAY_config);$s++){
			my $set = $ARRAY_config[$s];
			$set =~ s/^ //;							# no empty space at beginning
			$set .= "--".$set if($set !~ /^\-\-/) ; 	# add parameter prefix structure
			if($set =~ /size/){
				$set .= "G" if(($set !~ /M$/)&&($set !~ /G$/));
			}
			$ARRAY_config[$s] = $set;
		}
		$setting = join (" ", @ARRAY_config);
		print "\n using modified user settings : ".$setting."\n";
	}
	
	#create command shell script for background
	my @ARRAY_input = split(" ", $seq);
	my @ARRAY_shell = ();
	foreach my $call (@ARRAY_input){
		push(@ARRAY_shell, $path_fastqstats." ".$call." >".$call.".kmstats ");
	}
	
	#BUILD JELLY index	
	#if($seq_type eq "reads"){
	#	$setting .= " -C"
	#}
	
	print "\n ... start construction of kindex with the following parameters ".$setting." \n";
	my $FILE_jelly = "KINDEX_".$HASH_info{"kindex name"}.".jf";
	
	#CHECK suffix 
	my ($ext) = $seq =~ /(\.[^.]+)$/;
	if($ext =~ /gz$/){
		push(@ARRAY_shell, "zcat ".$seq." | $path_jellyfish count -m ".$k." ".$setting." -o ".$FILE_jelly);
	}else{
		push(@ARRAY_shell, "$path_jellyfish count -m ".$k." ".$setting." -o ".$FILE_jelly." ".$seq);
	}
		
	#CREATE BUILD script
	#system("cp ".$path."/.kmasker_background_process .");
	
	my $BASH 	= new IO::File(".kmasker_background_process", '>') or die "could not write shell script: $!\n";
	print $BASH "#!/bin/sh\n";
	foreach my $call (@ARRAY_shell){
		print $BASH "{ ".$call." ;} &\n";
	}
	print $BASH "wait\n";
	$BASH->close();
	chmod(0775, ".kmasker_background_process");
	system("./.kmasker_background_process");	
		
	system("cp ./.kmasker_background_process build_command.sh") if(defined $HASH_info{"verbose"});
	print "\n ... finished kindex construction!\n";
	
	
	#UPDATE repository information	
	$HASH_info{"absolut_path"} 			= getcwd if(!exists $HASH_info{"absolut_path"});	
	$HASH_info{"version BUILD"}			= $version_PM_build;
	$HASH_info{"sequencing depth"} 		= &read_stats($HASH_info{"genome size"});
	$HASH_info{"file"}					= $FILE_jelly;
	$HASH_info{"created"}				= $loctime2;
	&update_repository_information(\%HASH_info);	
	
	#STORE constructed kindex and repository.info
	my $PATH_final = $HASH_info{"PATH_kindex_private"};
	$PATH_final = $HASH_info{"PATH_kindex_external"} if($HASH_info{"status"} eq "external");
	if(!(-d $PATH_final)){
		system("mkdir ".$PATH_final);
	}
	if(-e $FILE_jelly){
		system("mv $FILE_jelly ".$PATH_final);
		system("mv $build_config ".$PATH_final);
		
		print "\n --> KINDEX ".$HASH_info{"kindex name"}." was sucessfully constructed!";
		print "\n --> [path] = ".$PATH_final.$FILE_jelly."\n\n";
	}
	$HASH_info{"absolut_path"} = $PATH_final.$FILE_jelly;
}


## subroutine
#
sub update_repository_information{
	#subroutine aims to fill repository information after build process is completed.
	#
	my $href 		= $_[0];
	my %HASH_info 	= %{$href};
	
	my $REPO 		= new IO::File("repository_".$HASH_info{"kindex name"}.".info", "r") or die "could not read file for repository_".$HASH_info{"kindex name"}.".info file: $!\n";
	my $REPO_update = new IO::File("repository_".$HASH_info{"kindex name"}.".info_update", "w") or die "could not write file for repository_".$HASH_info{"kindex name"}.".info_update file: $!\n";
	
	#1
	#READ available info from repository.info into HASH_info (ONLY if HASH is empty)
	while(<$REPO>){
		my $line 	= $_;
		next if($line =~ /^$/);
		$line		=~ s/\n$//;
		my $key	 	= +(split(":", $line))[0];
		my $value 	= +(split("\t", $line))[1];
		$value		= "" if(!defined $value);
			
		if($key =~ /^kindex name/)		{ 	$HASH_info{"kindex name"} 		= $value if(($HASH_info{"kindex name"} eq "")		||(!exists $HASH_info{"kindex name"})); 	}
		if($key =~ /^status/)     		{	$HASH_info{"status"}    		= $value if(($HASH_info{"status"} eq "")			||(!exists $HASH_info{"status"})); 	}
		if($key =~ /^common name/)		{	$HASH_info{"common name"}   	= $value if(($HASH_info{"common name"} eq "")		||(!exists $HASH_info{"common name"})); 	}
		if($key =~ /^scientific name/)	{	$HASH_info{"scientific name"}  	= $value if(($HASH_info{"scientific name"} eq "")	||(!exists $HASH_info{"scientific name"})); 	}
		if($key =~ /^type/)				{	$HASH_info{"type"}  			= $value if(($HASH_info{"type"} eq "")				||(!exists $HASH_info{"type"})); 	}
		if($key =~ /^genome size/)		{	$HASH_info{"genome size"}  		= $value if(($HASH_info{"genome size"} eq "")		||(!exists $HASH_info{"genome size"})); 	}
		if($key =~ /^sequencing depth/)	{	$HASH_info{"sequencing depth"}  = $value if(($HASH_info{"sequencing depth"} eq "")	||(!exists $HASH_info{"sequencing depth"})); 	}
		if($key =~ /^k-mer/)			{	$HASH_info{"k-mer"}				= $value if(($HASH_info{"k-mer"} eq "")				||(!exists $HASH_info{"k-mer"})); 	}
		if($key =~ /^sequence type/)	{	$HASH_info{"sequence type"}  	= $value if(($HASH_info{"sequence type"} eq "")		||(!exists $HASH_info{"sequence type"})); 	}
		if($key =~ /^general notes/)	{	$HASH_info{"general notes"}  	= $value if(($HASH_info{"general notes"} eq "")		||(!exists $HASH_info{"general notes"})); 	}
		if($key =~ /^file/)				{	$HASH_info{"file"}  			= $value if(($HASH_info{"file"} eq "")				||(!exists $HASH_info{"file"})); 	}
		if($key =~ /^created/)			{	$HASH_info{"created"}  			= $value if(($HASH_info{"created"} eq "")			||(!exists $HASH_info{"created"})); 	}
		if($key =~ /^version KMASKER/)	{	$HASH_info{"version KMASKER"}  	= $value if(($HASH_info{"version KMASKER"} eq "")	||(!exists $HASH_info{"version KMASKER"})); 	}
		if($key =~ /^version BUILD/)	{	$HASH_info{"version BUILD"}  	= $value if(($HASH_info{"version BUILD"} eq "")		||(!exists $HASH_info{"version BUILD"})); 	}
			
	}
	
	#2
	#WRITE FULL HASH into repository.info_update
	print $REPO_update sprintf("%-15s %s", "kindex name :", " ")."\t". 		$HASH_info{"kindex name"} ."\n";
	print $REPO_update sprintf("%-15s %s", "status :", " ")."\t".			$HASH_info{"status"}."\n";
	print $REPO_update sprintf("%-15s %s", "common name :", " ")."\t".		$HASH_info{"common name"}."\n";
	print $REPO_update sprintf("%-15s %s", "scientific name :", " ")."\t".	$HASH_info{"scientific name"}."\n";
	print $REPO_update sprintf("%-15s %s", "type :", " ")."\t".				$HASH_info{"type"}."\n";
	print $REPO_update sprintf("%-15s %s", "genome size :", " ")."\t".		$HASH_info{"genome size"} ."\n";
	print $REPO_update sprintf("%-15s %s", "sequencing depth :", " ")."\t".	$HASH_info{"sequencing depth"}."\n";
	print $REPO_update sprintf("%-15s %s", "k-mer :", " ")."\t".			$HASH_info{"k-mer"}	."\n";
	print $REPO_update sprintf("%-15s %s", "sequence type :", " ")."\t".	$HASH_info{"sequence type"}."\n";
	print $REPO_update sprintf("%-15s %s", "general notes :", " ")."\t".	$HASH_info{"general notes"}."\n";
	print $REPO_update sprintf("%-15s %s", "file :", " ")."\t".				$HASH_info{"file"}."\n";
	print $REPO_update sprintf("%-15s %s", "created :", " ")."\t".			$HASH_info{"created"}."\n";
	print $REPO_update sprintf("%-15s %s", "version KMASKER :", " ")."\t".	$HASH_info{"version KMASKER"}."\n";
	print $REPO_update sprintf("%-15s %s", "version BUILD :", " ")."\t".	$HASH_info{"version BUILD"}."\n";
	
	$REPO->close();
	$REPO_update->close();
	
	system("mv repository_".$HASH_info{"kindex name"}.".info_update repository_".$HASH_info{"kindex name"}.".info");	
}

## subroutine
#
sub make_config(){
	
	my $PATH_kindex_private = $_[0];
	my $PATH_kindex_external= $_[1];
	my $href				= $_[2];
	my %HASH_info			= %{$href};
		
		#Start interactiv mode
		
		print "\n Setup of configuration file\n";
		print   " ---------------------------\n";
		print "\n To create the KINDEX please provide:\n";
		print "\n Written in '[]' you find the default values. Type <enter> to skip!\n\n";
		
		#1
		my $default = "";
		$HASH_info{"kindex name"} = &promptUser("\n\tName of k-mer index (e.g. 'HvMRX' for Hordeum vulgare cultivar Morex)\n\t", $loctime);
		#check if name of index is unique (private and global)
		my $unique = 1;
		
		opendir( my $DIR_P, $PATH_kindex_private );
		while ( my $entry = readdir $DIR_P ) {
			$entry =~ s/KINDEX_//;
			$entry =~ s/\.jf$//;
			$entry = lc($entry);
			$unique = 0 if($entry eq lc($HASH_info{"kindex name"}));
		}
		close $DIR_P;
		
		opendir( my $DIR_E, $PATH_kindex_external );
		while ( my $entry = readdir $DIR_E ) {
			$entry =~ s/KINDEX_//;
			$entry =~ s/\.jf$//;
			$entry = lc($entry);
			$unique = 0 if($entry eq lc($HASH_info{"kindex name"}));
		}
		close $DIR_E;
		
		if($unique eq 0){
			print "\n name of index is not unique. Kmasker is stopped!\n\n";
			exit();
		}
		
		#2
		#$HASH_info{"status"} = &promptUser("\n\tStatus of index (private or external)\n\t", "private");
		$HASH_info{"status"} = "private";
			
		#3
		$default = "unknown";
		$default = $HASH_info{"common name"} if(exists $HASH_info{"common name"});
		$HASH_info{"common name"} = &promptUser("\n\tCommon name (e.g. 'barley')\n\t", $default);
		if($HASH_info{"common name"} eq ""){
			print "\n [WARNING]";
			print "\n Required value of 'common name ' is missing. Kmasker is stopped!\n\n";
			exit();
		}
		
		#4
		$HASH_info{"scientific name"} = &promptUser("\n\tScientific_name (e.g. 'Hordeum vulgare')\n\t", "");
		
		#5
		$HASH_info{"type"} = &promptUser("\n\tType (cultivar; genotype; etc. e.g. 'Morex'))\n\t", "");
		
		#6
		$default = "";
		$default = $HASH_info{"genome size"} if(exists $HASH_info{"genome size"});
		$HASH_info{"genome size"} = &promptUser("\n\tGenome size in Mb (e.g. '5100' for barley)\n\t", $default);
		if(($HASH_info{"genome size"} eq "")||($HASH_info{"genome size"} !~ (/^\d+$/))){
			print "\n [WARNING]";
			print "\n Required value of 'genome size' is missing or is not a number. Kmasker was stopped!\n\n";
			exit();
		}
		
		#7
		$default = 21;
		$default = $HASH_info{"k-mer"} if(exists $HASH_info{"k-mer"});
		$HASH_info{"k-mer"} = &promptUser("\n\tk-mer length \n\t", $default);
		
		#8
		$HASH_info{"sequence type"} = &promptUser("\n\tSequence_type (e.g. WGS, RNAseq, assembly)\n\t", "WGS");
		
		#9
		$HASH_info{"general notes"} = &promptUser("\n\tGeneral notes\n\t", "");
		
		#CREATE configuratio file
		my $USER_kindex_info 	= new IO::File("repository_".$HASH_info{"kindex name"}.".info", "w") or die "could not write file for repository_".$HASH_info{"kindex name"}.".info file: $!\n";
		print $USER_kindex_info sprintf("%-15s %s", "kindex name :", " ")."\t".$HASH_info{"kindex name"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "status :", " ")."\t".$HASH_info{"status"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "common name :", " ")."\t".$HASH_info{"common name"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "scientific name :", " ")."\t".$HASH_info{"scientific name"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "type :", " ")."\t".$HASH_info{"type"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "genome size :", " ")."\t".$HASH_info{"genome size"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "sequencing depth :", " ")."\t\n";
		print $USER_kindex_info sprintf("%-15s %s", "k-mer :", " ")."\t".$HASH_info{"k-mer"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "sequence type :", " ")."\t".$HASH_info{"sequence type"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "general notes :", " ")."\t".$HASH_info{"general notes"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "file :", " ")."\t\n";
		print $USER_kindex_info sprintf("%-15s %s", "created :", " ")."\t".$loctime2."\n";
		print $USER_kindex_info sprintf("%-15s %s", "version KMASKER :", " ")."\t".$HASH_info{"version KMASKER"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "version BUILD :", " ")."\t".$version_PM_build."\n";
				
		print $USER_kindex_info "\n";	
		
		$USER_kindex_info->close();
		
		#START
		my $start;
		$start = &promptUser("\n\n Would you like to start building the KINDEX '".$HASH_info{"kindex name"}." (type 'no' for exit)'\n\t", "yes");
		print "\n\n";
		
		#EXIT
		if((lc($start) ne "yes")&&(lc($start) ne "y")){
			exit();
		}

		return \%HASH_info;
}


## subroutine
#
sub make_minimal_info_config(){
	my $href 				= $_[0];
	
	my %HASH_info			= %{$href};
	my $USER_kindex_info 	= new IO::File("repository_".$HASH_info{"kindex name"}.".info", "w") or die "could not write file for repository_".$HASH_info{"kindex name"}.".info file: $!\n";		
	
	print "\n Creating configuration file\n";		
	#CREATE configuratio file
		print $USER_kindex_info sprintf("%-15s %s", "kindex name :", " ")."\t".$HASH_info{"kindex name"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "status :", " ")."\t".$HASH_info{"status"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "common name :", " ")."\t".$HASH_info{"common name"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "scientific name :", " ")."\t".$HASH_info{"scientific name"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "type :", " ")."\t".$HASH_info{"type"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "genome size :", " ")."\t".$HASH_info{"genome size"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "sequencing depth :", " ")."\t\n";
		print $USER_kindex_info sprintf("%-15s %s", "k-mer :", " ")."\t".$HASH_info{"k-mer"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "sequence type :", " ")."\t".$HASH_info{"sequence type"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "general notes :", " ")."\t".$HASH_info{"general notes"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "file :", " ")."\t\n";
		print $USER_kindex_info sprintf("%-15s %s", "created :", " ")."\t".$loctime."\n";
		print $USER_kindex_info sprintf("%-15s %s", "version KMASKER :", " ")."\t".$HASH_info{"version KMASKER"}."\n";
		print $USER_kindex_info sprintf("%-15s %s", "version BUILD :", " ")."\t".$version_PM_build."\n";
		print $USER_kindex_info "\n";	
	$USER_kindex_info->close();

}

## subroutine
#
sub read_config(){
		my $build_config 	= $_[0];
		my $href_info		= $_[1];
		my $href_repos_this = $_[2];
		my $call_process	= $_[3];
		my %HASH_info_this 	= %{$href_info};
		my %HASH_repo_this	= %{$href_repos_this};	
		
		#HASH
		my %HASH_code_words = ();
		$HASH_code_words{"kindex name"}		= 1;
		$HASH_code_words{"common name"}		= 1;
		$HASH_code_words{"scientific name"}	= 1;
		$HASH_code_words{"genome size"}		= 1;
		$HASH_code_words{"type"}			= 1;
		$HASH_code_words{"sequencing depth"}= 1;
		$HASH_code_words{"k-mer"}			= 1;
		$HASH_code_words{"sequence type"}	= 1;
		$HASH_code_words{"general notes"}	= 1;
		$HASH_code_words{"status"}			= 1;
		
		my $INPUT_kindex_info 	= new IO::File($build_config, "r") or die "could not read file for repository information (parameter '--config'): $!\n";
		while(<$INPUT_kindex_info>){
			my $line 	= $_;
			next if($line =~ /^$/);
			$line =~ s/\n$//;
			my $key	 	= +(split(":", $line))[0];
			$key		=~ s/\s+$//;
			my $value 	= +(split("\t", $line))[1];
			$value = "" if(!defined $value);
			if(exists $HASH_code_words{$key}){
				$HASH_info_this{$key} = $value if((!exists $HASH_info_this{$key}) || ($HASH_info_this{$key} eq ""));	#mean that input from command line has priority (config is second)
			}
		}		
		
		if($call_process eq "build"){
			#CHECK if mandatory field are set correctly
			my $check = 1;
			
			$check = 0 if((!exists $HASH_info_this{"kindex name"}) 	|| ($HASH_info_this{"kindex name"} eq ""));
			$check = 0 if((!exists $HASH_info_this{"k-mer"}) 		|| ($HASH_info_this{"k-mer"} eq ""));
			$check = 0 if((!exists $HASH_info_this{"common name"}) 	|| ($HASH_info_this{"common name"} eq ""));
			$check = 0 if((!exists $HASH_info_this{"genome size"})	|| ($HASH_info_this{"genome size"} eq ""));
			if($check == 0){
				#STOP
				print "\n\n WARNING: Kmasker (build) was stopped!!!\
					     \n Missing information in configuration!\n\n";
				exit(0);
			}
			
			if($HASH_info_this{"sequence type"} ne "assembly"){
				$HASH_info_this{"sequence type"} = "reads";
			}
			
			#CEHCK for extisting entry
			my $short_tag = $HASH_info_this{"kindex name"};
			
			if(exists $HASH_repo_this{$short_tag}){
				#kindex with identical name exists, kmasker needs to be stopped to avoid overwriting!
				print "\n\n WARNING: Kmasker (build) was stopped!!!\
					     \n A kindex with identical name already exists in repository\n\n";
				exit(0);
			}
		}		
		
		#RETURN
		return \%HASH_info_this;			
}


## subroutine
#
sub read_stats(){
	#INPUT
	my $gs 			= $_[0];
	#VAR
	my $calculation = 1;
	my $total_bases	= 0;
	
	opendir(Dir, ".") or die "cannot open directory .";
	@STATS = grep(/\.kmstats$/,readdir(Dir));
	foreach my $file (@STATS) {
		my $INPUT_stats = new IO::File($file, "r") or die "could not read ".$file." $!\n";
		while(<$INPUT_stats>){
			my $line = $_;
			next if($line =~ /^$/);
			$line =~ s/\n$//;
			if($line =~ /^total bases/){
				my @ARRAY_tmp 	= split("\t", $line);
				$total_bases += $ARRAY_tmp[1];
			}
		}
		$INPUT_stats->close();
	}
	system("rm *.kmstats");
	
	$calculation 	= sprintf("%.1f", $total_bases / ($gs * 1000000));
	if($calculation < 1){
		print "\n Notification:";
		print "\n The calculated sequencing depth of your dataset is below 1-fold!";
		print "\n This is expected if you used an assembled dataset as input for the index construction,";
		print "\n but might be too low if you used WGS data as input. It is possible to detect";
		print "\n abundant sequences but its not recommended to use this dataset for detection";
		print "\n of low-coverage regions. The normalisation factor of the constrcuted index is set to 1x\n";
		$calculation 	= 1;		
	}	
	return $calculation;
}

## subroutine
#
sub remove_kindex(){	
	my $kindex_shorttag			= $_[0];
	my $href_this 				= $_[1];
	my %HASH_repository_kindex	= %{$href_this};	
	
	if(exists $HASH_repository_kindex{$kindex_shorttag}){
		my @ARRAY_kindex_overview	= split("\t", $HASH_repository_kindex{$kindex_shorttag});
		my $status 					= $ARRAY_kindex_overview[3];
		my $abs_path_kindex_dir		= $ARRAY_kindex_overview[4];
		
		if(($status eq "global")||($status eq "external")){
			print "\n WARNING: the kindex ".$kindex_shorttag." is external or global (not permitted to remove)!\n\n";
		}else{
			#REMOVE DATA
			print "\n REMOVING the index structure for KINDEX '$kindex_shorttag' from :\n ".$abs_path_kindex_dir."\n\n";
			system("rm -r ".$abs_path_kindex_dir."*".$kindex_shorttag."*");
		}			
	}else{
		print "\n The kindex ".$kindex_shorttag." does not exist! Nothing to remove!\n\n";
	}
}


## subroutine
#
sub set_kindex_external(){
	
	my $kindex_shorttag		= $_[0];
	my $href_this 			= $_[1];
	my $path_private		= $_[2];
	my $path_global			= $_[3];
	
	my %HASH_repository_kindex_this		= %{$href_this};	
	my @ARRAY_repository				= split("\t", $HASH_repository_kindex_this{$kindex_shorttag});
	my $path 							= $ARRAY_repository[2];
	
	#MOVE directory into global path
	system("mv ".$path_private."KINDEX_".$kindex_shorttag." ".$path_global);
	
	#READ & EDIT repository.info
	my $utmp = "tmp_".$loctime.".txt";
	my $REPO_private 		= new IO::File($path_global."repository.info", "r") or die "\n unable to read repository.info list $!";
	my $REPO_private_edit 	= new IO::File($utmp, "w") or die "\n unable to write tmp repository.info $!";
	
	while(<$REPO_private>){
		if($_ =~ /^status/){
			$_ =~ s/private/global/;
			print $REPO_private_edit $_;
		}else{
			print $REPO_private_edit $_;
		}
	}
	$REPO_private->close();
	$REPO_private_edit->close();
	
	#RENAME file
	system("mv ".$utmp." ".$path_global."repository.info");			
}


## subroutine
#
sub set_kindex_name(){
	#RENAME existing KINDEX	
	
	#rename in folder
	
	#rename in file
	
	#rename in repository.info
}


## subroutine
#
sub clean_repository_directory(){
	#routinly check if KINDEX datasets exist, that are not complete
	# e.g. repository.file is missing or is empty
	
	#general cleaning
	system("rm \.kmasker_background_process") if(-e "\.kmasker_background_process");
	system("rm repository.info") if(-e "repository.info");
}

## subroutine
#
sub set_private_path(){
	my $private_path 	= $_[0];
	
	my $old_private_path="";
	$private_path 		.= "/" if($private_path !~ /\/$/ );
	my $uconf 			= $ENV{"HOME"}."/.kmasker_user.config";
	my $USER_CONF_OLD 	= new IO::File($uconf, "r") or die "could not read old user conf : $!\n";
	my $USER_CONF 		= new IO::File($uconf.".tmp", "w") or die "could not write new user conf : $!\n";
	while(<$USER_CONF_OLD>){
		my $line = $_;
		$line =~ s/\n//;
		if($line =~ /^PATH_kindex_private/){
			print $USER_CONF "PATH_kindex_private=".$private_path."\n";
			$old_private_path = +(split("=", $line))[1];
		}else{
			print $USER_CONF $line."\n";
		}		
	}
	system("mv ".$uconf.".tmp ".$uconf);
	$USER_CONF->close();
	$USER_CONF_OLD->close();
	
}

## subroutine
#  
sub set_external_path(){
	#INPUT
	my $new_path		= $_[0];
	my $href_info		= $_[1];
	
	if($new_path =~ /^\./){
		print "\n .. please use absolut path! Kmasker was stopped!\n";
		exit();
	}
	
	my $old_path		="";
	my $conf 			= $ENV{"HOME"}."/.kmasker_user.config";
	
	my $CONF_OLD 	= new IO::File($conf, "r") or die "could not read old conf : $!\n";
	my $CONF 		= new IO::File($conf.".tmp", "w") or die "could not write new conf : $!\n";

	# It is not allowed to change the global path by any user (only root)
	# Therefor, the set_global_path function only changes the external directory called 'PATH_kindex_external'	
	while(<$CONF_OLD>){
		my $line = $_;
		$line =~ s/\n//;
		my $pattern = "PATH_kindex_external";
		if($line =~ /^($pattern)/){
			print $CONF $pattern."=".$new_path."\n";
			$old_path = +(split("=", $line))[1];
			$old_path = "\'\'" if(!defined $old_path);
			print "\n OLD path: ".$old_path;
			print "\n NEW path: ".$new_path."\n\n";
		}else{
			print $CONF $line."\n";
		}		
	}
	system("mv ".$conf.".tmp ".$conf);
	$CONF->close();
	$CONF_OLD->close();
	
}


## subroutine
#
sub show_path_infos(){
	my $uconf 			= $ENV{"HOME"}."/.kmasker_user.config";
	my $USER_CONF_OLD 	= new IO::File($uconf, "r") or die "could not read old user conf : $!\n";
	print "\nKmasker path info:\n";
	while(<$USER_CONF_OLD>){
		my $line = $_;
		$line =~ s/\n//;
		if($line =~ /^PATH_kindex_private/){
			print "\n".$line;
		}		
		if($line =~ /^PATH_kindex_external/){
			print "\n".$line;
		}		
	}
	print "\n\n";
}


## subroutine
#
sub move_kindex_structures(){
	my $PATH_kindex_NEW = $_[0];
	my $PATH_kindex_OLD = $_[1];
	my $href_repository = $_[2];
	
	#quit if same
	exit if($PATH_kindex_NEW eq $PATH_kindex_OLD);
	
	my %HASH_repository_kindex_this = %{$href_repository};	
	#integrate into HASH
	#$HASH_repository_kindex{$kindex_name} = $kindex_name."\t".$common_name."\t".$status."\t".$kindex_path;
	
	foreach my $kindex_name (keys %HASH_repository_kindex_this){
		my @ARRAY_info = split("\t", $HASH_repository_kindex_this{$kindex_name});
		#ALL KINDEX that start with $PATH_kindex_OLD
		if($ARRAY_info[3] =~ /^$PATH_kindex_OLD/){
			print "This is detected to be moved to new location";
		}
	}
}


## subroutine
#
sub promptUser {

#-------------------------------------------------------------------------#
# promptUser, a Perl subroutine to prompt a user for input.
# Copyright 2010 Alvin Alexander, http://www.devdaily.com
# This code is shared here under the 
# Creative Commons Attribution-ShareAlike Unported 3.0 license.
# See http://creativecommons.org/licenses/by-sa/3.0/ for more information.
#-------------------------------------------------------------------------#
#LINK: https://alvinalexander.com/perl/edu/articles/pl010005

   #-------------------------------------------------------------------#
   #  two possible input arguments - $promptString, and $defaultValue  #
   #  make the input arguments local variables.                        #
   #-------------------------------------------------------------------#

   local($promptString,$defaultValue) = @_;

   #-------------------------------------------------------------------#
   #  if there is a default value, use the first print statement; if   #
   #  no default is provided, print the second string.                 #
   #-------------------------------------------------------------------#

   if ($defaultValue) {
      print $promptString, "[", $defaultValue, "]: ";
   } else {
      print $promptString, ": ";
   }

   $| = 1;               # force a flush after our print
   $_ = <STDIN>;         # get the input from STDIN (presumably the keyboard)


   #------------------------------------------------------------------#
   # remove the newline character from the end of the input the user  #
   # gave us.                                                         #
   #------------------------------------------------------------------#

   chomp;

   #-----------------------------------------------------------------#
   #  if we had a $default value, and the user gave us input, then   #
   #  return the input; if we had a default, and they gave us no     #
   #  no input, return the $defaultValue.                            #
   #                                                                 # 
   #  if we did not have a default value, then just return whatever  #
   #  the user gave us.  if they just hit the <enter> key,           #
   #  the calling routine will have to deal with that.               #
   #-----------------------------------------------------------------#

   if ("$defaultValue") {
      return $_ ? $_ : $defaultValue;    # return $_ if it has a value
   } else {
      return $_;
   }
}




1;
