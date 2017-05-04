package kmasker::kmasker_build;
use Exporter qw(import);
use File::Basename;
#use strict;
use warnings;
use Cwd;

our @ISA = qw(Exporter);
our @EXPORT = qw(
setup_index
build_kindex_jelly
make_config
remove_repository_entry
);
our @EXPORT_OK = qw(build_kindex_jelly make_config remove_repository_entry set_kindex_global clean_repository_directory set_private_path);

## VERSION
my $version_PM_build 	= "0.0.3 rc170504";

sub build_kindex_jelly{	
	my $href_info		= $_[0];
	my $build_config 	= $_[1];
	my $href_repos		= $_[2];
	my $href_path		= $_[3];
	my $store_input		= $_[4];	
	my %HASH_info 		= %{$href_info};
	my %HASH_repo		= %{$href_repos};
	my %HASH_path		= %{$href_path};
	my $path_fastqstats	= $HASH_path{"fastq-stats"};
		
	#LOAD info
	if(defined $build_config){
		#READ config for build
		my $href_this = &read_config($build_config, \%HASH_info, $href_repos);	
		%HASH_info = %{$href_this};	
	}
	
	#READ info
	my $user_name 		= $HASH_info{"user_name"};
	my $seq				= $HASH_info{"seq"};
	my $k				= $HASH_info{"k-mer"};
	
	
	#LOAD expert setting for build
	my $parameter_extern = $HASH_info{"expert_setting"};
	#my $setting = "-s 2G -t 1";	#notebook setting
	my $setting = "-s 24G -t 10";	#server setting
	if($parameter_extern ne ""){
		$setting = $parameter_extern;
	}
	
	#Create single input file and calculate md5sum for repository
	system("cat ".$seq." >INPUT.fastq");
	
	#make stats
	if($HASH_info{"genome_size"} !~ /^N/){
		#if genome size is provided that sequenincg depth is calculated automatically
		system($path_fastqstats." INPUT.fastq >INPUT.fastq.stats &");
	}
	
	#create meta data
	my @ARRAY_help 		= split(/\./, $seq);
	my $end 			= pop(@ARRAY_help);
	$end				=~ s/ //g;
	my $RESULT_md5_name = `md5sum INPUT.fastq`;	
	$RESULT_md5_name 	=~ s/\n//g;
	$RESULT_md5_name 	=~ s/  /\t/g;
	$RESULT_md5_name 	=~ s/ //g;
	my @ARRAY_tmp 		= split("\t", $RESULT_md5_name);
	my $md5sum 			= $ARRAY_tmp[0];
	$HASH_info{"seq"}	= "INPUT_".$md5sum.".".$end;
	$HASH_info{"md5sum"}= $md5sum;
	$HASH_info{"status"}= "private";
	$seq				= "INPUT_".$md5sum.".".$end;
	system("mv INPUT.fastq INPUT_".$md5sum.".".$end);
		
	
	#BUILD JELLY index	
	print "\n ... start construction of kindex with the following parameters ".$setting." \n"; 
	system("jellyfish count -m ".$k." ".$setting." -o KINDEX_".$HASH_info{"short_tag"}."_".$md5sum."_k".$k.".jf ".$seq."");
	print "\n ... finished kindex construction!\n";
	
	#FILL with default if information is not provided	
	$HASH_info{"absolut_path"} = getcwd if(!exists $HASH_info{"absolut_path"});	
	#
	if($HASH_info{"genome_size"} !~ /^N/){
		#calculate sequening depth
		my $caculated_sequencing_depth 	= 1;
		$caculated_sequencing_depth 	= &read_stats($HASH_info{"genome_size"});
		$HASH_info{"sequencing_depth"} 	= $caculated_sequencing_depth;
	}
	$HASH_info{"sequencing_depth"} = 1 if((!exists $HASH_info{"sequencing_depth"})|| ($HASH_info{"sequencing_depth"} < 1));
	
	
	#MOVE: make folder in directory and move jelly index
	system("mkdir ".$HASH_info{"PATH_kindex_private"}."KINDEX_".$HASH_info{"short_tag"});
	if(-e "KINDEX_".$HASH_info{"short_tag"}."_".$md5sum."_k".$k.".jf"){
		system("mv KINDEX_".$HASH_info{"short_tag"}."_".$md5sum."_k".$k.".jf ".$HASH_info{"PATH_kindex_private"}."KINDEX_".$HASH_info{"short_tag"});
		system("cp ".$build_config." ".$HASH_info{"PATH_kindex_private"}."KINDEX_".$HASH_info{"short_tag"});
	}
	$HASH_info{"absolut_path"} = $HASH_info{"PATH_kindex_private"}."KINDEX_".$HASH_info{"short_tag"}."/";
	
	#SAVE info of build
	&write_repository_entry(\%HASH_info);
	
	#STORE input sequence as gzip -9 compressed files???
	if(defined $store_input){
		system("gzip -9 ".$seq);
		system("mv ".$seq.".gz ".$HASH_info{"PATH_kindex_private"}."KINDEX_".$HASH_info{"short_tag"});		
	}else{
		system("rm ".$seq);
	}
}

## subroutine
#
sub setup_index{
	my ($file, $k) = @_;
	system("head -n 100 ".$file." >input_head.txt");
	system("tail -n 100 ".$file." >input_tail.txt");
	system("cat input_head.txt input_tail.txt >input.txt");	
	my $md5_name = `md5sum input.txt`;
	system("rm input_head.txt input_tail.txt input.txt");
	
	#create index
	&make_kindex_jelly($file, $md5_name, $k);	
}


## subroutine
#
sub extend_repository_information{
	#subroutine aims to fill repository information (e.g missing md5sum).
	#read private repository of user and check if md5sum information is missing. 
	#If yes, fill this gap by using the fastq or fastq files used to construct the index.
}

## subroutine
#
sub make_config(){
		my $USER_kindex_info 	= new IO::File("build_kindex.config", "w") or die "could not write file for repository information: $!\n";
		
		print $USER_kindex_info "##Kmasker build config file\n";
		#1
		print $USER_kindex_info "#short_tag [mandatory]\n\n";
		#2
		print $USER_kindex_info "#common_name [obligarory]\n\n";
		#3
		print $USER_kindex_info "#scientific_name [obligarory]\n\n";
		#4
		print $USER_kindex_info "#genome_size [mandatory]\n\n";
		#5
		print $USER_kindex_info "#type (cultivar;genotype;etc.) [obligarory]\n\n";
		#6
		print $USER_kindex_info "#sequencing_depth [mandatory]\n\n";
		#7
		print $USER_kindex_info "#k-mer [obligatory]\n\n";
		#8
		print $USER_kindex_info "#sequence_type (reads or assembly) [obligatory]\n\n";
		#9
		print $USER_kindex_info "#note [obligarory]\n\n";

		print $USER_kindex_info "\n";		
}


## subroutine
#
sub read_config(){
		my $build_config 	= $_[0];
		my $href_info		= $_[1];
		my $href_repos_this = $_[2];
		my %HASH_info_this 	= %{$href_info};
		my %HASH_repo_this	= %{$href_repos_this};	
		
		#HASH
		my %HASH_code_words = ();
		$HASH_code_words{"short_tag"}		= 1;
		$HASH_code_words{"common_name"}		= 1;
		$HASH_code_words{"scientific_name"}	= 1;
		$HASH_code_words{"type"}			= 1;
		$HASH_code_words{"sequencing_depth"}= 1;
		$HASH_code_words{"k-mer"}			= 1;
		$HASH_code_words{"sequence_type"}	= 1;
		$HASH_code_words{"note"}			= 1;
		$HASH_code_words{"genome_size"}		= 1;
		
		my $INPUT_kindex_info 	= new IO::File($build_config, "r") or die "could not read file for repository information (parameter '--config'): $!\n";
		my $code_word 	= "";
		my $status 		= 0;
		while(<$INPUT_kindex_info>){
			my $line = $_;
			next if($line =~ /^$/);
			$line =~ s/\n$//;
			if($line =~ /^#/){
				$line =~ s/^#//;
				my @ARRAY_help = split(" ", $line);
				if(exists $HASH_code_words{$ARRAY_help[0]}){
					$code_word = $ARRAY_help[0];
					$HASH_info_this{$code_word} = "";
					$status = 1;
				}
			}else{
				$HASH_info_this{$code_word} = $line;
				$status = 0;
			}
		}
		
		#CHECK if mandatory field are set correctly
		my $check = 1;
		$check = 0 if($HASH_info_this{"short_tag"} eq "");
		$check = 0 if($HASH_info_this{"k-mer"} eq "");
		$check = 0 if($HASH_info_this{"common_name"} eq "");
		if($check == 0){
			#STOP
			print "\n\n WARNING: Kmasker (build) was stopped!!!\
				     \n Missing information in configuration!\n\n";
			exit(0);
		}
		
		#CEHCK for extisting entry
		my $short_tag = $HASH_info_this{"short_tag"};
		
		if(exists $HASH_repo_this{$short_tag}){
			#kindex with identical name exists, kmasker needs to be stopped to avoid overwriting!
			print "\n\n WARNING: Kmasker (build) was stopped!!!\
				     \n A kindex with identical name already exists in repository\n\n";
			exit(0);
		}		
		
		#RETURN
		return \%HASH_info_this;			
}


## subroutine
#
sub read_stats(){
	my $gs = $_[0];
	my $INPUT_stats 	= new IO::File("INPUT.fastq.stats", "r") or die "could not read INPUT.fastq.stats $!\n";
	my $calculation 	= 1;
	while(<$INPUT_stats>){
		my $line = $_;
		next if($line =~ /^$/);
		$line =~ s/\n$//;
		if($line =~ /^total bases/){
			my @ARRAY_tmp 	= split("\t", $line);
			$calculation 	= $gs * 1000000 / $ARRAY_tmp[1];
			print "\n CALC 	= ".$gs * 1000000 / $ARRAY_tmp[1];
			print "\n RES  	= ".$calculation;
		}
	}
	
	return $calculation;
}


## subroutine
#
sub write_repository_entry(){
	my $href_this 			= $_[0];
	my %HASH_info_this 		= %{$href_this};
	my @ARRAY_repository 	= ();
	$ARRAY_repository[0]	= $HASH_info_this{"short_tag"};
	$ARRAY_repository[1]	= $HASH_info_this{"common_name"};
	$ARRAY_repository[2]	= $HASH_info_this{"scientific_name"};
	$ARRAY_repository[3]	= $HASH_info_this{"type"};
	$ARRAY_repository[4]	= $HASH_info_this{"sequencing_depth"};
	$ARRAY_repository[5]	= $HASH_info_this{"k-mer"};
	$ARRAY_repository[6]	= $HASH_info_this{"sequence_type"};
	$ARRAY_repository[7]	= $HASH_info_this{"note"};
	$ARRAY_repository[8]	= $HASH_info_this{"md5sum"};
	$ARRAY_repository[9]	= $HASH_info_this{"seq"};
	$ARRAY_repository[10]	= $HASH_info_this{"absolut_path"};
	$ARRAY_repository[11]	= $HASH_info_this{"status"};
	
	my $new_kindex_entry_in_repository = join("\t", @ARRAY_repository);
	my $urepositories 		= "/home/".$HASH_info_this{"user_name"}."/.user_repositories.kmasker";
	open(my $RH, '>>', $urepositories) or die "Could not open file '$urepositories' $!";
	say $RH $new_kindex_entry_in_repository;
	close $RH;	
}


## subroutine
#
sub remove_repository_entry(){
	
	my $kindex_shorttag					= $_[0];
	my $href_this 						= $_[1];
	my %HASH_info_this 					= %{$href_this};
	my @ARRAY_repository_entries		= ();
	my $urepositories 					= $ENV{"HOME"}."/.user_repositories.kmasker";
	my $utmp							= "kmasker.tmp.file";
	my $target_entry_details			= "";
	
	#READ
	open(my $RH, '<', $urepositories) or die "Could not open file '$urepositories' $!";
	open(my $FH, '>', $utmp) or die "Could not write file ".$utmp." $!";
	while (<$RH>) {
		my @ARRAY_tmp = split("\t", $_);
		if($ARRAY_tmp[0] ne $kindex_shorttag){
  			print $FH $_;
		}else{
			$target_entry_details = $_;
		}
	}
	
	if($target_entry_details eq ""){
		print "\n WARNING: the kindex ".$kindex_shorttag." does not exist or is global (not permitted to remove)!\n\n";
		exit();
	}
	
	#MOVE
	system("mv ".$utmp." ".$urepositories);
	
	#REMOVE DATA
	chomp $target_entry_details;
	$target_entry_details =~ s/\n$//;	
	my @ARRAY_entry_details = split("\t", $target_entry_details);
	my $absolut_path = $ARRAY_entry_details[10];
	print "\n REMOVING kindex folder : ".$absolut_path."\n\n";
	system("rm -r ".$absolut_path);
	
	#CLEAN
	&clean_repository_directory($href_this);
	
	close $RH;	
	close $FH;
}


## subroutine
#
sub set_kindex_global(){
	
	my $kindex_shorttag				= $_[0];
	my $href_this 					= $_[1];
	my %HASH_info_this 				= %{$href_this};
	my @ARRAY_repository_entries	= ();	
	my $urepositories 				= $ENV{"HOME"}."/.user_repositories.kmasker";
	my $path 						= $HASH_info_this{"path_bin"};
	my $grepositories				= $path."/repositories.kmasker";
	my $path_kindex_global			= $HASH_info_this{"PATH_kindex_global"};
	my $utmp						= "kmasker_urep.tmp.file";
	
	#READ & EDIT user and global repository
	my $REPO_private 		= new IO::File($urepositories, "r") or die "\n unable to read global repository list $!";
	my $REPO_private_edit 	= new IO::File($utmp, "w") or die "\n unable to write tmp repository list $!";
	open(my $REPO_global, ">>", $grepositories) or die "\n unable to open global repository list in append mode $!";	
	while(<$REPO_private>){
		next if($_ =~ /^$/);
		if($_ =~ /^#/){
			print $REPO_private_edit $_;
			next;
		}
		my $line = $_;
		$line =~ s/\n$//;
		my @ARRAY_line = split("\t", $line);
		if($ARRAY_line[0] eq $kindex_shorttag){
			#ADD to GLOBAL
			my $old_path 	= $ARRAY_line[10];
			$ARRAY_line[10] = $path_kindex_global."KINDEX_".$kindex_shorttag."/";
			my $new_path 	= $ARRAY_line[10];
			$ARRAY_line[11] = "global";	#status
			$line 			= join("\t", @ARRAY_line);
			
			#ENTER to global repository
			print $REPO_global $line."\n";
			
			#COPY and REMOVE KINDEX
			system("mkdir ".$new_path);
			system("cp ".$old_path."* ".$new_path);
			system("rm -r ".$old_path);		
			
			print "\n Your requested kindex ".." was moved from private path ".$old_path." to global ".$new_path."\n";
		}else{
			#WRITE to tmp
			print $REPO_private_edit $line."\n";
		}
	}
	
	close $REPO_private;	
	close $REPO_private_edit;
	close $REPO_global;
	
	#MOVE
	system("mv ".$utmp." ".$urepositories);			
}


## subroutine
#
sub clean_repository_directory(){
	my $href_this 			= $_[0];
	my $href_repo_this 		= $_[1];
	my %HASH_info_this 		= %{$href_this};
	my %HASH_repo_info_this = %{$href_repo_this};
	my $path_kindex_private = $HASH_info_this{"PATH_kindex_private"};
	my $get = `ls $path_kindex_private`;
	my @ARRAY_list_of_kindex = split("\n", $get);
	foreach my $entry(@ARRAY_list_of_kindex){
		$entry =~ s/KINDEX_//;
		if(!(exists $HASH_repo_info_this{$entry})){
		#ask for cleaning
			system("rm -r ".$path_kindex_private."KINDEX_".$entry);
		}
	}
}

## subroutine
#
sub set_private_path(){
	my $private_path 	= $_[0];
	my $href_repo		= $_[1];
	my $old_private_path="";
	my $uconf 		= $ENV{"HOME"}."/.user_config.kmasker";
	my $USER_CONF_OLD 	= new IO::File($uconf, "r") or die "could not read old user conf : $!\n";
	while(<$USER_CONF_OLD>){
		my $line = $_;
		$line =~ s/\n//;
		my @ARRAY_tmp = split("\t", $line);
		$old_private_path	= $ARRAY_tmp[1] if($ARRAY_tmp[0] eq "PATH_kindex_private");
	}
	
	#SETUP user conf
	my $USER_CONF 	= new IO::File($uconf, "w") or die "could not write new user conf : $!\n";
	print $USER_CONF "PATH_kindex_private\t".$private_path."";
	close $USER_CONF;
	
	#MOVE and EDIT
	&move_private_structures($private_path, $old_private_path, $href_repo);
}

## subroutine
#
sub move_private_structures(){
	my $PPN = $_[0];
	my $PPO = $_[1];
	my $REP = $_[2];
	#EDIT PRIVATE REPOSITORY
	
	exit if($PPN eq $PPO);
	
	#load priavte user repositories
	my $REPO_private 		= new IO::File($ENV{"HOME"}."/.user_repositories.kmasker", "r") or die "\n unable to read user repository list $!";	
	my $REPO_private_NEW 	= new IO::File($ENV{"HOME"}."/.user_repositories.kmasker.new", "w") or die "\n unable to write user repository list $!";	
	while(<$REPO_private>){
		next if($_ =~ /^$/);
		if($_ =~ /^#/){
			print $REPO_private_NEW $_;
			next;
		}
		my $line = $_;
		$line =~ s/\n//;
		my @ARRAY_line = split("\t", $line);
		$ARRAY_line[10] = $PPN."KINDEX_".$ARRAY_line[0]."/";
		print $REPO_private_NEW join("\t", @ARRAY_line)."\n";
	}
	system("mv ".$ENV{"HOME"}."/.user_repositories.kmasker.new ".$ENV{"HOME"}."/.user_repositories.kmasker");	
	
	#MOVE data
	if(-d $PPN){
		print "\n ".$PPN." already exists and will be into ".$ENV{"HOME"}."/TMP_storage/";
		system("mv ".$PPN." ".$ENV{"HOME"}."/TMP_storage/");
		system("cp -r ".$PPO." ".$PPN);
	}else{
		system("cp -r ".$PPO." ".$PPN);
	}	
	system("rm -r ".$PPO);
		
}



1;