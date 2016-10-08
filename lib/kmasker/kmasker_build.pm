package kmasker::kmasker_build;
use Exporter qw(import);
use File::Basename;
use strict;
use warnings;
use Cwd;

our @ISA = qw(Exporter);
our @EXPORT = qw(
setup_index
build_kindex_jelly
make_config
);
our @EXPORT_OK = qw(build_kindex_jelly make_config);


sub build_kindex_jelly{	
	my $href_info		= $_[0];
	my $build_config 	= $_[1];
	my $href_repos		= $_[2];
	my $store_input		= $_[3];	
	my %HASH_info 		= %{$href_info};
	my %HASH_repo		= %{$href_repos};
	my $user_name 		= $HASH_info{"user_name"};
	my $seq				= $HASH_info{"seq"};
	my $k				= $HASH_info{"k-mer"};
	
	#LOAD info
	if(defined $build_config){
		#READ config for build
		my $href_this = &read_config($build_config, \%HASH_info, $href_repos);	
		%HASH_info = %{$href_this};	
	}
	
	#Create single input file and calculate md5sum for repository
	system("cat ".$seq." >INPUT.fastq");
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
	$seq				= "INPUT_".$md5sum.".".$end;
	system("mv INPUT.fastq INPUT_".$md5sum.".".$end);
		
	#FILL with default if inforomation is not provided	
	$HASH_info{"absolut_path"} = getcwd if(!exists $HASH_info{"absolut_path"});
	$HASH_info{"sequencing_depth"} = 1 if(!exists $HASH_info{"sequencing_depth"});
	
	#BUILD JELLY index	
	print "\n ... start construction of kindex\n"; 
	system("jellyfish count -m ".$k." -s 24G -t 10 -o KINDEX_".$HASH_info{"short_tag"}."_".$md5sum."_k".$k.".jf ".$seq."");
	print "\n ... finished kindex construction!\n";
	
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
		print $USER_kindex_info "#species [mandatory]\n\n";
		#3
		print $USER_kindex_info "#botanic_name [obligarory]\n\n";
		#4
		print $USER_kindex_info "#type (cultivar;genotype;etc.) [obligarory]\n\n";
		#5
		print $USER_kindex_info "#sequencing_depth [mandatory]\n\n";
		#6
		print $USER_kindex_info "#k-mer [obligatory]\n\n";
		#7
		print $USER_kindex_info "#sequence_type (reads or assembly) [obligatory]\n\n";
		#8
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
		$HASH_code_words{"species"}			= 1;
		$HASH_code_words{"botanic_name"}	= 1;
		$HASH_code_words{"type"}			= 1;
		$HASH_code_words{"sequencing_depth"}= 1;
		$HASH_code_words{"k-mer"}			= 1;
		$HASH_code_words{"sequence_type"}	= 1;
		$HASH_code_words{"note"}			= 1;
		
		my $INPUT_kindex_info 	= new IO::File($build_config, "r") or die "could not write file for repository information: $!\n";
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
		$check = 0 if($HASH_info_this{"species"} eq "");
		if($check == 0){
			#STOP
			print "\n\n WARNING: Kmasker (build) was stopped!!!\
				     \n Missing informatio nin configuration!\n\n";
			exit(0);
		}
		
		#CEHCK for extisting entry
		my $short_tag = $HASH_info_this{"short_tag"};
		
		#TEST
#		foreach(keys %HASH_repo_this){
#			print "\n ENTRY in repo = ".$_;
#		}		
		#END TEST
		
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
sub write_repository_entry(){
	my $href_this 			= $_[0];
	my %HASH_info_this 		= %{$href_this};
	my @ARRAY_repository 	= ();
	$ARRAY_repository[0]	= $HASH_info_this{"short_tag"};
	$ARRAY_repository[1]	= $HASH_info_this{"species"};
	$ARRAY_repository[2]	= $HASH_info_this{"botanic_name"};
	$ARRAY_repository[3]	= $HASH_info_this{"type"};
	$ARRAY_repository[4]	= $HASH_info_this{"sequencing_depth"};
	$ARRAY_repository[5]	= $HASH_info_this{"k-mer"};
	$ARRAY_repository[6]	= $HASH_info_this{"sequence_type"};
	$ARRAY_repository[7]	= $HASH_info_this{"note"};
	$ARRAY_repository[8]	= $HASH_info_this{"md5sum"};
	$ARRAY_repository[9]	= $HASH_info_this{"seq"};
	$ARRAY_repository[10]	= $HASH_info_this{"absolut_path"};
	
	my $new_kindex_entry_in_repository = join("\t", @ARRAY_repository);
	my $urepositories 		= "/home/".$HASH_info_this{"user_name"}."/.user_repositories.kmasker";
	open(my $RH, '>>', $urepositories) or die "Could not open file '$urepositories' $!";
	say $RH $new_kindex_entry_in_repository;
	close $RH;	
}

1;