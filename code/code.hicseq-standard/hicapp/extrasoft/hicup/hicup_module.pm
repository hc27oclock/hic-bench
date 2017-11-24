#########################################################
#A collection of Perl subroutines for the HiCUP pipeline#
#########################################################

my $VERSION = '0.5.2';

use Data::Dumper;
use strict;
use warnings;
use FindBin '$Bin';
use lib $Bin;

###################################################################################
###################################################################################
##This file is Copyright (C) 2014, Steven Wingett (steven.wingett@babraham.ac.uk)##
##                                                                               ##
##                                                                               ##
##This file is part of HiCUP.                                                    ##
##                                                                               ##
##HiCUP is free software: you can redistribute it and/or modify                  ##
##it under the terms of the GNU General Public License as published by           ##
##the Free Software Foundation, either version 3 of the License, or              ##
##(at your option) any later version.                                            ##
##                                                                               ##
##HiCUP is distributed in the hope that it will be useful,                       ##
##but WITHOUT ANY WARRANTY; without even the implied warranty of                 ##
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  ##
##GNU General Public License for more details.                                   ##
##                                                                               ##
##You should have received a copy of the GNU General Public License              ##
##along with HiCUP.  If not, see <http://www.gnu.org/licenses/>.                 ##
###################################################################################
###################################################################################

###############################################################
#Sub: get_version
#Returns the version number of HiCUP
sub get_version {
    return $VERSION;
}

###############################################################
#Sub: hasval
#Takes a string and returns true (i.e. '1') if the string has a value
#(i.e. is not equal to nothing (''). This is useful since some
#variables may be set to nothing allowing them to be evaluated
#without initialisation errors, while simultaneously containing
#no information.
sub hasval {
    if ( $_[0] ne '' ) {
        return 1;
    } else {
        return 0;
    }
}


###############################################################
#Sub: hashVal
#Takes a hash and returns 1 if any of the hash's keys has a 
#value ne '' associated with it, else returns 0
sub hashVal {
    my %hash = @_;
    my $valueFound = 0;

    foreach my $key (keys %hash){
        $valueFound = 1 if hasval( $hash{$key} );
    }
    return $valueFound;
}


#Sub: deduplicate_array
#Takes and array and returns the array with duplicates removed
#(keeping 1 copy of each unique entry).
sub deduplicate_array {
    my @array = @_;
    my %uniques;

    foreach (@array) {
        $uniques{$_} = '';
    }
    my @uniques_array = keys %uniques;
    return @uniques_array;
}



#Sub: checkR
#Takes the config hash and and modifies $config{r} directly
#Script returns '0' if no path to R found
#The script checks whether the user specified path to R is valid and if 
#not tries to locate R automatically
sub checkR { 

    my $configHashRef = $_[0];

    return if($$configHashRef{r} eq '0');   #Already determined R not present, so don't repeat warnings

    if( hasval ( $$configHashRef{r} )   ){    #Check whether user specified a path
        if(-e $$configHashRef{r}){      #Check if R path exists

            #Check R runs   
            my $versionR = `$$configHashRef{r} --version 2>/dev/null`;
            unless($versionR =~ /^R version/){
                warn "R not found at '$$configHashRef{r}'\n";
                $$configHashRef{r} = '';
            }

        }else{
            warn "'$$configHashRef{r}' is not an R executable file\n";
            $$configHashRef{r} = '';
        }
    }

    unless(hasval $$configHashRef{r} ){
        warn "Detecting R automatically\n";
        if ( !system "which R >/dev/null 2>&1" ) {
            $$configHashRef{r}  = `which R`;
            chomp $$configHashRef{r};
            warn "Found R at '$$configHashRef{r}'\n";
        } else {
            warn "Could not find R (www.r-project.org), please install if graphs are required\n";
            $$configHashRef{r} = 0;    #Tells later scripts that R is not installed
        }
    }
}


#Sub: outdirFileNamer
#Takes an array of filenames as a reference and the output file directory name/path
#returns a hash of %{path/filename} = outdir/filename
sub outdirFileNamer {
    my $fileArrayRef = $_[0];    #Passed by reference
    my $outdir       = $_[1];
    my %inOutFilenames;

    $outdir .= '/' unless ( $outdir =~ /\/$/ );    #Ensure outdir ends with forward slash

    foreach my $file (@$fileArrayRef) {
        my @elements = split( '/', $file );
        my $outFile = $outdir . $elements[-1];
        $inOutFilenames{$file} = $outFile;
    }
    return %inOutFilenames;
}

############################
#Subroutine "process_config":
#Takes i) configuration file name and ii) %config hash (as a reference).
#The script then uses the configuration file to populate the hash as
#appropriate. Parameters passed via the command line take priority
#over those defined in the configuration file.
#The script modifies the hash directly, but returns as a hash the filename pairs (i.e the 
#lines in the configuration #file that could did not correspond configuration parameters
#Each file is on a separate line with pairs on adjacent lines, or a pair may be placed on
#the same line separated by pipe ('\')
sub process_config {
    my ( $config_file, $config_hash_ref ) = @_;
    my @non_parameters;    #Stores lines in the configuration file not defined as parameters

    #Open configuration file
    open( CONF, "$config_file" ) or die "Can't read $config_file: $!";

    while (<CONF>) {

        my $line = $_;
        chomp $line;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;    #Remove starting/trailing white spaces

        next if $line =~ /^\s*\#/;    #Ignore comments
        next if $line =~ /^\s*$/;     #Ignore whitespace

        #Check if this is a parameter
        my ( $parameter, $setting ) = split( /:/, $line );
        $parameter =~ s/^\s+//;
        $parameter =~ s/\s+$//;       #Remove starting/trailing white spaces
        $parameter = lc $parameter;
        $setting =~ s/^\s+// if defined ($setting);
        $setting =~ s/\s+$// if defined ($setting);         #Remove starting/trailing white spaces

        if ( exists $$config_hash_ref{$parameter} ) {
            if ( $$config_hash_ref{$parameter} eq '' ) {    #Check parameter not assigned value in command line
                $$config_hash_ref{$parameter} = $setting;    #Edit the configuration hash
            }
        } else {
            my @lineElements = split(/\|/, $line);    #There may be a pipe separating two files
            foreach my $element (@lineElements){
                $element =~ s/^\s+//;
                $element =~ s/\s+$//;       #Remove starting/trailing white spaces
                push( @non_parameters, $element );
            }
        }
    }
    close CONF or die "Could not close filhandle on configuration file: '$config_file'\n";

    if(scalar @non_parameters % 2){
      die "There needs to be an even number of files in the configuration file, see hicup --help for more details.\n";
    }
    return @non_parameters;
}

# Subroutine: check_filenames_ok
# Receives an array of filename pairs delimited by comma ','
# Checks the files exist and returns a hash of %{forward file} = reverse file1F
sub check_filenames_ok {
    my @files = @_;
    my %paired_filenames;

    foreach (@files) {
        my @file_combination = split /,/;
        if ( scalar @file_combination != 2 ) {
            die "Files need to be paired in the configuration file an/or command-line, see hicup --help for more details.\n";
        }

        foreach my $file (@file_combination) {
            if ( $file eq '' ) {
                die "Files need to be paired in the configuration file and/or command-line, see hicup --help for more details.\n";
            }
            $file =~ s/^\s+//;    #Remove white space at start and end
            $file =~ s/\s+$//;
        }
        $paired_filenames{ $file_combination[0] } = $file_combination[1];
    }
    return %paired_filenames;
}



# Subroutine: check_filenames_ok
# Receives a hash of filename pairs %{forward file} = reverse file1F
# Checks that no filename occurs more than once, irrespective of path
sub check_no_duplicate_filename{
	my %filePair = @_;
	my %uniqueNames;
	my %duplicateNames;    #Write here names that occurred multiple times
	my $ok = 1;    #Is the configuration acceptable?
	
	foreach my $key (keys %filePair){
		my $fileF = (split(/\//, $key))[-1];
		my $fileR = $filePair{$key};
		$fileR = (split(/\//, $fileR))[-1];
		
		foreach my $file ($fileF, $fileR){
			if(exists $uniqueNames{$file}){
				$duplicateNames{$file} = '';
				$ok = 0;
			}else{
				$uniqueNames{$file} = '';
			}			
		}	
	}
	
	unless($ok){
		foreach my $duplicateName (keys %duplicateNames){
			warn "Filename '$duplicateName' occurs multiple times\n";
		}	
	}

	return $ok;
}



# Subroutine: checkAligner
# Receives the config hash and determines whether the specified aligners
# are present, if not the aligners are searched for automatically and 
# the config hash is adjusted accordingly
# Returns 1 if successful or 0 if no valid aligner found
sub checkAligner{

	my $configRef = $_[0];
    my @aligners = ( 'bowtie', 'bowtie2' );    #List of aligners used by HiCUP
	my $parameters_ok = 1;

    #Check aligners
    my $found_aligner_flag = 0;
    my $aligner_count = 0;
    foreach my $aligner_name (@aligners) {
        if ( hasval( $$configRef{$aligner_name}) ) {
            $$configRef{aligner} = $aligner_name;
            $aligner_count++;
        }
    }

    if ( ( $aligner_count == 0 ) ) {    #Find an aligner if none specified
        warn "No aligner specified, searching for aligner\n";
        foreach my $aligner_name (@aligners) {
            if ( !system "which $aligner_name >/dev/null 2>&1" ) {
                my $aligner_path = `which $aligner_name`;
                chomp $aligner_path;
                warn "Path to $aligner_name found at: $aligner_path\n";
                $found_aligner_flag = 1;
				$$configRef{$aligner_name} = $aligner_path;
				$$configRef{aligner} = $aligner_name;
                last;
            } else {
                warn "Could not find path to '$aligner_name'\n";
            }
        }
    }

    if ( $aligner_count == 1 ) {    #Correct number (i.e. one) of aligners specified, check path correct
        my $aligner_path;
		my $aligner_name = $$configRef{aligner};
        if ( !system "which $aligner_name >/dev/null 2>&1" ) {
            my $aligner_path = `which $aligner_name`;
            chomp $aligner_path;
            print "Using $aligner_name at '$aligner_path'\n";
			$$configRef{$aligner_name} = $aligner_path;
            $found_aligner_flag = 1;
        }

        unless ($found_aligner_flag) {
			my $aligner_name = $$configRef{aligner};
            warn "Could not find $aligner_name at: $$configRef{$aligner_name}\n";
            $$configRef{$aligner_name} =~ s/\/$//;    #Remove final '/' from path
            $$configRef{$aligner_name} .= '/$aligner_name';
            warn "Trying $$configRef{$aligner_name}\n";
            if ( !system "which $aligner_name >/dev/null 2>&1" ) {
                $$configRef{$aligner_name} = `which $aligner_name`;
                print "Print found $aligner_name at '$$configRef{$aligner_name}'\n";
                $found_aligner_flag = 1;
            }
        }

        unless ($found_aligner_flag) {
			my $aligner_name = $$configRef{aligner};
            warn "Could not find $aligner_name at: $$configRef{$aligner_name}\n";
            if ( !system "which $aligner_name >/dev/null 2>&1" ) {
                $$configRef{$aligner_name} = `which $aligner_name`;
                chomp $$configRef{$aligner_name};
                warn "Path to $aligner_name found at: $$configRef{$aligner_name}\n";
                $found_aligner_flag = 1;
            } else {
                warn "Could not find $aligner_name using '$$configRef{$aligner_name}'\n";
            }
        }
    }

    if ( $aligner_count > 1 ) {    #Too many aligners specified (i.e. more than 1)
        warn "Please only specify only one aligner: either --bowtie or --bowtie2.\n";
        $parameters_ok = 0;
    }

    unless ($found_aligner_flag) {
        warn "Please specify a link to one valid aligner\n";
        $parameters_ok = 0;
    }


    if ( $$configRef{ambiguous} ) {
        unless( hasval($$configRef{bowtie2}) ) {
           warn "Option 'ambiguous' is only compatible wtih Bowtie2\n";
            $parameters_ok = 0;
        }
    }
	return $parameters_ok;
}



# Subroutine: checkAlignerIndices
# Receives the config hash and determines whether the specified indices
# are present.
# Returns 1 if successful or 0 if no valid aligner found	
sub checkAlignerIndices{
	my $configRef = $_[0];
	my $parameters_ok = 1;

	 #Check the index files exist
	if( hasval $$configRef{index} ){		
		my @index_suffixes;
		if ( $$configRef{aligner} eq 'bowtie' ) {
			@index_suffixes = ( '.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt' );
		} elsif ( $$configRef{aligner} eq 'bowtie2' ) {
			@index_suffixes = ( '.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2' );
		}

		foreach my $suffix (@index_suffixes) {
			my $indexFilename = $$configRef{index} . $suffix;
			unless( ( -e $indexFilename ) or ( -e $indexFilename . 'l' ) ){    #Bowtie2 also has larger indices
				warn "$$configRef{aligner} index file '$indexFilename' does not exist\n";
				$parameters_ok = 0;
			}
		}
	}else{
		warn "Please specify alinger indices (--index)\n";
		$parameters_ok = 0;
	}
	return $parameters_ok;	
}



###################################################################
#check_files_exist:
#Takes a reference to an array containing paths to filenames and verifies they exist
#Warns of files that do no exit. Returns 1 if all files exist but 0 if this is not
#the case.
#
#Also, takes a second argument:
#$_[1] should be 'EXISTS' or 'NOT_EXISTS'
#If 'NOT_EXIST' warns if file already exists.  Returns '1' if none of the
#files exists and '0' if one or multiple files already exist
sub check_files_exist {
    my $files      = $_[0];    #Reference to array
    my $check_for  = $_[1];
    my $all_exist  = 1;
    my $not_exists = 1;

    if ( $check_for eq 'EXISTS' ) {
        foreach my $file (@$files) {
            unless ( -e $file ) {
                warn "$file does not exist\n";
                $all_exist = 0;
            }
        }
    } elsif ( $check_for eq 'NOT_EXISTS' ) {
        foreach my $file (@$files) {
            if ( -e $file ) {
                warn "$file already exists\n";
                $not_exists = 0;
            }
        }
    } else {
        die "Subroutine 'check_files_exist' requires argument 'EXISTS' or 'NOT_EXISTS'.\n";
    }

    if ( $check_for eq 'EXISTS' ) {
        return $all_exist;
    } else {
        return $not_exists;
    }
}

#############################################################
#datestampGenerator
#Returns a suitably formatted datestamp
sub datestampGenerator {
    my @now       = localtime();
    my $datestamp = sprintf(
        "%02d-%02d-%02d_%02d-%02d-%04d",

        $now[2], $now[1],     $now[0],
        $now[3], $now[4] + 1, $now[5] + 1900
    );
    return $datestamp;
}


##################################################################################################
#Renaming Files
##################################################################################################



# Subroutine 'determine_truncater_outfiles':
# Takes as references: 1) hash of forward/reverse file pairs (i.e. %{fileF} = fileR)
# and 2) the configuration hash and 3) flag (1 or 0) stating whether to return summary file names etc
# Returns an array of the file that will be generated by hicup_truncater
sub determine_truncater_outfiles {
    my (%input_files) = @_;
    my %output_files;
	
	my $fileHashRef   = $_[0];
    my $configHashRef = $_[1];
    my @outfiles;

    for my $filename1 ( keys %input_files ) {
        my $filename2 = $input_files{$filename1};

        $filename1 =~ s/\.gz$|\.bz2$//;
        $filename1 =~ s/^.+\///;            #Remove folder references
        $filename1 =~ s/\.fq|\.fastq$//;    #Remove .fq or .fastq file extension
        $filename1 .= "_trunc";

        $filename2 =~ s/\.gz$|\.bz2$//;
        $filename2 =~ s/^.+\///;            #Remove folder references
        $filename2 =~ s/\.fq|\.fastq$//;    #Remove .fq or .fastq file extension
        $filename2 .= "_trunc";

        if ( $$configHashRef{zip} and $$configHashRef{keep} ) {
            $filename1 .= ".gz";
            $filename2 .= ".gz";
        }
        $output_files{$filename1} = $filename2;
    }

    print Dumper \%output_files;

    return %output_files;

}



# Subroutine 'determine_mapper_outfiles':
# Takes as references: 1) hash of forward/reverse file pairs (i.e. %{fileF} = fileR)
# and 2) the configuration hash
# Returns an array of the file that will be generated by hicup_mapper
sub determine_mapper_outfiles {
    my $fileHashRef   = $_[0];
    my $configHashRef = $_[1];
    my @outfiles;

    foreach my $fileForward ( keys %$fileHashRef ) {
        my $fileReverse = $$fileHashRef{$fileForward};
        foreach my $file ( $fileForward, $fileReverse ) {

            #Determine the .map files
            $file =~ s/^.+\///;    #Remove folder references
            $file =~ s/\.gz$//;
            $file .= ".map";
            $file = $$configHashRef{outdir} . $file;
            push( @outfiles, $file );
            $file =~ s/^.+\///;    #Remove folder path again to assist with .pair filename generation below
        }

        #Determine the .pair files
        my $pair_filename = "$$configHashRef{outdir}" . "$fileForward.$fileReverse.pair";
        $pair_filename .= '.gz' if $$configHashRef{zip};
        push( @outfiles, $pair_filename );

        if ( $$configHashRef{ambiguous} ) {
            my $ambiguous_pair_filename = "$$configHashRef{outdir}.$fileForward.$fileReverse.ambiguous.pair";
            $ambiguous_pair_filename .= '.gz' if $$configHashRef{zip};
            push( @outfiles, $ambiguous_pair_filename );
        }
    }

    my $summaryfile = "$$configHashRef{outdir}" . "hicup_mapper_summary.$$configHashRef{datestamp}.txt";
    push( @outfiles, ( $summaryfile, "$summaryfile.temp" ) );    #Check for summary files
    return @outfiles;
}

# Subroutine 'determine_filter_outfiles':
# Takes as references: 1) array of input filenames
# and 2) the configuration hash
# 3) flag (1 or ) stating whether to return summary file names etc
# Returns an array of the file that will be generated by hicup_filter
sub determine_filter_outfiles {
    my $fileArrayRef          = $_[0];
    my $configHashRef         = $_[1];
    my $reportAdditionalFiles = $_[2];
    my @outfiles;

    foreach (@$fileArrayRef) {

        my $file = $_;

        $file =~ s/^.+\///;    #Remove folder references
        if ( hasval $$configHashRef{outdir} ) {
            $file = $$configHashRef{outdir} . $file;
        }

        if ( $$configHashRef{zip} and $$configHashRef{samtools} ) {
            push( @outfiles, $file . '.bam' );
        } elsif ( $$configHashRef{zip} ) {
            push( @outfiles, $file . '.sam.gz' );
        } else {
            push( @outfiles, $file . '.sam' );
        }

        if ($reportAdditionalFiles) {    #Report summary files etc.
            push( @outfiles, $file . '_ditag_classification.png' );
            push( @outfiles, $file . '_ditag_size_distribution.png' );
            push( @outfiles, $file . 'ditag_lengths.' . $$configHashRef{datestamp} . '.temp' );
        }
    }

    push( @outfiles, $$configHashRef{outdir} . 'hicup_filter_summary_' . $$configHashRef{datestamp} . '.txt' ) if ($reportAdditionalFiles);
    return @outfiles;
}

# Subroutine 'determine_deduplicator_outfiles':
# Takes as references: 1) array of input filenames
# and 2) the configuration hash
# 3) flag (1 or ) stating whether to return summary file names etc
# Returns an array of the file that will be generated by hicup_filter
sub determine_deduplicator_outfiles {

    my $fileArrayRef          = $_[0];
    my $configHashRef         = $_[1];
    my $reportAdditionalFiles = $_[2];
    my @outfiles;

    foreach (@$fileArrayRef) {
        my $file = $_;         #Don't modify original array value
        $file =~ s/^.+\///;    #Remove folder references
        $file =~ s/\.gz$//;    #Remove the gzip file extension from the ouptfile name
        $file =~ s/.sam$//;    #Remove '.sam' file extension
        $file =~ s/.bam$//;    #Remove '.bam' file extension

        if ( hasval $$configHashRef{outdir} ) {
            $file = $$configHashRef{outdir} . $file;
        }

        if ( $$configHashRef{zip} and $$configHashRef{samtools} ) {
            $file .= '.bam';
        } elsif ( $$configHashRef{zip} ) {
            $file .= '.sam.gz';
        } else {
            $file .= '.sam';
        }

        $file = 'uniques_' . $file;
		
        push( @outfiles, $file );
    }
	
	push( @outfiles, $$configHashRef{outdir} . 'hicup_deduplicator_summary_' . $$configHashRef{datestamp} . '.txt' ) if ($reportAdditionalFiles);
    
	return @outfiles;
}



######################
#Subroutine "newopen":
#links a file to a filehandle
sub newopen {
    my $path = shift;
    my $fh;

    open( $fh, '>', $path ) or die "\nCould not create filehandles in subroutine \'newopen\'\n";

    return $fh;
}


##############################
#Subroutine 'quality_checker':
#determines the FASTQ format of a sequence file
#
#FASTQ FORMAT OVERVIEW
#---------------------
#Sanger: ASCII 33 to 126
#Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126 
#(although in raw read data the Phred quality score rarely exceeds 60, higher 
#scores are possible in assemblies or read maps)
#
#Solexa/Illumina 1.0 format: ASCII 59 to 126
#-5 to 62 using ASCII 59 to 126 (although in raw read data Solexa scores from -5 
#to 40 only are expected)
#
#Illumina 1.3 and before Illumina 1.8: ASCII 64 to 126
#the format encoded a Phred quality score from 0 to 62 using ASCII 64 to 126 
#(although in raw read data Phred scores from 0 to 40 only are expected). 
#
#Illumina 1.5 and before Illumina 1.8: ASCII 66 to 126
#the Phred scores 0 to 2 have a slightly different meaning. The values 0 and 
#1 are no longer used and the value 2, encoded by ASCII 66 "B", is used also 
#at the end of reads as a Read Segment Quality Control Indicator.
#
#phred64-quals: 
#ASCII chars begin at 64
#
#Starting in Illumina 1.8, the quality scores have basically returned to
#Sanger format (Phred+33)
#
#solexa-quals: ASCII chars begin at 59
#integer-qual: quality values integers separated by spaces
sub quality_checker {
    my $file       = $_[0];
    my $score_min  = 999;     #Initialise at off-the-scale values
    my $read_count = 1;

    if ( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Could not read file '$file' : $!";
    } else {
        open( IN, $file ) or die "Could not read file '$file' : $!";
    }

    while (<IN>) {

        next if (/^\s$/);     #Ignore blank lines

        if (/^@/) {
            scalar <IN>;
            scalar <IN>;

            my $quality_line = scalar <IN>;
            chomp $quality_line;
            my @scores = split( //, $quality_line );

            foreach (@scores) {
                my $score = ord $_;    #Determine the value of the ASCII character

                if ( $score < $score_min ) {
                    $score_min = $score;
                }

                if($score_min < 59){
                    return 'Sanger';
                }
            }
        }
        $read_count++;
    }

    close IN or die "Could not clode filehandle on '$file' : $!";

    if ( $read_count < 1_000_000 ) {
        return 0;    #File did not contain enough lines to make a decision on quality
    } 

    if($score_min < 64){
        return 'Solexa_Illumina_1.0'
    }elsif($score_min < 66){
        return 'Illumina_1.3'
    }else{
        return 'Illumina_1.5';
    }

}



################################
#Subroutine
#Receives the FASTQ fomrat and the aligner and determines the aligner-specific format flag
#Input values are Sanger, Solexa_Illumina_1.0, Illumina_1.3, Illumina_1.5 for the FASTQ fromat
#and bowtie or bowtie2 for the aligner
#If only the FASTQ format is specified 'NO_ALIGNER' will be returned, so the subroutine can
#be used to check whether the FASTQ format is valid.
sub determineAlignerFormat{

    my ($fastqFormat, $aligner) = @_;
    $fastqFormat = uc $fastqFormat;

    unless($fastqFormat =~ /^SANGER$|^SOLEXA_ILLUMINA_1.0$|^ILLUMINA_1.3$|^ILLUMINA_1.5$/){
        warn "'$fastqFormat' is not a valid FASTQ format (valid formats changed in HiCUP v0.5.2)\n";
        warn "Valid formats are: 'Sanger', 'Solexa_Illumina_1.0', 'Illumina_1.3' or 'Illumina_1.5'\n";
        return 0;
    }

    unless(defined $aligner){    #By returning this message if no aligner specified, the 
        return 'NO_ALIGNER';
    }

    if($aligner eq 'bowtie' ){
        if ($fastqFormat eq 'SANGER'){
            return 'phred33-quals';
        }elsif($fastqFormat eq 'SOLEXA_ILLUMINA_1.0'){
            return 'solexa-quals';
        }elsif($fastqFormat  eq 'ILLUMINA_1.3' ){
            return 'phred64-quals';
        }elsif($fastqFormat  eq 'ILLUMINA_1.5' ){
            return 'phred64-quals';
        }
    }

    if($aligner eq 'bowtie2'){
        if ($fastqFormat eq 'SANGER'){
            return 'phred33';
        }elsif($fastqFormat eq 'SOLEXA_ILLUMINA_1.0'){
            return 'solexa-quals';
        }elsif($fastqFormat  eq 'ILLUMINA_1.3' ){
            return 'phred64';
        }elsif($fastqFormat  eq 'ILLUMINA_1.5' ){
            return 'phred64';
        }
    }



}





#Subroutine "print_example_config_file"
#Takes the name of the config file and then copies to the current working directory
sub print_example_config_file{
    my $file = $_[0];
    my $fileAndPath = "$Bin/config_files/$file";
    !system("cp $fileAndPath .") or die "Could not create example configuratation file '$file' : $!";
    print "Created example configuration file '$file'\n";

}




1
