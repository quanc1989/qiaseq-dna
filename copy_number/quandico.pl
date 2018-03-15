#!/usr/bin/perl -w
use strict;
use Config::IniFiles qw<ParseConfig SaveConfig>;
use Data::Dump::Color qw<dd>;
use Data::Table::Excel qw<tables2xlsx>;
use Getopt::Long::Descriptive;
use Hash::Merge qw<merge>;
use File::chdir;
use File::Copy qw<copy move>;
use File::Slurp qw<read_file>;
use File::Spec::Functions qw<rel2abs catfile tmpdir>;
use version 0.77;
our $VERSION = version->declare(0.8.5);    ## GIT_VERSION ## our $VERSION = version->declare(<%>);
our $BRANCH  = 'develop';                  ## GIT_BRANCH ## our $BRANCH = '<%>';

# run the demo input like this
# quandico -s data=data/demo_sample_male.cov -s x=1 -s y=1 -r data=data/demo_reference_female.cov -r x=2 -r y=0
#
# define variables and some defaults
my ( $folder, $config, %reference, %sample, %global, %options );

# Getopt::Long::Descriptive - input parameters
my ( $opt, $usage ) = describe_options(
	'quandico %o <options>',
	['
  input file and global options:'
	],
	['sample|s=s%',    'sample details (hash: data, name, x, y)'],
	['reference|r=s%', 'reference details (hash: data, name, x, y)'],
	['genome|G=s%', 'genome details (hash: fai, data, name)', {default => {}}],
	['global|g=s%', 'global options (title, see manual for more)'],
	[
		'cluster|l!',
		'perform clustering of regions (on), using <-E> with <-a>, used per default, say --no-cluster to switch off',
		{default => 1}
	],    # this is a new feature from 0.8.1+
	['
  output files:'
	],
	['destdir|d=s',  'destination directory',     {default => rel2abs($CWD)}],
	['tmpdir|t=s',   'working directory',         {default => tmpdir()}],
	['basename|b=s', 'basename for output files', {default => '<generate>'}],
	['
  config file / batch mode:'
	],
	['config|c=s@', 'configuration file (.ini) with parameters'],
	['folder|f=s',  'folder containing .ini files (batch mode)'],
	['
  clustering parameters:'
	],    # this group of options is new from 0.8.1+
	['engine|E=s', 'script to run for clustering (cluster)', {default => 'cluster'}],
	['cp|C=s%',    'hash with parameters for clustering (see manual of engine)'],
	['
  internals:'
	],
	['rexe|x=s',  'R executable (default: /usr/bin/R)',      {default => '/usr/bin/R'}],
	['rpath|y=s', 'path where the R script resides in',      {default => '/srv/qiagen/software/quandico/R/'}],
	['rcall|z=s', 'R script to run (default: < quandico.R)', {default => ' < /srv/qgen/code/qiaseq-dna/copy_number/R/quandico.R'}],
	['mock|m',    'do not run R, only print the call'],
	['
  standard options:'
	],
	['verbose+', 'control verbosity (off)', {default => 0}],
	['version',  'show the version and exit'],
	['history',  'show the release notes and changes'],
	['quiet|q', 'suppress all output except errors', {implies     => {verbose => 0}}],
	['debug',   'print debug data (development only)'],
	['help|h',  'show this help screen'],
	['
  
# usage with demo data (included)
./quandico -s data=data/demo_sample_male.cov -s x=1 -s y=1 -r data=data/demo_reference_female.cov -r x=2 -r y=0 --cp names=data/refGene.txt
'
	]
);
print( $usage->text, 1 ), exit if $opt->help;
if ( $opt->version or $opt->history ) {
	my $name = $0;
	$name =~ s!^\.\/!!;
	printf "This is %s %s %s\nUse flags -h or --help for help.\n", $name, $VERSION->stringify, $BRANCH ? qq~($BRANCH)~ : '';
	print q~
Version  Notes

0.8.6    

         The chromosome details can now be provided as .fai index file:
         --genome data=/path/to/mygenome.fa.fai. The filename will be used 
         as default assembly ("mygenome" in this case) but setting a name
         with --genome name=newname has precedence over this.

0.8.4    author  Frank Reinecke <frank.reinecke@qiagen.com>	
         date    Fri, 28 Mar 2014 15:06:16 +0000 (16:06 +0100)
 
         Fixed bug that failed to write Excel files to the correct place.
         Added a new set of options to use custom chromsome files.
         Changed the weighted variance code to a better suited one.

0.8.2    author Frank Reinecke <frank.reinecke@qiagen.com>	
         date   Tue, 25 Mar 2014 13:43:42 +0000 (14:43 +0100)

         Updates to output formats: Added xlsx export. New functionality to
         preprocess input data using a newly developed clustering script.
		
0.8.0    author Frank Reinecke <frank.reinecke@qiagen.com>	
         date   Thu, 13 Mar 2014 09:15:32 +0000 (10:15 +0100)

         Initial release with basic functionaly implemented.

~ if ( $opt->history );
	exit;
} ## end if ( $opt->version or $opt->...)
;

# read command line
process_options();

# collect input
my @INI = ( '<args>', $opt->{config} );
if ( $opt->{folder} and not $opt->{config} ) {
	opendir( my $dir, $opt->{folder} );
	@INI = map { $opt->{folder} . '/' . $_ } grep { /\.ini$/ } readdir($dir);
}

# generate a name from samples names unless something was set
if ( $opt->{basename} eq '<generate>' ) {
	$opt->{basename} = join( '_', $opt->{sample}{name}, $opt->{reference}{name} );
}

# loop through all config files collected (or use command line parameters)
foreach my $config_file (@INI) {
	last unless $config_file;
	if ( $config_file ne '<args>' and -e $config_file ) {
		tie %{$opt}, 'Config::IniFiles', ( -file => $config_file, -fallback => 'global' );
		copy( $config_file, $opt->{global}{title} . '.ini' );
		dd $opt;
		exit;
	}
	## SANITY CHECKS BELOW:
	##
	## 1. check for all required input FILES
	foreach my $group (qw<sample reference>) {
		foreach my $type (qw<data>) {
			my $file = $opt->{$group}{$type} || '';
			if ( not defined $file or not -e $file ) {
				print STDERR "
ERROR: Input data file not readable.
   Group: <$group>
   Type:  <$type>
   Value: <$file>
   
";
				$usage->die;
			} ## end if ( not defined $file or not...)

			# clustering to be performed here
			if ( $opt->{cluster} ) {
				$file = rel2abs($file);
				my $sink           = $file . '.clustered';
				my $clust_cmd      = join( " ", $opt->{clustering}, '--input', $file, '>', $sink );
				my $cluster_output = `perl $clust_cmd`;
				$opt->{$group}{$type} = $sink;
			}
		} ## end foreach my $type (qw<data>)
	} ## end foreach my $group (qw<sample reference>)
	#my $temp_basename = join( '_', 'quandico', $$, random_chars(5) );
   my $temp_basename = $opt->{basename};   # hack 2018.02.18 JDiCarlo - derive readset name from prefix option
   $temp_basename =~ s/\.copy-number//;    # hack 2018.02.18 JDiCarlo - derive readset name from prefix option
	my $chromosomes = $opt->{genome}{data} ? rel2abs($opt->{genome}{data}) : '';

	# compile arguments for R
	my @RArgs = (
		map { qq~"$_"~ } (
			rel2abs( $opt->{sample}{data} ),       #  4 = absolute path of input for sample
			$opt->{sample}{x},                     #  5 = x of sample
			$opt->{sample}{y},                     #  6 = y of sample
			rel2abs( $opt->{reference}{data} ),    #  7 = absolute path of input for reference
			$opt->{reference}{x},                  #  8 = x of reference
			$opt->{reference}{y},                  #  9 = y of reference
			$opt->{tmpdir},                        # 10 = temporary directory to work in
			$temp_basename,                        # 11 = basename to use for files
			$chromosomes,                          # 12 = file with chromosome details (optional)
			$opt->{genome}{name}                   # 13 = config.assembly.version
		)
	);
	my $redirect = $opt->{verbose} ? '' : '2>/dev/null';
	#my $Rcall = join( " ", $opt->{rexe}, '--vanilla', '--slave', @RArgs, $opt->{rcall}, $redirect );
        my $Rcall = join( " ", $opt->{rexe}, '--vanilla', '--slave', @RArgs, $opt->{rcall} );

	# prepare the call to R
	RUN_R: {
		#local $CWD = $opt->{rpath};
		if ( $opt->{mock} ) {
			print $Rcall, "\n";
		}
		else {
			# run the R tool
			print "\nCalling R:\n---\n$Rcall\n---\n...\n" ;#if $opt->{verbose};
			my $Routput = `$Rcall`;
			print $Routput, "\n" if $opt->{verbose} > 1;
			#      if (not $Routput =~ /Done/m) {
			#        warn $Routput;
			#        exit;
			#      }
		} ## end else [ if ( $opt->{mock} ) ]
	} ## end RUN_R:

	# post process output
	# create the folder if it is not already present
	$opt->{destdir} = rel2abs( $opt->{destdir} );
	mkdir( $opt->{destdir} ) unless -d $opt->{destdir};
	opendir( my $tmpdir, $opt->{tmpdir} );

	#  my @Files = grep { /^$temp_basename/ } readdir ($tmpdir);
	foreach my $file ( readdir($tmpdir) ) {
		next if -d $file;
		my $target = $file;
		my $rename = $opt->{basename};
		if ( $target =~ s/$temp_basename/$rename/ ) {
			move( catfile( $opt->{tmpdir}, $file ), catfile( $opt->{destdir}, $target ) );
			printf "move %s to %s\n", $file, $target if $opt->{verbose};
			#if ( $target =~ /\.csv$/ ) {
			#	csv2xlsx( catfile( $opt->{destdir}, $target ) );
			#}
		}
	} ## end foreach my $file ( readdir($tmpdir...))
} ## end foreach my $config_file (@INI)
print "\nDone.\n" if $opt->{verbose};
exit;

# convert csv to xlsx with same basename (.csv > .xlsx)
sub csv2xlsx {
	my $csv  = shift;
	my $tab  = Data::Table::fromFile($csv);
	my $xlsx = $csv;
	$xlsx =~ s/\.csv$/.xlsx/;

	# colors: even = 22 (silver), odd = 9 (white), head = 30 (dark blue)
	tables2xlsx( $xlsx, [$tab], ["copy-number"], [[22, 9, 30]] );
} ## end sub csv2xlsx

sub random_chars {
	my $n = shift || 3;
	$n = 3  if $n < 1;
	$n = 32 if $n > 32;
	my @chars = ( 1 .. 7, 'a' .. 'z', 'A' .. 'Z' );
	my $out = '';
	for ( 1 .. $n ) {
		$out .= $chars[rand(@chars)];
	}
	return $out;
} ## end sub random_chars

sub process_options {

	# set some defaults if these are mising
	$opt->{global}{title}   = "untitled"  unless $opt->{global}{title};
	$opt->{sample}{name}    = "sample"    unless $opt->{sample}{name};
	$opt->{reference}{name} = "reference" unless $opt->{reference}{name};

	#	if (not defined $opt->{genome}) {
	#	$opt->{genome}{data} = '';
	#	$opt->{genome}{name} = '';
	#	}
	$opt->{genome}{data} = '' if not defined $opt->{genome}{data};
	if ( not defined $opt->{genome}{name} ) {
		$opt->{genome}{name} = '';
		if ( $opt->{genome}{data} =~ /fai$/ ) {
			my @parts = split( /\//, $opt->{genome}{data} );
			my @name = split(/\./, $parts[-1]);
			$opt->{genome}{name} = $name[0] if $name[0];
		}
	} ## end if ( not defined $opt->{genome...})
	#
	_calc_xy();

	# preset clustering parameters
	my $preset = {
		#  above    => 100000,                                                    # { 0}
		#  annot    => 12,                                                        # { 1}
		#  bascol   => 3,                                                         # { 2}
		#  below    => 1000,                                                      # { 3}
		#  chrcol   => 0,                                                         # { 4}
		#  cntcol   => 5,                                                         # { 5}
		#  dbh      => ":memory:",                                                # { 6}
		#  deadband => 1000,                                                      # { 7}
		#  dircol   => 2,                                                         # { 8}
		#  distal   => 250000,                                                    # { 9}
		#  dump     => 1,                                                         # {10}
		#  export   => ["chr", "pos", "dir", "base", "name", "count", "orgname"], # {11}
		#  force    => 200,                                                       # {12}
		#  gaplogic => "OR",                                                      # {13}
		#  input    => "demo_sample_male.txt",                                    # {14}
		#  logic    => "AND",                                                     # {15}
		#  match    => "",                                                        # {16}
		#  min      => 20,                                                        # {17}
		#  namcol   => 4,                                                         # {18}
		names => "/srv/qgen/data/annotation/refGene.txt",                         # {19}

		#  nmbegin  => 6,                                                         # {20}
		#  nmchr    => 2,                                                         # {21}
		#  nmend    => 7,                                                         # {22}
		#  poscol   => 1,                                                         # {23}
		#  relative => 0.25,                                                      # {24}
		#  sep      => "\t",                                                      # {25}
		#  span     => 0.5,                                                       # {26}
		#  suffix   => "",                                                        # {27}
		#  type     => "refgene",                                                 # {28}
	};
	$opt->{cp} = {} if not defined $opt->{cp};
	my $merged = merge( $opt->{cp}, $preset );

	#	dd $merged;
	#	exit;
	$opt->{clustering} = join( " ", $opt->{engine}, map { qq~--$_ "$merged->{$_}"~ } keys %{$merged} );
	
#	if ($opt->{debug}) {
#		dd $opt;
#		exit;
#	}
	
} ## end sub process_options

sub _calc_xy {
	for my $who ( 'sample', 'reference' ) {
		if ( not $opt->{$who}{x} ) {

			# no parameter explicitly set
			if ( is_male( $opt->{$who}{sex} ) ) {
				$opt->{$who}{x} = 1;
			}
			else {
				$opt->{$who}{x} = 2;
			}
		} ## end if ( not $opt->{$who}{x} )
		if ( not $opt->{$who}{y} ) {

			# no parameter explicitly set
			if ( is_male( $opt->{$who}{sex} ) ) {
				$opt->{$who}{y} = 1;
			}
			else {
				$opt->{$who}{y} = 0;
			}
		} ## end if ( not $opt->{$who}{y} )
	} ## end for my $who ( 'sample', 'reference')
} ## end sub _calc_xy

sub is_male {
	my $input = shift;
	return 0 unless $input;
	return 1 if $input =~ /^m/io;
	return 1 if $input =~ /y/io;
	return 0;
} ## end sub is_male

# example .ini file
__DATA__
; set some global parameters or settings documented elsewhere
[global] 
title=untitled

; details on the sample
[sample]
data=filename
name=string
sex=male

; details on the reference
[reference]
data=filename
name=string
sex=male



