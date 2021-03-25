#!/usr/local/bin/perl -w 
# -*- perl -*-
######################################################################
# gffcmp.pl -- Simple program to compare BioWare Gff files
# Copyright (c) 2004 Tero Kivinen
# All Rights Reserved.
######################################################################
#         Program: gffcmp.pl
#	  $Source: /u/samba/nwn/bin/RCS/gffcmp.pl,v $
#	  Author : $Author: kivinen $
#
#	  (C) Tero Kivinen 2004 <kivinen@iki.fi>
#
#	  Creation          : 20:10 Aug  6 2004 kivinen
#	  Last Modification : 01:26 May 24 2007 kivinen
#	  Last check in     : $Date: 2007/05/23 22:26:22 $
#	  Revision number   : $Revision: 1.7 $
#	  State             : $State: Exp $
#	  Version	    : 1.96
#	  Edit time	    : 23 min
#
#	  Description       : Simple program to compare BioWare Gff files
#
#	  $Log: gffcmp.pl,v $
#	  Revision 1.7  2007/05/23 22:26:22  kivinen
#	  	No changes.
#
#	  Revision 1.6  2007/05/23 22:26:12  kivinen
#	  	Fixed path splitting to accept windows paths.
#
#	  Revision 1.5  2005/02/05 18:09:47  kivinen
#	  	Fixed =from to =item.
#
#	  Revision 1.4  2005/02/05 17:50:26  kivinen
#	  	Added documentation.
#
#	  Revision 1.3  2005/02/05 14:34:15  kivinen
#	  	Added --diff option.
#
#	  Revision 1.2  2004/09/20 11:46:41  kivinen
#	  	Added internal globbing.
#
#	  Revision 1.1  2004/08/15 12:36:24  kivinen
#	  	Created.
#
#	  $EndLog$
#
#
#
######################################################################
# initialization

require 5.6.0;
package GffCmp;
use strict;
use Getopt::Long;
use File::Glob ':glob';
use Gff;
use GffRead;
use Pod::Usage;

$Opt::verbose = 0;
$Opt::exclude = undef;
$Opt::include = undef;
$Opt::exclude_field = undef;
$Opt::include_field = undef;
$Opt::diff = 0;

######################################################################
# Get version information

open(PROGRAM, "<$0") || die "Cannot open myself from $0 : $!";
undef $/;
$Prog::program = <PROGRAM>;
$/ = "\n";
close(PROGRAM);
if ($Prog::program =~ /\$revision:\s*([\d.]*)\s*\$/i) {
    $Prog::revision = $1;
} else {
    $Prog::revision = "?.?";
}

if ($Prog::program =~ /version\s*:\s*([\d.]*\.)*([\d]*)\s/mi) {
    $Prog::save_version = $2;
} else {
    $Prog::save_version = "??";
}

if ($Prog::program =~ /edit\s*time\s*:\s*([\d]*)\s*min\s*$/mi) {
    $Prog::edit_time = $1;
} else {
    $Prog::edit_time = "??";
}

$Prog::version = "$Prog::revision.$Prog::save_version.$Prog::edit_time";
$Prog::progname = $0;
$Prog::progname =~ s/^.*[\/\\]//g;

$| = 1;

######################################################################
# Read rc-file

if (defined($ENV{'HOME'})) {
    read_rc_file("$ENV{'HOME'}/.gffcmprc");
}

######################################################################
# Option handling

Getopt::Long::Configure("no_ignore_case");

if ($0 =~ /diff/) {
    $Opt::diff = 1;
}

if (!GetOptions("config=s" => \$Opt::config,
		"verbose|v+" => \$Opt::verbose,
		"help|h" => \$Opt::help,
		"exclude|e=s" => \$Opt::exclude,
		"include|i=s" => \$Opt::include,
		"exclude-field=s" => \$Opt::exclude_field,
		"include-field=s" => \$Opt::include_field,
		"diff|d" => sub { $Opt::diff = 1; },
		"cmp|c" => sub { $Opt::diff = 0; },
		"version|V" => \$Opt::version) || defined($Opt::help)) {
    usage();
}

if (defined($Opt::version)) {
    print("\u$Prog::progname version $Prog::version by Tero Kivinen.\n");
    exit(0);
}

while (defined($Opt::config)) {
    my($tmp);
    $tmp = $Opt::config;
    undef $Opt::config;
    if (-f $tmp) {
	read_rc_file($tmp);
    } else {
	die "Config file $Opt::config not found: $!";
    }
}

######################################################################
# Main loop

$| = 1;

my($i, %args, $gff, $exitval);

%args = (include => $Opt::include,
	 exclude => $Opt::exclude,
	 include_field => $Opt::include_field,
	 exclude_field => $Opt::exclude_field);

$args{'filename'} = shift;

if (join(";", @ARGV) =~ /[*?]/) {
    my(@argv);
    foreach $i (@ARGV) {
	push(@argv, bsd_glob($i));
    }
    @ARGV = @argv;
}

if ($Opt::verbose) {
    print("Reading basefile $args{filename}...\n");
}
$gff = GffRead::read(%args);

$exitval = 0;
foreach $i (@ARGV) {
    my($gff1, $ret);

    if ($Opt::verbose) {
	print("Reading file $i...\n");
    }
    $args{'filename'} = $i;
    $gff1 = GffRead::read(%args);
    if ($Opt::diff) {
	my(@ret);
	@ret = $gff->diff($gff1);
	$ret = join("\n", @ret);
	if ($ret eq '') {
	    $ret = undef;
	} else {
	    $ret = "\n" . $ret;
	}
    } else {
	$ret = $gff->match($gff1);
    }
    if ($ret) {
	print("File $i do not match to base file : $ret\n");
	$exitval = 1;
	exit $exitval if (!$Opt::diff);
    } else {
	if ($Opt::verbose) {
	    print("File $i matches...\n");
	}
    }
}

exit $exitval;

######################################################################
# Read rc file

sub read_rc_file {
    my($file) = @_;
    my($next, $space);
    
    if (open(RCFILE, "<$file")) {
	while (<RCFILE>) {
	    chomp;
	    while (/\\$/) {
		$space = 0;
		if (/\s+\\$/) {
		    $space = 1;
		}
		s/\s*\\$//g;
		$next = <RCFILE>;
		chomp $next;
		if ($next =~ s/^\s+//g) {
		    $space = 1;
		}
		if ($space) {
		    $_ .= " " . $next;
		} else {
		    $_ .= $next;
		}
	    }
	    if (/^\s*([a-zA-Z0-9_]+)\s*$/) {
		eval('$Opt::' . lc($1) . ' = 1;');
	    } elsif (/^\s*([a-zA-Z0-9_]+)\s*=\s*\"([^\"]*)\"\s*$/) {
		my($key, $value) = ($1, $2);
		$value =~ s/\\n/\n/g;
		$value =~ s/\\t/\t/g;
		eval('$Opt::' . lc($key) . ' = $value;');
	    } elsif (/^\s*([a-zA-Z0-9_]+)\s*=\s*(.*)\s*$/) {
		my($key, $value) = ($1, $2);
		$value =~ s/\\n/\n/g;
		$value =~ s/\\t/\t/g;
		eval('$Opt::' . lc($key) . ' = $value;');
	    }
	}
	close(RCFILE);
    }
}


######################################################################
# Usage

sub usage {
    Pod::Usage::pod2usage(0);
}

=head1 NAME

gffcmp - Compare gff files

gffdiff - Print diff of gff files

=head1 SYNOPSIS

gffcmp [B<--help>|B<-h>] [B<--version>|B<-V>] [B<--verbose>|B<-v>]
    [B<--config> I<config-file>]
    [B<--exclude>|B<-e> I<exclude-regexp>]
    [B<--include>|B<-i> I<include-regexp>]
    [B<--exclude-field> I<exclude-regexp>]
    [B<--include-field> I<include-regexp>]
    [B<--diff>|B<-d>]
    [B<--cmp>|B<-c>]
    I<filename> ...

gffdiff [B<--help>|B<-h>] [B<--version>|B<-V>] [B<--verbose>|B<-v>]
    [B<--config> I<config-file>]
    [B<--exclude>|B<-e> I<exclude-regexp>]
    [B<--include>|B<-i> I<include-regexp>]
    [B<--exclude-field> I<exclude-regexp>]
    [B<--include-field> I<include-regexp>]
    [B<--diff>|B<-d>]
    [B<--cmp>|B<-c>]
    I<filename> ...

gffcmp B<--help>

gffdiff B<--help>

=head1 DESCRIPTION

B<gffcmp> checks if two or more gffs are identical and if so returns
0. In case the files are different it prints out the first difference
and returns 1.

B<gffdiff> checks if two or more gffs are identical and if so returns
0. In case the files are different it prints out differences between
files and returns 1. All files compared against the first file. 

=head1 OPTIONS

=over 4

=item B<--help> B<-h>

Prints out the usage information.

=item B<--version> B<-V>

Prints out the version information. 

=item B<--verbose> B<-v>

Enables the verbose prints. This option can be given multiple times,
and each time it enables more verbose prints. 

=item B<--config> I<config-file>

All options given by the command line can also be given in the
configuration file. This option is used to read another configuration
file in addition to the default configuration file. 

=item B<--exclude> B<-e> I<exclude-regexp>

Exclude the given regexp when reading the data in. This will skip the
whole structure behind the given structure, meaning that B<--include>
cannot be used to get parts of that back. This can be used to speed up
the processing if only specific parts of the tree is required.
Normally this should be something like I<^/Creature List> meaning that
all creature list information is skipped when reading gff.

=item B<--include> B<-i> I<include-regexp>

Only include the given regexp when reading the data in. This will skip
all other structures which do not match the regexp. This can be used
to speed up the processing if only specific parts of the tree is
required. Normally this should be something like I<^/Creature List>
meaning that only  creature list information is read in. 

=item B<--exclude-field> I<exclude-regexp>

Exclude given fields to be read in in case their labels match the
given regexp. This only matches the end labels, not intermediate
structure labels. 

=item B<--include-field> I<include-regexp>

Only include given fields matching the given regexp to the structures.
This only matches the end labels, not intermediate structure labels.

=item B<--diff> B<-d>

Print multiple differences, i.e. do not stop on first difference. 

=item B<--cmp> B<-c>

Print only first difference, and exit.

=back

=head1 EXAMPLES

    gffcmp file1.git file2.git
    gffcmp --exclude-field TemplateResRef file1.uti file2.uti
    gffdiff ztk_ada_spear.uti ztk_ada_scythe.uti

=head1 FILES

=over 6

=item ~/.gffcmprc

Default configuration file.

=back

=head1 SEE ALSO

gffencode(1), gffmodify(1), gffprint(1), Gff(3), and GffRead(3).

=head1 AUTHOR

Tero Kivinen <kivinen@iki.fi>.

=head1 HISTORY

This program first apperead as B<gffcmp>, but the B<gffdiff> support
was added later to print out more than difference between files. 

