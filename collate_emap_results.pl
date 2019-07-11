#!/usr/bin/env perl

# PODNAME: collate_emap_results.pl
# ABSTRACT: Description

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use Readonly;

# get options
my %options;
get_and_check_options();

# read in config file
my $config_file = $ARGV[0];
open my $c_fh, '<', $config_file;
if( $options{'header'} ){
    my $line = <$c_fh>;
}
# for each line, check file exists, output results to all file
my $header;
while( my $line = <$c_fh> ){
    chomp $line;
    my ( $mutant, $comparison, $results_set, $file ) = split /\t/, $line;
    if( !-e $file ){
        # check error file
        $file =~ s/results.tsv \z/emap.e/xms;
        open my $err_fh, '<', $file;
        my $skipped = 0;
        while( my $err_line = <$err_fh> ){
            if( $err_line =~ m/\A NO\sSIGNIFICANT\sGENES/xms ) {
                $skipped = 1;
            }
        }
        if( !$skipped ) {
            my $err_msg = "File does not exist: $file\n";
            pod2usage($err_msg);
        }
    } else {
        open my $results_fh, '<', $file;
        if( $options{'header'} ){
            my $results_line = <$results_fh>;
            if (!$header) {
                chomp $results_line;
                my @cols = split /\t/, $results_line;
                print join("\t", qw{Gene Comparison Set}, @cols, "-log10(pvalue)"), "\n";
                $header = 1;
            }
        }
        # for each line, check file exists, output results to all file
        while( my $results_line = <$results_fh> ){
            chomp $results_line;
            my @info = split /\t/, $results_line;
            # This line assumes that the pvalue is in column 7
            Readonly my $P_VALUE_COLUMN => 6;
            Readonly my $BASE_LOG => 10; # to calculate log10 we need to do log(PVALUE)/log(10)
            print join("\t", $mutant, $comparison, $results_set, @info,
                       -log($info[$P_VALUE_COLUMN])/log($BASE_LOG) ), "\n";
        }
    }
}


################################################################################
# SUBROUTINES

#get_and_check_options
#
#  Usage       : get_and_check_options()
#  Purpose     : parse the options supplied to the script using GetOpt::Long
#  Returns     : None
#  Parameters  : None
#  Throws      : 
#  Comments    : The default option are
#                help which print a SYNOPSIS
#                man which prints the full man page
#                debug
#                verbose

sub get_and_check_options {
    
    GetOptions(
        \%options,
        'header',
        'help',
        'man',
        'debug+',
        'verbose',
    ) or pod2usage(2);
    
    # Documentation
    if( $options{help} ) {
        pod2usage( -verbose => 0, -exitval => 1, );
    }
    elsif( $options{man} ) {
        pod2usage( -verbose => 2 );
    }
    
    $options{debug} = $options{debug} ? $options{debug} : 0;
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

collate_emap_results.pl

=head1 DESCRIPTION

Description

=cut

=head1 SYNOPSIS

    collate_emap_results.pl [options] config file | STDIN
        --header                input file have a header line
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

=item config_file

A tab-separated file detailing the results files to use.
Columns are:

=over

=item Mutant - Name of the KO line

=item Comparison - DESeq2 comparison (e.g. hom_vs_het_wt)

=item Set - name of results set

=item File - path of results file

=back

=back

=head1 OPTIONS

=over

=item B<--header>

Input files contain a header line

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=head1 DEPENDENCIES

None

=head1 AUTHOR

=over 4

=item *

Richard White <richard.white@sanger.ac.uk>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2018 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut