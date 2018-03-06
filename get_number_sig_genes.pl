#!/usr/bin/env perl

# PODNAME: get_number_sig_genes.pl
# ABSTRACT: get the number of DE genes from DESeq2 results

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use Readonly;

# get options
my %options;
get_and_check_options();

# set directory to find DESeq results
Readonly my $RESULTS_DIR =>
    'deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers';

while(<>){
    chomp;
    my $gene = $_;
    
    my $gene_count;
    my $type = 'unfiltered';
    my $deseq_file =
        File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
                            $options{'comparison'} . '.sig.tsv' );
    
    if( -e $deseq_file ) {
        # this assumes that the file has a header line and
        # all the rest are significant genes
        my $lines = `wc -l $deseq_file`;
        chomp $lines;
        my ( $num_genes, undef ) = split / /, $lines;
        warn $num_genes, "\n" if $options{'debug'};
        $gene_count = --$num_genes;
    }
    # print output
    $gene_count = defined $gene_count ? $gene_count : 'NA';
    print join("\t", $gene, $options{'comparison'}, 'ko_response', $type, $gene_count, ), "\n";
    
    $type = 'filtered';
    my $baseline_comparison_dir =
        File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
            $options{'comparison'} . '.baseline-comp' );
    $gene_count = undef;
    if( -e $baseline_comparison_dir ) {
        # open baseline overlaps file
        my $overlaps_file = File::Spec->catfile($baseline_comparison_dir,
                                                'overlaps.txt');
        open my $overlap_fh, '<', $overlaps_file;
        while( my $line = <$overlap_fh> ){
            chomp $line;
            my ( $category, $num_genes, ) = split /: /, $line;
            # print out line
            if ($category ne 'ko_response') {
                print join("\t", $gene, $options{'comparison'}, $category, $type, $num_genes, ), "\n";
            } else {
                $gene_count = $num_genes;
            }
        }
    }
    $gene_count = defined $gene_count ? $gene_count : 'NA';
    print join("\t", $gene, $options{'comparison'}, 'ko_response', $type, $gene_count, ), "\n";
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
        'dir=s',
        'comparison=s',
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
    
    # defaults
    $options{'debug'} = $options{'debug'} ? $options{'debug'} : 0;
    $options{'dir'} = $options{'dir'} ? $options{'dir'}
        :   './';
    $options{'comparison'} = $options{'comparison'} ? $options{'comparison'}
        :   'hom_vs_het_wt';
    
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

get_number_sig_genes.pl

=head1 DESCRIPTION

Takes a list of directory names on stdin and counts up the genes called as
differentially expressed for that DESeq2 run.

=cut

=head1 SYNOPSIS

    get_number_sig_genes.pl [options] input file | STDIN
        --dir                   working directory   default: cwd
        --comparison            comparison to use   default: hom_vs_het_wt
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

=item <input_file>

A list of directory names to get the results from. Can also be supplied on STDIN.
The DESeq results are expected to be in

    [working_directory]/[gene]/[RESULTS_DIR]/[comparison].sig.tsv

    see --dir, --comparison options and RESULTS_DIR (set at top of the script)
    
=back

=head1 OPTIONS

=over

=item B<--dir>

Working directory to search for the directory names supplied in the input.
default: current working directory

=item B<--comparison>

DESeq comparison to get the results from.
default: hom_vs_het_wt

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

This software is Copyright (c) 2017 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut