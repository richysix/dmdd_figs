#!/usr/bin/env perl

# PODNAME: get_number_sig_genes.pl
# ABSTRACT: get the number of DE genes from DESeq2 results

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;

# get options
my %options;
get_and_check_options();

while(<>){
    chomp;
    my $gene = $_;
    
    my $gene_count;
    my $type = 'unfiltered';
    my $deseq_file =
        File::Spec->catfile($options{'base_dir'}, $gene,
            'deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers',
            $options{'comparison'} . '.sig.tsv' );
    
    if( -e $deseq_file ) {
        # this assumes that the file has a header line and
        # all the rest are significant genes
        my $lines = `/usr/bin/wc -l $deseq_file`;
        chomp $lines;
        my ( $num_genes, undef ) = split / /, $lines;
        warn $num_genes, "\n" if $options{'debug'};
        $gene_count = --$num_genes;
    }
    # print output
    $gene_count = defined $gene_count ? $gene_count : 'NA';
    print join("\t", $gene, $options{'comparison'}, $type, $gene_count, ), "\n";
    
    $gene_count = undef;
    $type = 'filtered';
    my $baseline_comparison_dir =
        File::Spec->catfile($options{'base_dir'}, $gene,
            'deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers',
            $options{'comparison'} . '.baseline-comp' );
    if( -e $baseline_comparison_dir ) {
        # open baseline overlaps file
        my $overlaps_file = File::Spec->catfile($baseline_comparison_dir,
                                                'overlaps.txt');
        open my $overlap_fh, '<', $overlaps_file;
        while( my $line = <$overlap_fh> ){
            my ( $category, $num_genes, ) = split /: /, $line;
            if( $category eq 'plus_baseline.plus_stage' ||
                   $category eq 'no_baseline.plus_baseline.plus_stage' ){
                $gene_count += $num_genes;
            }
        }
    }
    # print out line
    $gene_count = defined $gene_count ? $gene_count : 'NA';
    print join("\t", $gene, $options{'comparison'}, $type, $gene_count, ), "\n";
}


#
#Is there a baseline comparison directory?
#
#YES -> open baseline-comp.tsv
#	count the number of no_baseline.plus_baseline.plus_stage and plus_baseline.plus_stage sig genes
#
#NO -> does $comparison.sig.tsv exist
#	YES -> open $comparison.sig.tsv
#		count the number of sig genes
#	NO -> output null record
#


################################################################################
# SUBROUTINES

#subroutine_name
#
#  Usage       : subroutine_name( arguments )
#  Purpose     : 
#  Returns     : 
#  Parameters  : 
#  Throws      : 
#  Comments    : None

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
        'base_dir=s',
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
    $options{'base_dir'} = $options{'base_dir'} ? $options{'base_dir'}
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
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

arguments

=back

=head1 OPTIONS

**Same for optional arguments.

=over

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