#!/usr/bin/env perl

# PODNAME: get_go_results.pl
# ABSTRACT: Description

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use Readonly;
#use JSON;
#use REST::Client;
use Data::Dumper;

# get options
my %options;
my $check_baseline;
get_and_check_options();

# set directory to find DESeq results
Readonly my $RESULTS_DIR => defined $options{'results_dir'}
    ?   $options{'results_dir'}
    :   'deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers';
Readonly my $FISHER_THRESHOLD => 100;

# print header
print join("\t", qw{Gene Comparison Domain GO.ID Term Annotated Significant
                    Expected pval -log10(pvalue) Fold.Enrichment log10(FE) }, ), "\n";

my @comparisons = (qw{hom_vs_het_wt het_vs_wt hom_vs_het});
my @baseline_sets = (qw{mrna_abnormal mrna_as_wt not_used ko_response});

# open baseline comparison file
open my $baseline_fh, '>', $options{'baseline_file'};
print {$baseline_fh}
    join("\t", qw{Gene Comparison Set Domain GO.ID Term Annotated Significant
                    Expected pval }, ), "\n";

# also keep a track of the terms (and desc) and output the most common ones
my %times_seen;
my %description_for;
# go through each directory name supplied
while(<>){
    chomp;
    my $gene = $_;
    
    foreach my $comparison ( @comparisons ) {
        
        if( $check_baseline ) {
            # Check for baseline analysis. Add later
            my $baseline_dir = 
                File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
                    $comparison . '.baseline-comp' );
            if( -e $baseline_dir ){
                foreach my $set ( @baseline_sets ){
                    my $go_results_file =
                        File::Spec->catfile($baseline_dir, $set . '.fisher',
                            $options{'go_domain'} . '.sig.tsv', );
                    
                    my $num_sig_terms = 0;
                    if( -e $go_results_file ){
                        open my $go_fh, '<', $go_results_file;
                        while( my $line = <$go_fh> ){
                            next if $line =~ m/\A GO.ID/xms;
                            $num_sig_terms++;
                            my ($go_term, $description, $annotated, $significant,
                                $expected, $pvalue, undef ) = split /\t/, $line;
                            # filter for observed gene number and fold enrichment
                            my ( $filter, $fe ) = filter_term( $significant, $expected, \%options, ); 
                            next if $filter;
                            
                            if ($set eq 'ko_response'){
                                $times_seen{$go_term}++;
                                $description_for{$go_term} = $description;
                            }
                            
                            # calculate -log10 pvalue and log10 fold enrichment
                            # convert '< 1e-30' to number
                            $pvalue = $pvalue eq '< 1e-30' ? 1e-30 : $pvalue;
                            my $log_p = -log($pvalue) / log(10);
                            my $log_fe = log($fe) / log(10);
                            
                            print {$baseline_fh}
                                join("\t", $gene, $comparison, $set, $options{'go_domain'},
                                           $go_term, $description, $annotated, $significant,
                                           $expected, $pvalue, $log_p, $fe, $log_fe, ), "\n";
                            if( $set eq 'ko_response' ){
                                print join("\t", $gene, $comparison, $options{'go_domain'},
                                           $go_term, $description, $annotated, $significant,
                                           $expected, $pvalue, $log_p, $fe, $log_fe, ), "\n";
                            }
                        }
                        close $go_fh;
                    }
                    if( $num_sig_terms == 0 && $set eq 'ko_response' ){
                        warn join(q{ }, 'GENE:', $gene, 'COMPARISON:', $comparison,
                                  '- NO GO ENRICHMENTS'), "\n";
                    }
                }
                # if we're using the baseline results don't get the unfiltered results
                last; # go to next comparison
            }
        }
        
        # Check existence of sig file
        my $deseq_file =
            File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
                $comparison . '.sig.tsv' );
        if( -e $deseq_file ) {
            # Get number of sig genes
            my $lines = `wc -l $deseq_file`;
            chomp $lines;
            my ( $num_genes, undef ) = split / /, $lines;
            $num_genes--;
            if( $num_genes == 0 ) {
                warn join(q{ }, 'GENE:', $gene, 'COMPARISON:', $comparison,
                          '- NO SIG DE GENES'), "\n";
                last;
            }
            # If there are less than 100 DE genes use fisher
            my $go_results_file;
            if( $num_genes < $FISHER_THRESHOLD ) {
                $go_results_file =
                    File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
                        $comparison . '.fisher',
                        $options{'go_domain'} . '.sig.tsv', );
            } else {
                # else use KS
                $go_results_file =
                    File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
                        $comparison . '.go',
                        $options{'go_domain'} . '.sig.tsv', );
            }
            # check for existence of results file
            if( !-e $go_results_file ) {
                warn join(q{ }, 'GENE:', $gene, 'COMPARISON:', $comparison,
                          '- NO GO FILE'), "\n";
                last;
            }
            
            # open results file
            open my $go_fh, '<', $go_results_file;
            my $num_sig_terms = 0;
            while( my $line = <$go_fh> ){
                next if $line =~ m/\A GO.ID/xms;
                $num_sig_terms++;
                my ($go_term, $description, $annotated, $significant,
                    $expected, $pvalue, undef ) = split /\t/, $line;
                # filter for observed gene number and fold enrichment
                my ( $filter, $fe ) = filter_term( $significant, $expected, \%options, ); 
                next if $filter;
                
                $times_seen{$go_term}++;
                $description_for{$go_term} = $description;
                
                # calculate -log10 pvalue and log10 fold enrichment
                # convert '< 1e-30' to number
                $pvalue = $pvalue eq '< 1e-30' ? 1e-30 : $pvalue;
                my $log_p = -log($pvalue) / log(10);
                my $log_fe = log($fe) / log(10);
                
                print join("\t", $gene, $comparison, $options{'go_domain'},
                           $go_term, $description, $annotated, $significant,
                           $expected, $pvalue, $log_p, $fe, $log_fe, ), "\n";
            }
            close $go_fh;
            if( $num_sig_terms == 0 ){
                warn join(q{ }, 'GENE:', $gene, 'COMPARISON:', $comparison,
                          '- NO GO ENRICHMENTS'), "\n";
                last;
            }
            last;
        }
    }
}

close $baseline_fh;

# output most common terms
open my $shared_fh, '>', $options{'shared_file'};
# sort terms by frequency and pick top 40 terms
my $num_terms;
my $freq;

foreach my $term ( sort { $times_seen{$b} <=> $times_seen{$a} } keys %times_seen ){
    $num_terms++;
    if( $num_terms > $options{'number_top_terms'} ){
        # check for ties
        if( $times_seen{$term} < $freq ) {
            last;
        } else {
            print {$shared_fh} join("\t", $term, $description_for{$term},
                                    $times_seen{$term}, ), "\n";
            $freq = $times_seen{$term};
        }
    } else {
        print {$shared_fh} join("\t", $term, $description_for{$term},
                                $times_seen{$term}, ), "\n";
        $freq = $times_seen{$term};
    }
}



################################################################################
# SUBROUTINES

# filter_term
#
#  Usage       : filter_term($significant, $expected, $options,)
#  Purpose     : Decied whether to filter a term
#  Returns     : filter 0 or 1
#              : fold enrichment Num
#  Parameters  : significant Int
#                expected Num
#                options HashRef
#  Throws      : 
#  Comments    : None

sub filter_term {
    my ( $significant, $expected, $options, ) = @_;
    
    my $filter = 0;
    # filter for observed gene number and fold enrichment
    my $fe = $expected == 0
        ?   $significant/($expected + 0.01)
        :   $significant/$expected;
    if( $fe <= $options->{'fe_threshold'} ) {
        $filter = 1;
    }
    if( $significant <= $options->{'observed_threshold'} ) {
        $filter = 1;
    }
    return( $filter, $fe );
}

# get_and_check_options
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
        'baseline_file=s',
        'shared_file=s',
        'number_top_terms=i',
        'results_dir=s',
        'observed_threshold=i',
        'fe_threshold=f',
        'go_domain=s',
        'skip_baseline',
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
    Readonly my @GO_DOMAINS => ('BP', 'MF', 'CC', );
    Readonly my %GO_DOMAINS => map { $_ => 1 } @GO_DOMAINS;
    if( !defined $options{'go_domain'} ) {
        $options{'go_domain'} = $GO_DOMAINS[0];
    } elsif( !exists $GO_DOMAINS{$options{'go_domain'}} ) {
        my $msg = "\nOption --go_domain must be one of " . join(q{, }, @GO_DOMAINS ) . "\n";
        pod2usage($msg);
    }

    # defaults
    $options{'debug'} = $options{'debug'} ? $options{'debug'} : 0;
    $options{'baseline_file'} =
        $options{'baseline_file'} ? $options{'baseline_file'}
        :   'topgo-baseline-results.tsv';
    $options{'shared_file'} =
        $options{'shared_file'} ? $options{'shared_file'}
        :   'topgo-shared-terms.tsv';
    
    $options{'number_top_terms'} = $options{'number_top_terms'} ? $options{'number_top_terms'} : 40;
    $options{'observed_threshold'} = $options{'observed_threshold'} ? $options{'observed_threshold'} : 2;
    $options{'fe_threshold'} = $options{'fe_threshold'} ? $options{'fe_threshold'} : 1;
    
    # turn skip_baseline option into check baseline
    $check_baseline = $options{'skip_baseline'} ? 0 : 1;
    
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

get_go_results.pl

=head1 DESCRIPTION

Get GO enrichment result for multiple analyses.

=cut

=head1 SYNOPSIS

    get_go_results.pl [options] input file | STDIN
        --dir                   working directory                       default: cwd
        --baseline_file         name of file for baseline results       default: topgo-baseline-results.tsv
        --shared_file           name of file for shared terms           default: topgo-shared-terms.tsv
        --number_top_terms      number of top shared terms to output    default: 40
        --results_dir           name of results directory
        --observed_threshold    minimum number of significant genes
        --fe_threshold          minimum fold enrichment
        --go_domain             GO domain (BP, MF, CC)                  default: BP
        --skip_baseline         switch to skip looking for baseline results
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

=item <input_file>

A list of directory names to get the results from. Can also be supplied on STDIN.
The DESeq results are expected to be in

    [working_directory]/[gene]/[results_dir]/[comparison].sig.tsv

    see --dir and --results_dir options
    comparisons used are hom_vs_het_wt, het_vs_wt, hom_vs_het in that order.
    Results are taken from the first results file that exists.

=back

=head1 OPTIONS

=over

=item B<--dir>

Working directory to search for the directory names supplied in the input.
default: current working directory

=item B<--baseline_file>

Name of file to output baseline results to. 
default: topgo-baseline-results.tsv

=item B<--baseline_file>

Name of file to output the top 40 most shared terms across the experiments. 
default: topgo-shared-terms.tsv

=item B<--results_dir>

Name of results directory inside each directory in the list
    default: 

=item B<--go_domain>

Gene Ontology domain to use (One of BP, MF or CC).
default: BP

=item B<--skip_baseline>

Skip looking for baseline results

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