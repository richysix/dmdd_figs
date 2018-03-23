#!/usr/bin/env perl

# PODNAME: blind_go_results.pl
# ABSTRACT: Description

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use Readonly;
use File::Spec;
use File::Path qw(make_path);
use File::Copy::Recursive qw(dircopy);

# get options
my %options;
get_and_check_options();

# set directory to find GO results in
Readonly my $RESULTS_DIR => defined $options{'results_dir'}
    ?   $options{'results_dir'}
    :   'deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers';
Readonly my $FISHER_THRESHOLD => 100;

# set up output dir
if ( !-e $options{'output_dir'} ){
    if ( !make_path($options{'output_dir'}) ) {
        die("Could not create directory, $options{'output_dir'}\n");
    }
}

my $id_for; #  %id_for->{ $gene }{ filtered => x, unfiltered => y }
my $ids; # %ids->{ id => 1 }

# input is a list of mutant names
while(<>){
    chomp;
    my $gene = $_;
    
    ## Filtered set
    my $type = 'filtered';
    # Check for baseline analysis.
    my $go_results_dir = 
        File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
            $options{'comparison'} . '.baseline-comp', 'ko_response.fisher' );
    
    ( $ids, $id_for ) = copy_dir( $go_results_dir, $ids, $id_for,
                                    $gene, $type, );
    
    
    ## Unfiltered set
    $type = 'unfiltered';
    
    # Check existence of sig file
    my $deseq_file =
        File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
            $options{'comparison'} . '.sig.tsv' );
    if( !-e $deseq_file ) {
        die("Results file, $deseq_file, does not exist\n");
    }
    
    # Get number of sig genes
    my $lines = `wc -l $deseq_file`;
    chomp $lines;
    my ( $num_genes, undef ) = split / /, $lines;
    $num_genes--;
    if( $num_genes == 0 ) {
        warn join(q{ }, 'GENE:', $gene, 'COMPARISON:', $options{'comparison'},
                  '- NO SIG DE GENES'), "\n";
        next;
    }
    # If there are less than 100 DE genes use fisher
    if( $num_genes < $FISHER_THRESHOLD ) {
        $go_results_dir =
            File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
                $options{'comparison'} . '.fisher', );
    } else {
        # else use KS
        $go_results_dir =
            File::Spec->catfile($options{'dir'}, $gene, $RESULTS_DIR,
                $options{'comparison'} . '.go', );
    }
    
    ( $ids, $id_for ) = copy_dir( $go_results_dir, $ids, $id_for,
                                    $gene, $type, );
    
    
}

# output a file of ids to gene and filtering type
my @set1;
my @set2;
my $key_file = $options{'key_file'}
        ? $options{'key_file'}
        : File::Spec->catfile($options{'output_dir'}, 'go_blind_key.tsv');
open my $key_fh, '>', $key_file;
foreach my $gene ( sort keys %{$id_for} ) {
    print $key_fh join("\n",
                       map { join("\t", $gene, $_, $id_for->{$gene}{$_} ); }
                            (qw{filtered unfiltered}) ), "\n";
    
    # pick one to assign to set1, assign other to set2
    my $coin_toss = int(rand(2));
    if($coin_toss){
        push @set1, $id_for->{$gene}{'unfiltered'};
        push @set2, $id_for->{$gene}{'filtered'};
    } else {
        push @set1, $id_for->{$gene}{'filtered'};
        push @set2, $id_for->{$gene}{'unfiltered'};
    }
}

print "Set1: ", join(q{ }, @set1, ), "\n";
print "Set2: ", join(q{ }, @set2, ), "\n";

################################################################################
# SUBROUTINES

# copy_dir
#
#  Usage       : copy_dir( $go_results_dir, $ids, $id_for, $gene, $type, )
#  Purpose     : copy the supplied directory to a new randomly numbered directory
#  Returns     : $ids, Hashref lookup table for already used ids
#                $id_for, Hashref mapping of gene and filtering type to id
#  Parameters  : $go_results_dir, Str Name of the directory to copy
#                $ids, Hashref lookup table for already used ids
#                $id_for, Hashref mapping of gene and filtering type to id
#                $gene, Str Gene name
#                $type, Str 'filtered' or 'unfiltered'
#  Throws      : If $go_results_dir does not exist
#                If output directory could not be created
#                If $go_results_dir could not be copied
#  Comments    : None

sub copy_dir {
    my ( $go_results_dir, $ids, $id_for, $gene, $type, ) = @_;
    # check GO results dir exists
    if ( !-e $go_results_dir ) {
        die("Results directory, $go_results_dir, does not exist!\n");
    }
    
    # pick a random number and make sure it hasn't been used yet
    (my $id, $ids, $id_for) = choose_id($ids, $id_for, $gene, $type, );
    
    # make a new directory with the id number
    my $dir = File::Spec->catfile($options{'output_dir'}, $id);
    if ( !make_path($dir) ) {
        die("Could not create directory, $dir\n");
    }
    
    $File::Copy::Recursive::MaxDepth = 1;
    # copy the results directory into the new directory
    my $num_of_files = dircopy($go_results_dir, $dir);
    if( !$num_of_files ){
        die("Could not copy directory, $go_results_dir\n");
    }
    
    return( $ids, $id_for );
}

# choose_id
#
#  Usage       : choose_id( $ids, $id_for, $gene, $type, )
#  Purpose     : copy the supplied directory to a new randomly numbered directory
#  Returns     : $ids, Hashref lookup table for already used ids
#                $id_for, Hashref mapping of gene and filtering type to id
#  Parameters  : $ids, Hashref lookup table for already used ids
#                $id_for, Hashref mapping of gene and filtering type to id
#                $gene, Str Gene name
#                $type, Str 'filtered' or 'unfiltered'
#  Throws      : 
#  Comments    : None

sub choose_id {
    my ($ids, $id_for, $gene, $type, ) = @_;
    
    my $id = 0;
    while(!$id){
        $id = int(rand(1000));
        $id = exists $ids->{$id} ? 0 : $id;
    }
    # add to hashes
    $ids->{$id} = 1;
    $id_for->{$gene}{$type} = $id;
    
    return($id, $ids, $id_for);
}

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
        'results_dir=s',
        'output_dir=s',
        'comparison=s',
        'key_file=s',
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
    
    # set defaults
    $options{'dir'} = $options{'dir'} ? $options{'dir'} : './';
    $options{'output_dir'} = $options{'output_dir'}
        ? $options{'output_dir'}
        : 'go_blind';
    $options{'comparison'} = $options{'comparison'}
        ? $options{'comparison'}
        : 'hom_vs_het_wt';
    
    $options{debug} = $options{debug} ? $options{debug} : 0;
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

blind_go_results.pl

=head1 DESCRIPTION

Description

=cut

=head1 SYNOPSIS

    blind_go_results.pl [options] input file | STDIN
        --dir                   working directory                       default: cwd
        --results_dir           name of results directory
        --output_dir            name of directory to copy results to    default: go_blind
        --comparison            name of DESeq2 comparison               default: hom_vs_het_wt
        --key_file              name of file that contains the mapping of ids to results set default: [output_dir]/go_blind_key.tsv
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

=item B<input_file>

A list of gene names to get the results from. Can also be supplied on STDIN.
The GO results are expected to be in

    [working_directory]/[gene]/[results_dir]/[comparison].sig.tsv

    see --dir, --comparison and --results_dir options


=back

=head1 OPTIONS

=over

=item B<--dir>

Working directory to search for the directory names supplied in the input.
default: current working directory

=item B<--results_dir>

Name of results directory inside each directory in the list
    default: deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers

=item B<--output_dir>

name of directory to copy results to
    default: go_blind
    
=item B<--comparison>

Name of the DESeq2 comparison. The GO results files contain the comparison name.
    default: hom_vs_het_wt

=item B<--key_file>

Name of file that contains the mapping of ids to results set.
Columns are Gene, Filtering type, id
    
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