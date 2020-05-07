# Jack Taylor and Kirk Phillips
# BIO 371
# Semester Project Group 5

use strict;
use warnings;

# Ask the user for the filename of our DNA sequence
print "Please enter the filename of the FASTA format DNA sequence.\n";
my $filename = <STDIN>;
chomp $filename;

# Open the file, or kill the program upon failure
open(FH, '<', $filename) or die $!;

# Remove FASTA header and store the remaining text into a string
my $fastaHeader = <FH>;
my $seq = "";
while (<FH>) {
    $seq .= <FH>;
    chomp $seq;
}
# Close the file
close(FH);

# Convert the DNA sequence into a protein-coding sequence
$seq =~ tr/T/U/;

# Print this sequence to the user so that they can verify
print ("\nThe protein coding sequence is: $seq\n");

# Now that the sequence is ready, we parse the codon usage table
# Ask the user for the filename of our codon usage table
print "Please enter the filename of the Codon Usage Table of a reference genome.\n";
$filename = <STDIN>;
chomp $filename;

# Open the file, or kill the program upon failure
open(FH, '<', $filename) or die $!;

# Extract only the triplets and the frequencies from the table
my @words;
my @codons;
my @freqs;
my $counter = 0;
while (<FH>) {
    chomp;
    @words = split(' ');
    $counter = 0;
    foreach my $word (@words) {
        if ($counter % 2 == 0) {
            push(@codons, $word);
        }
        else {
            push(@freqs, $word);
        }
        $counter++;
    }
}
# Close the file
close(FH);

foreach my $freq (@freqs) {
    # Uses regex for finding text matching format: (some number) and deleting it
    $freq =~ s/\(.*\)//;
}

# Hash the codons to their frequencies
my %table;
@table{@codons} = @freqs;

# Print the hash to the user so they can verify
print "\nThe resulting hash from the table is as shown:\n";
while (my ($k, $v) = each %table) {
    print "$k => $v\n";
}

# Create an hash of all codon to amino acid pairs
# All synonymous codons share the same number
my %codings = ("UUU", 1, "UUC", 1, #phe
               "UUA", 2, "UUG", 2, "CUU", 2, "CUC", 2, "CUA", 2, "CUG", 2, #leu
               "UCU", 3, "UCC", 3, "UCA", 3, "UCG", 3, "AGU", 3, "AGC", 3, #ser
               "UAU", 4, "UAC", 4, #tyr
               "UGU", 5, "UGC", 5, #cys
               "UGG", 6, #trp
               "CCU", 7, "CCC", 7, "CCA", 7, "CCG", 7, #pro
               "CAU", 8, "CAC", 8, #his
               "CAA", 9, "CAG", 9, #gin
               "CGU", 10, "CGC", 10, "CGA", 10, "CGG", 10, "AGA", 10, "AGG", 10, #arg
               "AUU", 11, "AUC", 11, "AUA", 11, #lle
               "AUG", 12, #met
               "ACU", 13, "ACC", 13, "ACA", 13, "ACG", 13, #thr
               "AAU", 14, "AAC", 14, #asn
               "AAA", 15, "AAG", 15, #lys
               "GUU", 16, "GUC", 16, "GUA", 16, "GUG", 16, #val
               "GCU", 17, "GCC", 17, "GCA", 17, "GCG", 17, #ala
               "GAU", 18, "GAC", 18, #asp
               "GAA", 19, "GAG", 19, #glu
               "GGU", 20, "GGC", 20, "GGA", 20, "GGG", 20, #gly
               "UAA", 21, "UAG", 21, "UGA", 21); #stop codons

# Now we can begin to calculate the CAI
print "\nNow calculating the CAI...\n";
# First step: get observed frequencies for each codon
# Split sequence into groups of 3 nucleotides for translation
my @seqCodons = ($seq =~ m/.{3}/g);
# Create variable for total codons in sequence
my $total = @seqCodons;

# Count the number of times each codon appears in the sequence
my %count;
$count{$_}++ foreach @seqCodons;

# Store these observed counts into arrays
my @observedCodons = ();
my @observedCounts = ();
while (my ($key, $value) = each(%count)) {
    push(@observedCodons, $key);
    push(@observedCounts, $value);
}
# Hash the codons to their observed counts
my %observedCounts;
@observedCounts{@observedCodons} = @observedCounts;

# Now, iterate through each codon to calculate a value for Wi
my $index = 0;
my $fi = -1;
my $fj = -1;
my $wi = 0;
my $AAfound = -1;
my $synTotal = 0.0;
my $maxSyn = -1;

foreach my $seqCodon (@seqCodons) {
    $index = 0;
    $fi = -1;
    $fj = -1;
    
    # First we need to find the associated AA code for our current codon
    $AAfound = -1;
    keys %codings; # reset the iterator so previous calls don't affect it
    while(my($k, $v) = each %codings) {
        if ($seqCodon eq $k) {
            $AAfound = $v;
        }
    }
    
    if ($AAfound < 0) {
        die "ERROR: codon frequency not found in table.\n";
    }
    
    
    # Find the observed count of the current codon and divide it by the total count
    # of all synonymous codons
    $synTotal = 0.0;
    keys %codings; # reset the iterator so previous calls don't affect it
    while(my($k, $v) = each %codings) {
        if ($AAfound == $v) {
            # Make sure this codon exists in our sequence before adding its value
            if (exists($observedCounts{$k})) {
                $synTotal = $synTotal + $observedCounts{$k};
            }
        }
    }
    
    # Checks to see if calculation of fi was done correctly
    if ($synTotal <= 0) {
        die "ERROR: fi could not be calculated.\n";
    }
    
    $fi = ($observedCounts{$seqCodon} / $synTotal);
    
    # Now we must calculate the max fj frequency (max freq for all synonymous codons in table)
    $synTotal = 0.0;
    $maxSyn = -1;
    keys %codings; # reset the iterator so previous calls don't affect it
    while(my($k, $v) = each %codings) {
        if ($AAfound == $v) {
            $synTotal = $synTotal + $table{$k};
            if ($table{$k} > $maxSyn) {
                $maxSyn = $table{$k};
            }
        }
    }
    
    # Checks to see if calculation of fj was done correctly
    if ($synTotal <= 0) {
        die "ERROR: fj could not be calculated.\n";
    }
    
    $fj = $maxSyn / $synTotal;
    
    # Use logarithm form to avoid underflow
    $wi = $wi + log($fi / $fj);
}

# Now that the summation is complete, finish the CAI calculation
my $CAI = (1 / $total) * $wi;
$CAI = exp $CAI;

print "\nThe calculated Codon Adaptation Index is: $CAI\n";
