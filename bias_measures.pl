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
print ("The protein coding sequence is:$seq\n");
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
my $count = 0;
while (<FH>) {
    chomp;
    @words = split(' ');
    $count = 0;
    foreach my $word (@words) {
        if ($count % 2 == 0) {
            push(@codons, $word);
        }
        else {
            push(@freqs, $word);
        }
        $count++;
    }
}
# Close the file
close(FH);

foreach my $freq (@freqs) {
    # Uses regex for finding text matching format: (some number) and deleting it
    $freq =~ s/\(.*\)//;
}

print "The respective codons and their frequencies (per thousand):\n";
for (my $i = 0; $i < 64; $i++) {
    print "($codons[$i], $freqs[$i])\n";
}


