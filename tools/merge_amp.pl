#!/usr/bin/perl

use warnings;
use strict;
use Scalar::Util qw(looks_like_number);

if ( scalar(@ARGV) != 1 and scalar(@ARGV) != 2 ) {
  die("Usage is perl program.pl [RDX SUM FILT FILE] [OPTIONAL CHROMOSOME PREFIX]\n");
}

my $file       = $ARGV[0];
my $chr_prefix = $ARGV[1];
if ( !defined($chr_prefix) ) {
  $chr_prefix = "";
}
open( IN, "<$file" ) or die("Could not find the input file!\n");

my @record  = ();
my @record2 = ();
my $line    = "";
my $line2   = "";
my $temp_pos;
my $new_start;
my $new_end;
my $line3;
my $line4;
my $chr;
my $merge_count = 0;
my @filerec = <IN>;
my $length = scalar(@filerec);
for ( my $i = 0 ; $i < (scalar(@filerec) - 1) ; $i++ ) {
  $line = $filerec[$i];
  chomp($line);
  @record = split( /\t/, $line );
  $chr   = $record[7];
  $line2 = $filerec[ $i + 1 ];
  chomp($line2);
  @record2 = split( /\t/, $line2 );
  if ( !looks_like_number( $record[9] ) || !looks_like_number( $record2[8] ) ) {
    next;
  }
  if ( ( abs( $record[9] - $record2[8] ) < 10000 && $chr eq $record2[7] ) ) {
    $new_start = int($record[8])
      ; #Start of merged amplicon will be the start of the first original amplicon
    for ( my $j = $i + 1 ; $j < scalar(@filerec) - 1 ; $j++ ) {
      $line3 = $filerec[$j];
      chomp($line3);
      $line4 = $filerec[ $j + 1 ];
      chomp($line4);
      @record  = split( /\t/, $line3 );
      @record2 = split( /\t/, $line4 );
      if ( !looks_like_number( $record[9] )
        || !looks_like_number( $record2[8] ) )
      {
        next;
      }
      if ( abs( int($record[9]) - int($record2[8]) ) < 10000 && $chr eq $record2[7] )
      {
        #Do nothing. We really want the point where this condition is no longer true.
      }
      else {
        $new_end = int($record[9]);
        $merge_count++;
        $i = $j + 1;
        # Set i to the place where the new end position was created. 
        # If we didn't do this, then the code would just
        # start over from the original position.

        print "$chr_prefix$chr\t$new_start\t$new_end\n";
        last;
      }
    }
  } else {
    print "$chr_prefix$chr\t" . int($record[8]) . "\t" . int($record[9]) . "\n";
  }
}

print STDERR  "$merge_count merge operations occurred\n";
close IN;
;

