#!/usr/bin/perl

#Merge adjancent amplicons in RDxplorer output if 1) they are less than 10,000 bp apart and 2) if the segment medians are with 2

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
my $sv_file = 0;    # no longer used
my $window  = 0;    # no longer used
#my $sv_file = $ARGV[1];
#my $window = $ARGV[2];
open( IN, "<$file" ) or die("Could not find the input file!\n");
#open(SV, "<$sv_file") or die("Could not open the SV link file\n");
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
#my @SV = <SV>;
sub checkLinks {
#This function sees if the two amplicons under consideration are proximal to any SV links.
#NOTE: no longer used
  my $chr     = $_[0];
  my $amp1    = $_[1];
  my $amp2    = $_[2];
  my $SVlinks = $_[3];
  my $line;
  my @rec = ();
  for ( my $i = 0 ; $i < scalar(@$SVlinks) ; $i++ ) {
    $line = $$SVlinks[$i];
    chomp($line);
    @rec = split( /\t/, $line );
    #print "CHR is $chr and rec[0] is ".$rec[0]."\n";
    if ( ( $chr eq $rec[0] && abs( $amp1 - $rec[1] ) <= $window )
      || ( $chr eq $rec[3] && abs( $amp1 - $rec[4] ) <= $window )
      || ( $chr eq $rec[0] && abs( $amp2 - $rec[1] ) <= $window )
      || ( $chr eq $rec[3] && abs( $amp2 - $rec[4] ) <= $window ) )
    {
      #print "IN IF\n";
      # print "CHR is $chr and rec[0] is ".$rec[0]."\n";
      return 1;
    }
  }
  return 0;
}

#while($line = <IN>)
for ( my $i = 0 ; $i < scalar(@filerec) - 1 ; $i++ ) {
  #$temp_pos = tell(IN);
A:
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
    $new_start = $record[8]
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
#if((abs($record[9] - $record2[8]) < 10000 && abs($record[4] - $record2[4]) <= 3))
      if ( ( abs( $record[9] - $record2[8] ) < 10000 && $chr eq $record2[7] ) )
      {
   #Do nothing. We really want the point where this condition is no longer true.
      }
      else {
        $new_end = $record[9];
        $merge_count++;
        $i = $j + 1
          ; #Set i to the place where the new end position was created. If we didn't do this, then the code would just
            #start over from the original position.

        print "$chr_prefix$chr\t$new_start\t$new_end\n";
        goto A;
      }
    }
  } else {
    print "$chr_prefix$chr\t" . $record[8] . "\t" . $record[9] . "\n";
  }
  #}
}

print "$merge_count merge operations occurred\n";
close IN;
#close SV;
1;

