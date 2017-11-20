#!/usr/bin/perl

use lib("/nfs/gpfs/PDS0257/cwr0353/SAM_REPOSITORY/SIM/TEST_TEMP/BIOINFORMATICS/NEW_DIRECTORY/DM_FIND/Graph-0.94/lib");
#Will change to "use Graph;"
use lib("/nfs/gpfs/PDS0257/cwr0353/SAM_REPOSITORY/SIM/TEST_TEMP/BIOINFORMATICS/NEW_DIRECTORY/DM_FIND/Graph-Easy-0.64/lib");

use Graph::Easy;
use Graph::Directed;
use Graph::Undirected;
use warnings;
use strict;
use List::Util qw (min max);

my $min_cov  = 100000000;
my $max_cov  = 0;
my $g        = Graph::Undirected->new;
my $h        = Graph::Directed->new;
my $viz      = Graph::Easy->new( { 'undirected', 'true' } );
my $dviz     = Graph::Easy->new( { 'undirected', 'false' } );
my $dm_count = 0;
if ( $viz->is_directed() ) {
  die;
}
use constant INFINITY => 99999999999999999999999;
my @visited_edges    = ();
my @visited_vertices = ();
my %explored       = ();    #explored vertices
my %explored_edges = ();
my @edges          = ();
my $count_edges_added = 0;
my $file_count        = 0;
my $tmp1              = "temp1.$$.txt";
my $tmp2              = "temp2.$$.txt";
my $tmp1_sampled      = "temp1.$$.sampled.txt";
my $tmp2_sampled      = "temp2.$$.sampled.txt";
my $tmp1_zeros        = "temp1.$$.zeros.txt";
#open(REPORT, ">report.txt") or die("Could not create the report file!\n");
sub sampleFiles {
  my $freq = $_[0];
  open( F1, "<$tmp1" ) or die("Could not find the first coverage file!\n");
  open( F2, "<$tmp2" ) or die("Could not find the second coverage file!\n");
  open( S1, ">$tmp1_sampled" )
    or die("Could not write to the sampled file!\n");
  open( S2, ">$tmp2_sampled" )
    or die("Could not write to the second samnpled file!\n");
  my $count = 0;
  my $line1;
  my $line2;
  while ( $line1 = <F1> ) {
    if ( $count % $freq == 0 ) {
      #print S1 "$line1";
    }
    $count++;
  }
  $count = 0;
  while ( $line2 = <F2> ) {
    if ( $count % $freq == 0 ) {
      #print S2 "$line2";
    }
    $count++;
  }
  close F1;
  close F2;
  close S1;
  close S2;
}

sub cleanup {
  system("rm -f $tmp1");
  system("rm -f $tmp2");
  system("rm -f $tmp1_sampled");
  system("rm -f $tmp2_sampled");
  system("rm -f $tmp1_zeros");
}

sub write_dm_segment_to_csv_file {
  my ( $output_csv, $dm_idx, $dm_segments_ref, $dm_type ) = @_;
  my @dm_segments = @{$dm_segments_ref};
  for my $dm_segment (@dm_segments) {
    my ( $chr, $start, $end ) = split( /\t/, $dm_segment );
    print $output_csv "DM$dm_idx,$chr,$start,$end,$dm_type\n";
  }
}

sub average {
  my $file = $_[0];
  my $ID   = $_[1];
  open( F1, "<$file" ) or die("average: Could not open the input file!\n");
#open(F2, ">$ID.txt") or die("average: Could not write to the average file!\n");
  my $count = 0;
  my $n     = 0;
  my $sum;
  my $line;
  my $val;
  my @record = ();
  while ( $line = <F1> ) {
    chomp($line);
    @record = split( /\t/, $line );
    $val = $record[2];
    if ( $count % 100 == 0 ) {
      $sum += $val;
      $n++;
      #print F2 "$sum\n";
    }
    $count++;
  }
  #print "Sample size: $n\n";
  #close F2;
  close F1;
  return $sum / $n;
}

sub dfs {
  my @adj = $g->neighbors( $_[0] );
  push @visited_vertices, $_[0];
  $explored{ $_[0] } = 1;
  for ( my $i = 0 ; $i < scalar(@adj) ; $i++ ) {
    if ( $explored{ $adj[$i] } == 0 ) {
      $h->add_weighted_edge( $_[0], $adj[$i], 1 );
      $dviz->add_edge_once( $_[0], $adj[$i] );
      #print "CASE 1: ADDED DIRECTED EDGE ".$_[0]." ---> ".$adj[$i]."\n";
      push @edges, [ $_[0], $adj[$i], 1 ];
      $count_edges_added++;
      $explored_edges{ $_[0] . " " . $adj[$i] } = 1;
      $explored_edges{ $adj[$i] . " " . $_[0] } = 1;
      dfs( $adj[$i] );
    }
    else {
      if (
        (
             ( not defined $explored_edges{ $adj[$i] . " " . $_[0] } )
          || ( not defined $explored_edges{ $_[0] . " " . $adj[$i] } )
        )
        || ( $explored_edges{ $adj[$i] . " " . $_[0] } == 0
          && $explored_edges{ $_[0] . " " . $adj[$i] } == 0 )
        )
      {
        $h->add_weighted_edge( $_[0], $adj[$i], 1 );
        $dviz->add_edge_once( $_[0], $adj[$i] );

        $explored_edges{ $adj[$i] . " " . $_[0] } = 1;
        push @edges, [ $_[0], $adj[$i], 1 ];
        #print "2nd PUSHED ".$_[0]." ".$adj[$i]."\n";
        $count_edges_added++;
      }
    }
  }
}

sub parse_vcf_record_info {
  my $vcf_info_line = $_[0];
  my %vcf_info_dictionary;
  my @vcf_info_records = split( /\;/, $vcf_info_line );
  foreach my $vcf_info_record (@vcf_info_records) {
    my @key_val = split( /\=/, $vcf_info_record );
    if ( scalar @key_val < 2 ) {
      next;
    }
    $vcf_info_dictionary{ uc $key_val[0] } = $key_val[1];
  }
  return %vcf_info_dictionary;
}
if ( scalar(@ARGV) != 9 ) {
  die(
"Usage is perl program.pl [SV FILE] [CN SEGMENT FILE] [WINDOW SIZE] [BAM FILE] [MINQUAL] [MIN CYCLIC] [MIN NON CYCLIC] [REPORT FILE] [GRAPH FILE]\n"
  );
}
open( SV, "<" . $ARGV[0] )
  or die("Could not open structural variant breakpoint file!\n");
open( CN, "<" . $ARGV[1] )
  or die("Could not open copy number amplification file!\n");
my $line      = "";
my $temp_line = "";
my $window = $ARGV[2];
my $bam_file       = $ARGV[3];
my $min_qual       = $ARGV[4];
my $min_cyclic     = $ARGV[5];
my $min_non_cyclic = $ARGV[6];
my $report_file    = $ARGV[7];
my $graph_file     = $ARGV[8];
my $other_chr      = "";
my $other_pos      = 0;
my $cn_chr         = "";
my $cn_start       = "";
my $cn_end         = "";
my $sv_chr_i    = "";
my $sv_chr_i_bp = "";
my $e           = "";
my $sv_chr_j    = "";
my $sv_chr_j_bp = "";
my @record     = ();
my @seg_record = ();    #will store all amplified segment records
my $startFlag =
  0;  #1 if graph constructor should look at start of a CN segment, 0 otherwise.
my @SV = ();
my @temp_rec       = ();
my @temp_rec_other = ();
@SV = <SV>;
my $seg_line = "";
my @seg_rec  = ();
while ( $line = <CN> ) {
  @record = split( /\t/, $line );
  push @seg_record, $line;
}
my $vcf_comment_line_pattern = "^[^#]*#[^#]*\$";
for ( my $i = 0 ; $i < scalar(@seg_record) - 1 ; $i++ ) {
  $seg_line = $seg_record[$i];
  if ( $seg_line =~ /$vcf_comment_line_pattern/ ) {
    next;
  }
  chomp($seg_line);
  @seg_rec = split( /\t/, $seg_line );
  for ( my $k = 0 ; $k < scalar(@SV) ; $k++ ) {
    #Search SV file
    $line = $SV[$k];
    if ( $line =~ /$vcf_comment_line_pattern/ ) {
      next;
    }
    #print "$line\n";
    chomp($line);
    #print "link: $line\n";
    @temp_rec = split( /\t/, $line );
    #Now get the other chromosome and position
    my $vcf_info_line       = $temp_rec[7];
    my %vcf_info_dictionary = parse_vcf_record_info($vcf_info_line);
    my $chr2     = lc $vcf_info_dictionary{"CHR2"};
    my $chr2_loc = int( $vcf_info_dictionary{"END"} );
    if ( ( $temp_rec[0] eq $seg_rec[0] || $temp_rec[0] eq $seg_rec[0] )
      && abs( $temp_rec[1] - $seg_rec[1] ) <= $window )
    {
      #See if the link goes to another segment
      $other_chr = $chr2;
      $other_pos = $chr2_loc;
      for ( my $j = $i + 1 ; $j < scalar(@seg_record) ; $j++ ) {
        #Check other CN segments for link
        $temp_line = $seg_record[$j];
        chomp($temp_line);
        @temp_rec_other = split( /\t/, $temp_line );
        if (
          (
               $other_chr eq $temp_rec_other[0]
            || $other_chr eq $temp_rec_other[0]
          )
          && abs( $other_pos - $temp_rec_other[1] ) <= $window
          )
        {
          #print "Edge added 1: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
          $g->add_edge( $seg_record[$i], $seg_record[$j] );
          $explored{ $seg_record[$i] } = 0;
          $explored{ $seg_record[$j] } = 0;
          $explored_edges{ $seg_record[$i] . " " . $seg_record[$j] } =
            0;
          $viz->add_edge_once( $seg_record[$i], $seg_record[$j] );
          goto A;
        }
        if (
               ( $other_chr eq $temp_rec_other[0]
              || $other_chr eq $temp_rec_other[0] )
            && abs( $other_pos - $temp_rec_other[2] ) <= $window
          )    #check other end of CN segment
        {
          #print "Edge added 2: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
          $g->add_edge( $seg_record[$i], $seg_record[$j] );
          $explored{ $seg_record[$i] } = 0;
          $explored{ $seg_record[$j] } = 0;
          $explored_edges{ $seg_record[$i] . " " . $seg_record[$j] } =
            0;
          $viz->add_edge_once( $seg_record[$i], $seg_record[$j] );
          goto A;
        }
      }
    }
    elsif ( ( $chr2 eq $seg_rec[0] || $chr2 eq $seg_rec[0] )
      && abs( $chr2_loc - $seg_rec[1] ) <= $window )
    {
      $other_chr = $temp_rec[0];
      $other_pos = $temp_rec[1];
      for ( my $j = $i + 1 ; $j < scalar(@seg_record) ; $j++ ) {
        #Check other CN segments for link
        $temp_line = $seg_record[$j];
        chomp($temp_line);
        @temp_rec_other = split( /\t/, $temp_line );
        if (
          (
               $other_chr eq $temp_rec_other[0]
            || $other_chr eq $temp_rec_other[0]
          )
          && abs( $other_pos - $temp_rec_other[1] ) <= $window
          )
        {
          #print "Edge added 3: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
          $g->add_edge( $seg_record[$i], $seg_record[$j] );

          $explored{ $seg_record[$i] } = 0;
          $explored{ $seg_record[$j] } = 0;
          $explored_edges{ $seg_record[$i] . " " . $seg_record[$j] } =
            0;
          $viz->add_edge_once( $seg_record[$i], $seg_record[$j] );
          goto A;
        }
        if (
             ( $other_chr eq $temp_rec_other[0]
            || $other_chr eq $temp_rec_other[0] )
          && abs( $other_pos - $temp_rec_other[2] ) <= $window)    #check other end of CN segment
        {
          #print "Edge added 4: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
          $g->add_edge( $seg_record[$i], $seg_record[$j] );
          $explored{ $seg_record[$i] } = 0;
          $explored{ $seg_record[$j] } = 0;
          $explored_edges{ $seg_record[$i] . " " . $seg_record[$j] } =
            0;
          $viz->add_edge_once( $seg_record[$i], $seg_record[$j] );
          goto A;
        }

      }
    }
    elsif (
      ( $temp_rec[0] eq $seg_rec[0] || $temp_rec[0] eq $seg_rec[0] )
      && abs( $temp_rec[1] - $seg_rec[2] ) <= $window )
    {
      #Start of CN segment has link
      $other_chr = $chr2;
      $other_pos = $chr2_loc;
      for ( my $j = $i + 1 ; $j < scalar(@seg_record) ; $j++ ) {
        #Check other CN segments for link
        $temp_line = $seg_record[$j];
        chomp($temp_line);
        @temp_rec_other = split( /\t/, $temp_line );
        if (
          (
               $other_chr eq $temp_rec_other[0]
            || $other_chr eq $temp_rec_other[0]
          )
          && abs( $other_pos - $temp_rec_other[1] ) <= $window
          )
        {
       #ADD EDGE
       #SV link at start of current CN segment goes to start of other CN segment
          #print "Edge added 5: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
          $g->add_edge( $seg_record[$i], $seg_record[$j] );
          $explored{ $seg_record[$i] } = 0;
          $explored{ $seg_record[$j] } = 0;
          $explored_edges{ $seg_record[$i] . " " . $seg_record[$j] } = 0;
          $viz->add_edge_once( $seg_record[$i], $seg_record[$j] );
          goto A;
        }
        if (
          (
               $other_chr eq $temp_rec_other[0]
            || $other_chr eq $temp_rec_other[0]
          )
          && abs( $other_pos - $temp_rec_other[2] ) <= $window
          )    #check other end of CN segment
        {
          #print "Edge added 6: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
          $g->add_edge( $seg_record[$i], $seg_record[$j] );
          $explored{ $seg_record[$i] } = 0;
          $explored{ $seg_record[$j] } = 0;
          $explored_edges{ $seg_record[$i] . " " . $seg_record[$j] } =
            0;
          $viz->add_edge_once( $seg_record[$i], $seg_record[$j] );
          goto A;
        }
      }
    }
    elsif ( ( $chr2 eq $seg_rec[0] || $chr2 eq $seg_rec[0] )
      && abs( $chr2_loc - $seg_rec[2] ) <= $window )
    {
      #End of CN segment has link
      #See if the link goes to another segment
      $other_chr = $temp_rec[0];
      $other_pos = $temp_rec[1];
      for ( my $j = $i + 1 ; $j < scalar(@seg_record) ; $j++ ) {
        #Check other CN segments for link
        $temp_line = $seg_record[$j];
        chomp($temp_line);
        @temp_rec_other = split( /\t/, $temp_line );
        if (
          (
               $other_chr eq $temp_rec_other[0]
            || $other_chr eq $temp_rec_other[0]
          )
          && abs( $other_pos - $temp_rec_other[1] ) <= $window
          )
        {
          #print "Edge added 7: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
          $g->add_edge( $seg_record[$i], $seg_record[$j] );
          $explored{ $seg_record[$i] } = 0;
          $explored{ $seg_record[$j] } = 0;
          $explored_edges{ $seg_record[$i] . " " . $seg_record[$j] } =
            0;
          $viz->add_edge_once( $seg_record[$i], $seg_record[$j] );
          goto A;
        }
        if (
          (
               $other_chr eq $temp_rec_other[0]
            || $other_chr eq $temp_rec_other[0]
          )
          && abs( $other_pos - $temp_rec_other[2] ) <= $window
          )    #check other end of CN segment
        {
          #print "Edge added 8: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
          $g->add_edge( $seg_record[$i], $seg_record[$j] );
          $explored{ $seg_record[$i] } = 0;
          $explored{ $seg_record[$j] } = 0;
          $explored_edges{ $seg_record[$i] . " " . $seg_record[$j] } =
            0;
          $viz->add_edge_once( $seg_record[$i], $seg_record[$j] );
          goto A;
        }
      }
    }
  A:
  }
}
my @V = $g->vertices;
#print "seg_record: ".$V[0]."\n";
foreach my $e (@V) {
  dfs($e);
}
my $u = 0;
my $v = 0;
for ( $e = 0 ; $e < scalar(@edges) ; $e++ ) {
  $u = $edges[$e][0];
  $v = $edges[$e][1];
}

for ( $e = 0 ; $e < scalar(@edges) ; $e++ ) {
  $u = $edges[$e][0];
  $v = $edges[$e][1];
}
close CN;
#print "The undirected graph is $g\n";
#print "The directed graph is $h\n";
my @adj = ();
@adj = $g->neighbors( $V[0] );
my $apsp = $h->APSP_Floyd_Warshall();    #All pairs shortest path now calculated
#Now calculate the SCCs
my @scc =
  $h->strongly_connected_components;     #scc contains list of anonymous arrays
#print "number of SCCs is ".scalar(@scc)."\n";
#print "\n";
#print "EDGES PUSHED ON STACK: $count_edges_added\n";
#print @scc."\n";
my @shortest_path;
my $path_total_weight = INFINITY;
my $temp_weight       = INFINITY;
open( OUTPUT_CVS, ">$report_file" )
  or die("Could not create the report file $report_file\n");
print OUTPUT_CVS "DM_index,chromosome,amplicon start,amplicon end,type\n";
my $dm_index = 0;
foreach my $e (@scc)  #Cycle through all SCCs in the graph to find potential DMs
{
  @shortest_path     = ();
  $temp_weight       = INFINITY;
  $path_total_weight = INFINITY;
  foreach my $n (@$e) {
    @adj = $h->neighbors($n);
    #We now have all the neighbors of the vertex n
    foreach my $m (@adj) {
      #for each neighbor of n
      if ( ( defined $h->get_edge_weight( $n, $m ) )
        && $h->get_edge_weight( $n, $m ) ne ''
        && ( defined $apsp->path_length( $m, $n ) )
        && $apsp->path_length( $m, $n ) ne '' )
      {
        #Now return shortest path from m to n
        $temp_weight =
          $h->get_edge_weight( $n, $m ) + $apsp->path_length( $m, $n );
    #print "THE EDGE WEIGHT OF INTEREST IS :".$h->get_edge_weight($n, $m)."\n";
    #print "temp_weight: $temp_weight\npath_total_weight: $path_total_weight\n";
        if ( $temp_weight < $path_total_weight ) {
          $path_total_weight = $temp_weight;
          @shortest_path = $apsp->path_vertices( $m, $n );
          #print "PATH VERTICES ARE @shortest_path\n";
        }
      }
    }
  }
  for ( my $i = 0 ; $i < scalar(@shortest_path) ; $i++ ) {
    chomp( $shortest_path[$i] );
    my @rec_i = split( /\t/, $shortest_path[$i] );
    my $chr   = $rec_i[0];
    my $start = $rec_i[1];
    my $end   = $rec_i[2];
    system(
      "samtools depth -r $chr:$start-$end -Q $min_qual $bam_file > $tmp1" );
    system(
"cat $tmp1 | awk 'BEGIN { prev_chr=\"\";prev_pos=0;} { if(\$1==prev_chr && prev_pos+1!=int(\$2)) {for(i=prev_pos+1;i<int(\$2);++i) {printf(\"%s\\t%d\\t0\\n\",\$1,i);}} print; prev_chr=\$1;prev_pos=int(\$2);}' > $tmp1_zeros"
    );
    #print average("$tmp1_zeros", $i)."\n";
    if ( average( "$tmp1_zeros", $i ) >= $max_cov ) {
      $max_cov = average( "$tmp1_zeros", $i );
    }
    if ( average( "$tmp1_zeros", $i ) <= $min_cov ) {
      $min_cov = average( "$tmp1_zeros", $i );
    }
  }
  if ( scalar(@shortest_path) > $min_cyclic
    && ( $max_cov - $min_cov ) < 35 )    #Enforce minimum length requirement
  {
    write_dm_segment_to_csv_file( *OUTPUT_CVS, ++$dm_index, \@shortest_path,
      'compete' );
#print OUTPUT_CVS "SEGMENTS FOR PREDICTED DOUBLE MINUTE: @shortest_path\n*************************************\n";
#print REPORT "SEGMENTS FOR PREDICTED DOUBLE MINUTE: @shortest_path\n*************************************\n";
    #print "SCC is @$e\n";
    #print "size of SCC: ".scalar(@$e)."\n";
    $dm_count++;
  }
  $min_cov = 100000000;
  $max_cov = 0;
  my $ssize = scalar(@shortest_path);
  $h = $h->delete_path(@shortest_path);
  @shortest_path = ();
  open( VIZ, ">$graph_file" )
    or die("Could not create the graph file $graph_file\n");
  if ( $viz->is_directed() ) {
    #die("GRAPH IS DIRECTED");
  }
  print VIZ $viz->as_graphviz();
  close VIZ;
B:
}
#Get connected subgraphs. These are amplicons connected by SV breakpoints that were not predicted in the first two steps.
##These are likely double minutes.
my @wcc = $h->weakly_connected_components();
foreach my $e (@wcc) {
  @shortest_path = @$e;
  for ( my $i = 0 ; $i < scalar(@shortest_path) ; $i++ ) {
    chomp( $shortest_path[$i] );
    my @rec_i = split( /\t/, $shortest_path[$i] );
    my $chr   = $rec_i[0];
    my $start = $rec_i[1];
    my $end   = $rec_i[2];
    system( "samtools depth -r $chr:$start-$end -Q $min_qual $bam_file > $tmp1" );
    system( "cat $tmp1 | awk 'BEGIN { prev_chr=\"\";prev_pos=0;} { if(\$1==prev_chr && prev_pos+1!=int(\$2)) {for(i=prev_pos+1;i<int(\$2);++i) {printf(\"%s\\t%d\\t0\\n\",\$1,i);}} print; prev_chr=\$1;prev_pos=int(\$2);}' > $tmp1_zeros" );
    #print average("$tmp1_zeros", $i)."\n";
    if ( average( "$tmp1_zeros", $i ) >= $max_cov ) {
      $max_cov = average( "$tmp1_zeros", $i );
    }
    if ( average( "$tmp1_zeros", $i ) <= $min_cov ) {
      $min_cov = average( "$tmp1_zeros", $i );
    }
  }
  if ( scalar(@shortest_path) > $min_non_cyclic
    && ( $max_cov - $min_cov ) < 35
    ) #Once again, this is to ensure we don't get a predicted dmin that has only one amplicon
  {
    write_dm_segment_to_csv_file( *OUTPUT_CVS, ++$dm_index, \@shortest_path,
      'incomplete' );
#print "SEGMENTS FOR PREDICTED INCOMPLETE DOUBLE MINUTE: @shortest_path\n*************************************\n";
#print REPORT "SEGMENTS FOR PREDICTED INCOMPLETE DOUBLE MINUTE: @shortest_path\n*************************************\n";
#print "SCC is @$e\n";
    $h = $h->delete_edges(@shortest_path);
    #print "size of SCC: ".scalar(@$e)."\n";
    $dm_count++;
  }
  $min_cov = 100000000;
  $max_cov = 0;
}
C:
cleanup();
print "There are $dm_count predicted double minutes\n";
#close REPORT;
close OUTPUT_CVS;
close SV;
