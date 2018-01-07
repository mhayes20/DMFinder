#!/usr/bin/perl

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
my $dm_complete_count = 0;
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
    $line =~ s/\R//g;
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
    } else {
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

sub read_amplicon_list
{
  my $path_cn_file = shift;
  open( CN, "<" . $path_cn_file)
    or die("Could not open copy number amplification file!\n");
  my $line;
  my @amplicon_list;
  while ( $line = <CN> ) {
    if ($line =~ /(\w+)\t(\w+)\t(\w+)/) {
      push @amplicon_list, "$1\t$2\t$3";
    }
  }
  return @amplicon_list;
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
#    $key_val[0] =~ s/\R//g;
#    $key_val[1] =~ s/\R//g;
    $vcf_info_dictionary{ uc $key_val[0] } = $key_val[1];
  }
  return %vcf_info_dictionary;
}

sub parse_vcf_line
{
  my $line = shift;
  chomp($line);
  my @temp_rec = split( /\t/, $line );
  my $vcf_info_line       = $temp_rec[7];
  my %vcf_info_dictionary = parse_vcf_record_info($vcf_info_line);
  my %edge;
  $edge{"start_chr"} = lc $temp_rec[0];
  $edge{"start_pos"} = int($temp_rec[1]);
  $edge{"end_chr"}   = lc $vcf_info_dictionary{"CHR2"};
  $edge{"end_pos"}   = int( $vcf_info_dictionary{"END"} );
  return %edge; 
}

sub add_edge_if_exist
{
  my $amplicon_list_ref  = shift;
  my $other_chr          = shift;
  my $other_pos          = shift;
  my $i                  = shift;
  my $window             = shift;
  my @amplicon_list = @$amplicon_list_ref;
  for ( my $j = $i + 1 ; $j < scalar(@amplicon_list) ; $j++ ) {
    my $temp_line = $amplicon_list[$j];
    chomp($temp_line);
    my @temp_rec_other = split( /\t/, $temp_line );
    if ( $other_chr eq $temp_rec_other[0]
      && abs( $other_pos - $temp_rec_other[1] ) <= $window )
    {
      $g->add_edge( $amplicon_list[$i], $amplicon_list[$j] );
      $explored{ $amplicon_list[$i] } = 0;
      $explored{ $amplicon_list[$j] } = 0;
      my $edge = $amplicon_list[$i] . " " . $amplicon_list[$j];
      $explored_edges{ $edge } = 0;
      $viz->add_edge_once( $amplicon_list[$i], $amplicon_list[$j] );
      
    } elsif ( $other_chr eq $temp_rec_other[0]
           && abs( $other_pos - $temp_rec_other[2] ) <= $window )
    {
      $g->add_edge( $amplicon_list[$i], $amplicon_list[$j] );
      $explored{ $amplicon_list[$i] } = 0;
      $explored{ $amplicon_list[$j] } = 0;
      my $edge = $amplicon_list[$i] . " " . $amplicon_list[$j];
      $explored_edges{ $edge } = 0;
      $viz->add_edge_once( $amplicon_list[$i], $amplicon_list[$j] );
    }
  }
}    

sub split_amplicon {
  my $amplicon_ref = shift;
  my $mid_pos      = shift;
  my $verbose      = shift;

  if ($verbose) {
    print join '', 
      "splitting the amplicon ", 
      $amplicon_ref->{"chr"}, ":", $amplicon_ref->{"start"}, 
      "-", $amplicon_ref->{"end"},  "mid=$mid_pos\n";
  }
  my $amp1 = join "\t", $amplicon_ref->{"chr"}, $amplicon_ref->{"start"}, $mid_pos;
  my $amp2 = join "\t", $amplicon_ref->{"chr"}, $mid_pos, $amplicon_ref->{"end"};
  my @splitted_amplicons;
  push @splitted_amplicons, $amp1;
  push @splitted_amplicons, $amp2;
  return @splitted_amplicons;
}

sub print_discarded_double_minute {
  my $shortest_path_ref = shift;
  my $average_cov_ref   = shift; 
  my $dm_type           = shift;
  my @shortest_path     = @$shortest_path_ref;
  my @average_cov       = @$average_cov_ref;
  print "Following possible $dm_type double minute has been discarded due to ";
  print "high variance of average mapping coverage among the amplicons:\n";
  my $shortest_path_str = join '->', @shortest_path;
  print "$shortest_path_str\n";
  print "The averave mapping coverages of the amplicons are:\n@average_cov\n\n";
}


if ( scalar(@ARGV) != 11 ) {
  die("Usage: perl program.pl [SV FILE] [CN SEGMENT FILE] [WINDOW SIZE] 
                              [BAM FILE] [MINQUAL] [MIN CYCLIC] [MIN NON CYCLIC]
                              [REPORT FILE] [GRAPH FILE] [SPLIT AMPLICONS] [VERBOSITY]\n");
}
open( SV, "<" . $ARGV[0] )
  or die("Could not open structural variant breakpoint file!\n");
my $line      = "";
my $temp_line = "";
my $path_to_cn_file = $ARGV[1];
my $window          = int($ARGV[2]);
my $bam_file        = $ARGV[3];
my $min_qual        = $ARGV[4];
my $min_cyclic      = $ARGV[5];
my $min_non_cyclic  = $ARGV[6];
my $report_file     = $ARGV[7];
my $graph_file      = $ARGV[8];
my $split_amplicons = $ARGV[9];
my $verbose         = $ARGV[10];
my $e               = "";
my @amplicon_list   = (); # will store all amplified segment records
my $startFlag       =  0; # 1 if graph constructor should look at start
                          # of a CN segment, 0 otherwise.

@amplicon_list = read_amplicon_list($path_to_cn_file);

my @SV = ();

my $vcf_comment_line_pattern = "^#.*";
@SV = <SV>;
if ($split_amplicons) {
  foreach my $vcf_line ( @SV ) {
    if ( $vcf_line =~ /$vcf_comment_line_pattern/ ) {
      next;
    }
    my %edge = parse_vcf_line($vcf_line);
    foreach my $amplicon ( @amplicon_list ) {
      my @seg_rec = split( /\t/, $amplicon );
      my %amplicon;
      $amplicon{"chr"}   = $seg_rec[0];
      $amplicon{"start"} = int($seg_rec[1]);
      $amplicon{"end"}   = int($seg_rec[2]);

     if ( $edge{"start_chr"} eq $amplicon{"chr"} 
         && $edge{"start_pos"} > $amplicon{"start"} + 3 * $window
         && $edge{"start_pos"} < $amplicon{"end"}   - 3 * $window
      ) 
      {
        my( $index )= grep { $amplicon_list[$_] eq $amplicon } 0..$#amplicon_list;
        splice @amplicon_list, $index, 1;
        push @amplicon_list, split_amplicon(\%amplicon, $edge{'start_pos'});
      }
      elsif ( $edge{"end_chr"} eq $amplicon{"chr"}
          && $edge{"end_pos"} > $amplicon{"start"} + 3 * $window
          && $edge{"end_pos"} < $amplicon{"end"}   - 3 * $window  
      )
      {
        my( $index )= grep { $amplicon_list[$_] eq $amplicon } 0..$#amplicon_list;
        splice @amplicon_list, $index, 1;
        push @amplicon_list, split_amplicon(\%amplicon, $edge{'end_pos'}, $verbose);
      }
    }
  }
  print "Amplicons splitting finished.\n";
}
for ( my $i = 0 ; $i < scalar(@amplicon_list) - 1 ; $i++ ) {
  my $seg_line = $amplicon_list[$i];
  if ( $seg_line =~ /$vcf_comment_line_pattern/ ) {
    next;
  }
  chomp($seg_line);
  my @seg_rec = split( /\t/, $seg_line );
  my %amplicon;
  $amplicon{"chr"}   = $seg_rec[0];
  $amplicon{"start"} = int($seg_rec[1]);
  $amplicon{"end"}   = int($seg_rec[2]);
  for ( my $k = 0 ; $k < scalar(@SV) ; $k++ ) {
    #Search SV file
    $line = $SV[$k];
    if ( $line =~ /$vcf_comment_line_pattern/ ) {
      next;
    }
    my %edge = parse_vcf_line($line);
    if ( $edge{"start_chr"} eq $amplicon{"chr"} )
    {
      if ( abs( $edge{"start_pos"} - $amplicon{"start"} ) <= $window
        || abs( $edge{"start_pos"} - $amplicon{"end"} )   <= $window )
      {
        add_edge_if_exist( \@amplicon_list, $edge{"end_chr"}, 
                           $edge{"end_pos"}, $i, $window);
      }
    }
    elsif ( $edge{"end_chr"} eq $amplicon{"chr"} )
    {
      if ( abs( $edge{"end_pos"} - $amplicon{"start"} ) <= $window 
        || abs( $edge{"end_pos"} - $amplicon{"end"} ) <= $window )
      {
        add_edge_if_exist( \@amplicon_list, $edge{"start_chr"}, 
                           $edge{"start_pos"}, $i, $window);
      }
    } 
  }
}


my @V = $g->vertices;
#print "amplicon_list: ".$V[0]."\n";
foreach my $e (@V) {
  dfs($e);
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
        if ($apsp->path_length( $m, $n ) == 1) {
          next;
        }
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
  if (scalar(@shortest_path) == 0) {
    next;
  }
  my @average_cov = ();
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
    my $amplicon_average_cov = average( "$tmp1_zeros", $i );
    if ( $amplicon_average_cov >= $max_cov ) {
      $max_cov = $amplicon_average_cov;
    }
    if ( $amplicon_average_cov <= $min_cov ) {
      $min_cov = $amplicon_average_cov;
    }
    push @average_cov, int($amplicon_average_cov);
  }
  if ( $verbose && $max_cov - $min_cov >= 35) {
    print_discarded_double_minute(\@shortest_path, \@average_cov, "complete");
  }
  if ( scalar(@shortest_path) > $min_cyclic  #Enforce minimum length requirement
    && ( $max_cov - $min_cov ) < 35
     )  
  {
    write_dm_segment_to_csv_file( *OUTPUT_CVS, 
                                  ++$dm_index, 
                                  \@shortest_path,
                                  'complete' );
#print OUTPUT_CVS "SEGMENTS FOR PREDICTED DOUBLE MINUTE: @shortest_path\n*************************************\n";
#print REPORT "SEGMENTS FOR PREDICTED DOUBLE MINUTE: @shortest_path\n*************************************\n";
#print "SCC is @$e\n";
#print "size of SCC: ".scalar(@$e)."\n";
    $dm_count++;
    $dm_complete_count++;
  }
  $min_cov = 100000000;
  $max_cov = 0;
  my $ssize = scalar(@shortest_path);
  $h->delete_vertices(@shortest_path);
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
#print "Looking for WCC...\n";
my @wcc = $h->weakly_connected_components();
foreach my $e (@wcc) {
  my @shortest_path = @$e;
  my @average_cov = ();
  for ( my $i = 0 ; $i < scalar(@shortest_path) ; $i++ ) {
    chomp( $shortest_path[$i] );
    my @rec_i = split( /\t/, $shortest_path[$i] );
    my $chr   = $rec_i[0];
    my $start = $rec_i[1];
    my $end   = $rec_i[2];
    system( "samtools depth -r $chr:$start-$end -Q $min_qual $bam_file > $tmp1" );
    system( "cat $tmp1 | awk 'BEGIN { prev_chr=\"\";prev_pos=0;} { if(\$1==prev_chr && prev_pos+1!=int(\$2)) {for(i=prev_pos+1;i<int(\$2);++i) {printf(\"%s\\t%d\\t0\\n\",\$1,i);}} print; prev_chr=\$1;prev_pos=int(\$2);}' > $tmp1_zeros" );
    my $amplicon_average_cov = average( "$tmp1_zeros", $i );
    if ( $amplicon_average_cov >= $max_cov ) {
      $max_cov = $amplicon_average_cov;
    }
    if ( $amplicon_average_cov <= $min_cov ) {
      $min_cov = $amplicon_average_cov;
    }
    push @average_cov, int($amplicon_average_cov);
  }
  if ( $verbose && $max_cov - $min_cov >= 35) {
    print_discarded_double_minute(\@shortest_path, \@average_cov, 'incomplete');
  }
 
#print "SCC is @$e\n";
#print "size of shortest_path: ".scalar(@shortest_path)."\n";
#print "max_cov=$max_cov, min_cov=$min_cov\n";
  if ( scalar(@shortest_path) >= $min_non_cyclic
    && ( $max_cov - $min_cov ) < 35
    ) #Once again, this is to ensure we don't get a predicted dmin that has only one amplicon
  {
    write_dm_segment_to_csv_file( *OUTPUT_CVS, 
                                  ++$dm_index,
                                  \@shortest_path,
                                 'incomplete' );
#print "SEGMENTS FOR PREDICTED INCOMPLETE DOUBLE MINUTE: @shortest_path\n*************************************\n";
#print REPORT "SEGMENTS FOR PREDICTED INCOMPLETE DOUBLE MINUTE: @shortest_path\n*************************************\n";
#print "SCC is @$e\n";
    $h->delete_vertices(@shortest_path);
#print "size of SCC: ".scalar(@$e)."\n";
    $dm_count++;
  }
  $min_cov = 100000000;
  $max_cov = 0;
}
C:
cleanup();
print "There are $dm_count predicted double minutes, including $dm_complete_count perfectly identified.\n";
#close REPORT;
close OUTPUT_CVS;
close SV;
