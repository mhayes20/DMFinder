#!/usr/bin/perl
#
#THIS VERSION WILL DUMP THE SEGMENT AVERAGES TO A FILE TO BE USED BY R


#9/19/14
	#RELAX THE PROXIMITY REQUIREMENT. If the breakpoint is ANYWHERE within the segment, then add the edge. Not just if the segment boundaries are near a breakpoint.


use lib("/fs/project/PDS0257/cwr0353/SAM_REPOSITORY/SIM/TEST_TEMP/BIOINFORMATICS/NEW_DIRECTORY/DM_FIND/Graph-0.94/lib");
use lib("/fs/project/PDS0257/cwr0353/SAM_REPOSITORY/SIM/TEST_TEMP/BIOINFORMATICS/NEW_DIRECTORY/DM_FIND/Graph-Easy-0.64/lib");

use Graph::Easy;
use Graph::Directed;
use Graph::Undirected;

system("rm graph.dot");

use warnings;
use strict;
use List::Util qw (min max);

#Take as input 1) SV breakpoint results and 2) copy number segments produced by RDxplorer
#

my $min_cov = 100000000;
my $max_cov = 0;
my $g = Graph::Undirected->new;
my $h = Graph::Directed->new;
my $viz = Graph::Easy->new({'undirected','true'});
my $dviz = Graph::Easy->new({'undirected','false'});
my $dm_count = 0;

if($viz->is_directed())
{
	die;
}
use constant INFINITY => 99999999999999999999999;
my @visited_edges = ();
my @visited_vertices = ();

my %explored = (); #explored vertices
my %explored_edges = ();
my @edges = ();

my $count_edges_added = 0;
my $file_count = 0;

#open(REPORT, ">report.txt") or die("Could not create the report file!\n");

sub sampleFiles
{

	my $freq = $_[0];
	open(F1, "<temp1.txt") or die("Could not find the first coverage file!\n");
	open(F2, "<temp2.txt") or die("Could not find the second coverage file!\n");

	open(S1, ">temp1.sampled.txt") or die ("Could not write to the sampled file!\n");
	open(S2, ">temp2.sampled.txt") or die("Could not write to the second samnpled file!\n");

	my $count = 0;

	
	my $line1;	
	my $line2;

	while($line1 = <F1>)
	{

		if($count % $freq == 0)
		{
			print S1 "$line1";
		}


		$count++;
	}
	$count = 0;

	while($line2 = <F2>)
        {

                if($count % $freq == 0)
                {
                        print S2 "$line2";
                }

		$count++;
        }

	close F1;
	close F2;
	close S1;
	close S2;
}

	

sub array_abs
{

	my $data1 = $_[0];
	my @ret = ();

	for (my $i = 0; $i < scalar(@$data1); $i++)
	{
		push @ret, abs($$data1[$i]);
	}

	return \@ret;

}
		
	
sub array_diff
{
	
	my $data1 = $_[0];
	my $data2 = $_[1];
	
	my @ret = ();
	#These arrays should be the same length

	for(my $i = 0; $i < scalar(@$data1); $i++)
	{
		push @ret, abs($$data1[$i] - $$data2[$i]);
	}
	
	return \@ret;

}


sub calc_cdf
{
	my $data_ref = $_[0];
	my $data_all_ref = $_[1];
	
	my @cdf = ();

	my $ge_count = 0;
	#For each element in data_all_ref, see how many elements in data_ref it is >= to
	my $iter_count = scalar(@$data_all_ref);

	#print "Outer iterations remaining...$iter_count\n";
	for(my $i = 0; $i < scalar(@$data_all_ref); $i++)
	{

		if($i % 1000 == 0)
		{
			print "Outer iterations remaining...".($iter_count - $i)."\n";
		}
		
		$ge_count = 0;
		for(my $j = 0; $j < scalar(@$data_ref); $j++)
		{
			if($$data_all_ref[$i] >= $$data_ref[$j])
			{
				$ge_count++;
			}



		}

		push @cdf, ($ge_count/scalar(@$data_ref));


	}	

	return \@cdf;


}

sub qks
{

	my $lambda = $_[0];
	my $iter = $_[1]; #Number of iterations to approximate Q function

	my $val = 0;
	for(my $j = 1; $j <= $iter; $j++)
	{

		$val = (-1**($j-1)) * exp(-2 * $j*$j * $lambda * $lambda);
		
		$val += $val;

	}

	$val *= 2;

	return $val;

}

sub cleanup
{
        system("rm temp1.txt");
        system("rm temp2.txt");
        system("rm temp1.sampled.txt");
        system("rm temp2.sampled.txt");
        system("rm temp1.zeros.txt");

}

sub average
{

	my $file = $_[0];
	my $ID = $_[1];

	open(F1, "<$file") or die("average: Could not open the input file!\n");
	#open(F2, ">$ID.txt") or die("average: Could not write to the average file!\n");
	
	my $count = 0;
	my $n = 0;
	my $sum;
	my $line;
	my $val;

	my @record = ();

	while($line = <F1>)
	{
		chomp($line);
		@record = split(/\t/, $line);
		$val = $record[2];

		if($count % 100 == 0)
		{	
			$sum += $val;
			$n++;	
			#print F2 "$sum\n";
		}
		$count++;
	}
	print "Sample size: $n\n";

	#close F2;
	close F1;
	return $sum/$n;

}

sub average1
{
	#This function calculates the average of an amplicon ONLY in betwen sub-amplicons originally identified by rdxplorer

	my $chr = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	my $file1 = $_[3];	
	my $bam_file = $_[4];
	my $ID = $_[5];


	open(F1, "<$file1") or die("average1: Could not open file1!\n");
	#open(F2, "<$file2") or die("average1: Could not open file2!\n");
	my $count = 0;

	#system("rm temp1.final.txt");
	my $sum;
	my $line;
	my $line2;
	my @record = ();
	
	my $chr1;	
	my $start1;
	my $end1;

	while($line = <F1>)
	{
		chomp($line);
		@record = split(/	/, $line);
		$chr1 = $record[7];
		$start1 = $record[8];
		$end1 = $record[9];
		
		#print "VALUES : $chr1	$start1	$end1\n";
		#print "VALUES1 : $chr   $start $end\n";
		if($chr eq "chr".$chr1 && $start1 >= $start && $end1 <= $end)
		{

			print "IN IF\n";
			#DUMP BAM COVERAGE
			system("~/samtools-0.1.18/samtools depth -r chr$chr1:$start1-$end1 -Q 30 $bam_file > temp1.txt");
			system("cat temp1.txt | awk 'BEGIN { prev_chr=\"\";prev_pos=0;} { if(\$1==prev_chr && prev_pos+1!=int(\$2)) {for(i=prev_pos+1;i<int(\$2);++i) {printf(\"%s\\t%d\\t0\\n\",\$1,i);}} print; prev_chr=\$1;prev_pos=int(\$2);}' > temp1.zeros.txt");		
			system("cat temp1.zeros.txt >> temp1.final.txt");
		}
	}
	
	close F1;

	#system("cat temp1.txt | awk 'BEGIN { prev_chr=\"\";prev_pos=0;} { if(\$1==prev_chr && prev_pos+1!=int(\$2)) {for(i=prev_pos+1;i<int(\$2);++i) {printf(\"%s\\t%d\\t0\\n\",\$1,i);}} print; prev_chr=\$1;prev_pos=int(\$2);}' > temp1.zeros.txt");	

	system("cp temp1.final.txt AMPLICON$file_count.txt");
	$file_count++;

	return average("temp1.final.txt", $ID);
	#system("cp temp1.final.txt AMPLICON$file_count.txt");
	#$file_count++;
	#print "AMPLICON:$file_count\n"

}

sub ks
{
	
	#This function should read the read depth files and perform the ks test

	my $file1 = $_[0];
	my $file2 = $_[1];

	my @data1 = ();
	my @data2 = ();

	my @data_all = ();
	my $length1;
	my $length2;
	my $line;

	my $cdf1 = ();	#reference to arrays
	my $cdf2 = ();
#
	open(D1, "<$file1") or die("Could not open the first data file!\n");
	open(D2, "<$file2") or die("Could not open the second data file!\n");

	print "Doing KS test\n";
	while($line = <D1>)
	{

		chomp($line);
		my @tmp = split(/	/, $line);
		push @data1, $tmp[2];

	}

	while($line = <D2>)
        {

                chomp($line);
                my @tmp = split(/	/, $line);
                push @data2, $tmp[2];

        }
		
	close D1;
	close D2;
	
	print "Data retrieved...\n";
	my $num1 = scalar(@data1);
	my $num2 = scalar(@data2);
	
	@data1 = sort {$a <=> $b} @data1;
	@data2 = sort {$a <=> $b} @data2;
	
	print "Sorting complete...\n";

	@data_all = @data1;

	push @data_all, @data2;

	print "Calculating CDFs...\n";
	$cdf1 = calc_cdf(\@data1, \@data_all);
	$cdf2 = calc_cdf(\@data2, \@data_all);
	print "done.\n";

	my $diff = array_diff(\@data1, \@data2);
	my $ab = array_abs($diff);
	
	my $d = max(@$ab);

	my $c = sqrt(($num1*$num2)/($num1 + $num2));

	print "C: $c\n";
	print "D: $d\n";
	#die;
		
	my $lambda = ($c + 0.12 + 0.11/$c) * $d;
	
	print "P-VALUE: ".qks($lambda, 1000)."\n";	
	return qks($lambda, 1000);


}

sub dfs
{

	
	my @adj = $g->neighbors($_[0]);
	push @visited_vertices, $_[0];

	$explored{$_[0]} = 1;

	for(my $i = 0; $i < scalar(@adj); $i++)
	{
		if($explored{$adj[$i]} == 0)
		{

			$h->add_weighted_edge($_[0], $adj[$i], 1);
			$dviz->add_edge_once($_[0], $adj[$i]);

			#print "CASE 1: ADDED DIRECTED EDGE ".$_[0]." ---> ".$adj[$i]."\n";
			push @edges, [$_[0], $adj[$i], 1];

			#print "1st PUSHED ".$_[0]."	".$adj[$i]."\n";
			$count_edges_added++;			

			$explored_edges{$_[0]." ".$adj[$i]} = 1;
			$explored_edges{$adj[$i]." ".$_[0]} = 1;
			#print $explored_edges{$_[0]." ".$adj[$i]}."\n";
			#die;

			dfs($adj[$i]);

		}

		else
		{
			if(( (not defined $explored_edges{$adj[$i]." ".$_[0]}) && (not defined $explored_edges{$_[0]." ".$adj[$i]})) || ($explored_edges{$adj[$i]." ".$_[0]} == 0 && $explored_edges{$_[0]." ".$adj[$i]} == 0))
			{
				#Edge does not exist. Add it.
				 $h->add_weighted_edge($_[0], $adj[$i], 1);			
				$dviz->add_edge_once($_[0], $adj[$i]);

				 $explored_edges{$adj[$i]." ".$_[0]} = 1;
				 push @edges, [$_[0], $adj[$i], 1];
				
				$count_edges_added++;
	
			}

			
			
		}

	}		

}
if(scalar(@ARGV) != 7)
{
        die("Usage is perl program.pl [SV FILE] [CN SEGMENT FILE] [WINDOW SIZE] [BAM FILE] [MINQUAL] [MIN CYCLIC] [MIN NON CYCLIC]\n");
}

open(SV, "<".$ARGV[0]) or die ("Could not open SV file!\n");
open(CN, "<".$ARGV[1]) or die ("Could not open CN file!\n");

my $line = "";
my $temp_line = "";

my $window = $ARGV[2];

my $bam_file = $ARGV[3];
my $min_qual = $ARGV[4];
my $min_cyclic = $ARGV[5];
my $min_non_cyclic = $ARGV[6];
my $other_chr = "";
my $other_pos = 0;
my $cn_chr = "";
my $cn_start = "";
my $cn_end = "";

my $sv_chr_i = "";
my $sv_chr_i_bp = "";

my $e = "";
my $sv_chr_j = "";
my $sv_chr_j_bp = "";

my @record = ();
my @seg_record = (); #will store all amplified segment records
my $startFlag = 0; #1 if graph constructor should look at start of a CN segment, 0 otherwise.
my @SV = ();

my @temp_rec = ();
my @temp_rec_other = ();

#Note the following property: all segments participating in a DM should have SV links at both ends that are proximal to other segments
#Thus, to preprocess, we can exclude CN segments that don't meet this criteria
#The graph construction algorithm works as follows: 
#	#For the initial segment, "side" is arbitrary (could be either start or end)
# 	For all CN segments i = 1 to N
# A:		Look for SV records (on both sides) that are proximal to side(i)
#			If none
#				then fetch next record (i.e. do nothing and wait for next iteration)
#			ELSE
#				# we found a potential link
#				Search for CN segments j that are proximal to the OTHER side of the SV prediction
#					If none found	
#						then get next segment (i.e. do nothing and wait for next iteration)
#					ELSE
#						if proximal CN side is at end of segment, then side = start, else side = end
#						add directed edge (i, j)
#						i = j
#						goto A
						
#				

@SV = <SV>;

my $seg_line = "";
my @seg_rec = ();
while($line = <CN>)
{
	@record = split(/\t/, $line);	
	push @seg_record, $line;
}


#@seg_record now contains ALL amplified segments


for (my $i = 0; $i < scalar(@seg_record)-1; $i++)
{
			$seg_line = $seg_record[$i];
			chomp($seg_line);
			@seg_rec = split(/\t/, $seg_line);

				
			for(my $k = 0; $k < scalar(@SV); $k++)
			{


				$line = $SV[$k];
				#print "$line\n";
				chomp($line);
		
				#print "link: $line\n";
				@temp_rec = split(/\t/, $line);
			
				#print $temp_rec[7]."\n";

				#print "FIRST: temp_rec[6]: ".$temp_rec[6]."\nseg_rec[7] ".$seg_rec[7]."\ntemp_rec[7]: ".$temp_rec[7]."\nseg_rec[8] ".$seg_rec[8]."\n";

				if(($temp_rec[0] eq $seg_rec[0] || $temp_rec[0] eq "chr".$seg_rec[0]) && abs($temp_rec[1] - $seg_rec[1]) <= $window)
				{
		
					#Start of CN segment has link 
					#See if the link goes to another segment
			

					#print "seg_rec[7]: ".$seg_rec[0]." seg_rec[8]: ".$seg_rec[1]."\n";	
					$other_chr = $temp_rec[3];
					$other_pos = $temp_rec[4];
				
					for(my $j = $i+1; $j < scalar(@seg_record); $j++)
				        {
						#Check other CN segments for link

						$temp_line = $seg_record[$j];

						chomp($temp_line);							
						@temp_rec_other = split(/\t/, $temp_line);

						#$print "temp_rec_other[7]: ".$temp_rec_other[7]."\n";
						#print "seg_line: $seg_line\n";
						#print "FIRST: other_chr: ".$other_chr." temp_rec_other[7] chr".$temp_rec_other[0]." other_pos: ".$other_pos." temp_rec_other[8]: ".$temp_rec_other[1]."\n";
						if(($other_chr eq "chr".$temp_rec_other[0] || $other_chr eq $temp_rec_other[0]) && abs($other_pos - $temp_rec_other[1]) <= $window)
        		                        {
									
							$g->add_edge($seg_record[$i], $seg_record[$j]);								
							$explored{$seg_record[$i]} = 0;
							$explored{$seg_record[$j]} = 0;
							$explored_edges{$seg_record[$i]." ".$seg_record[$j]} = 0;

							$viz->add_edge_once($seg_record[$i], $seg_record[$j]);
							goto A;


						}
						 #print "seg_line: $seg_line\n";
						# print "SECOND: other_chr: ".$other_chr." temp_rec_other[7] ".$temp_rec_other[0]." other_pos: ".$other_pos." temp_rec_other[9]: ".$temp_rec_other[2]."\n";
						if(($other_chr eq "chr".$temp_rec_other[0] || $other_chr eq $temp_rec_other[0]) && abs($other_pos - $temp_rec_other[2]) <= $window) #check other end of CN segment
						{

							#Adding edges is an idempotent operation, so I removed the "elsif"
							#SV link at start of current CN segment goes to end of other CN segment
							# print "Edge added 2: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
							$g->add_edge($seg_record[$i], $seg_record[$j]);

							 $explored{$seg_record[$i]} = 0;
                                                        $explored{$seg_record[$j]} = 0;
							$explored_edges{$seg_record[$i]." ".$seg_record[$j]} = 0;
							 $viz->add_edge_once($seg_record[$i], $seg_record[$j]);
							goto A;

						}

					}

				}
				elsif(($temp_rec[3] eq "chr".$seg_rec[0] || $temp_rec[3] eq $seg_rec[0]) && abs($temp_rec[4] - $seg_rec[1]) <= $window)
                                {

					#print "Edge added 4: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
                                        #End of CN segment has link
                                        #See if the link goes to another segment

					#print "seg_rec[7]: ".$seg_rec[0]." seg_rec[8]: ".$seg_rec[8]."\n";
                                        $other_chr = $temp_rec[0];
                                        $other_pos = $temp_rec[1];

                                        for(my $j = $i+1; $j < scalar(@seg_record); $j++)
                                        {
                                                #Check other CN segments for link

                                                $temp_line = $seg_record[$j];

                                                chomp($temp_line);
                                                @temp_rec_other = split(/\t/, $temp_line);

						
						 #print "seg_line: $seg_line\n";
					#	print "THIRD: other_chr: ".$other_chr." temp_rec_other[7] ".$temp_rec_other[7]." other_pos: ".$other_pos." temp_rec_other[8]: ".$temp_rec_other[8]."\n";
						 if(($other_chr eq "chr".$temp_rec_other[0] || $other_chr eq $temp_rec_other[0]) && abs($other_pos - $temp_rec_other[1]) <= $window)
                                                {

                                         	
				                       #ADD EDGE
                                                        #SV link at end of current CN segment goes to start of other CN segment

							 #print "Edge added 3: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
                                                        $g->add_edge($seg_record[$i], $seg_record[$j]);

							 $explored{$seg_record[$i]} = 0;
                                                        $explored{$seg_record[$j]} = 0;
							$explored_edges{$seg_record[$i]." ".$seg_record[$j]} = 0;
							 $viz->add_edge_once($seg_record[$i], $seg_record[$j]);
							goto A;
                                                }
						 #print "seg_line: $seg_line\n";
						#print "FOURTH: other_chr: ".$other_chr." temp_rec_other[7] ".$temp_rec_other[7]." other_pos: ".$other_pos." temp_rec_other[9]: ".$temp_rec_other[9]."\n";
						if(($other_chr eq "chr".$temp_rec_other[0] || $other_chr eq $temp_rec_other[0]) && abs($other_pos - $temp_rec_other[2]) <= $window) #check other end of CN segment
						{	
							# print "Edge added 4: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
                                                        $g->add_edge($seg_record[$i], $seg_record[$j]);
							$explored{$seg_record[$i]} = 0;
                                                        $explored{$seg_record[$j]} = 0;
							$explored_edges{$seg_record[$i]." ".$seg_record[$j]} = 0;
							 $viz->add_edge_once($seg_record[$i], $seg_record[$j]);
							goto A;
                                                }

                                        }

                                }

				 elsif(($temp_rec[0] eq "chr".$seg_rec[0] || $temp_rec[0] eq $seg_rec[0]) && abs($temp_rec[1] - $seg_rec[2]) <= $window)
                                {

                                        #Start of CN segment has link
                                        #See if the link goes to another segment


                                       # print "seg_rec[7]: ".$seg_rec[7]." seg_rec[9]: ".$seg_rec[9]."\n";
                                        $other_chr = $temp_rec[3];
                                        $other_pos = $temp_rec[4];

                                        for(my $j = $i+1; $j < scalar(@seg_record); $j++)
                                        {
                                                #Check other CN segments for link

                                                $temp_line = $seg_record[$j];

                                                chomp($temp_line);
                                                @temp_rec_other = split(/\t/, $temp_line);

                                                #$print "temp_rec_other[7]: ".$temp_rec_other[7]."\n";
                                                #print "seg_line: $seg_line\n";
                                               # print "FIFTH: other_chr: ".$other_chr." temp_rec_other[7] ".$temp_rec_other[7]." other_pos: ".$other_pos." temp_rec_other[8]: ".$temp_rec_other[8]."\n";
                                        	if(($other_chr eq "chr".$temp_rec_other[0] || $other_chr eq $temp_rec_other[0]) && abs($other_pos - $temp_rec_other[1]) <= $window)
                                                {
                                                

                                                                #ADD EDGE
                                                        #SV link at start of current CN segment goes to start of other CN segment

                                                #        print "Edge added 5: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
                                                        $g->add_edge($seg_record[$i], $seg_record[$j]);
                                                        $explored{$seg_record[$i]} = 0;
                                                        $explored{$seg_record[$j]} = 0;
							$explored_edges{$seg_record[$i]." ".$seg_record[$j]} = 0;
							 $viz->add_edge_once($seg_record[$i], $seg_record[$j]);
                                                        goto A;


                                                }
                                                 #print "seg_line: $seg_line\n";
                                           #      print "SIXTH: other_chr: ".$other_chr." temp_rec_other[7] ".$temp_rec_other[7]." other_pos: ".$other_pos." temp_rec_other[9]: ".$temp_rec_other[9]."\n";
                                
						 if(($other_chr eq "chr".$temp_rec_other[0] || $other_chr eq $temp_rec_other[0]) && abs($other_pos - $temp_rec_other[2]) <= $window) #check other end of CN segment
                                                {
                                                        #Adding edges is an idempotent operation, so I removed the "elsif"
                                                        #SV link at start of current CN segment goes to end of other CN segment
                                                        # print "Edge added 6: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
                                                        $g->add_edge($seg_record[$i], $seg_record[$j]);

                                                         $explored{$seg_record[$i]} = 0;
                                                        $explored{$seg_record[$j]} = 0;
							$explored_edges{$seg_record[$i]." ".$seg_record[$j]} = 0;
							 $viz->add_edge_once($seg_record[$i], $seg_record[$j]);
                                                        goto A;

                                                }

                                        }

                                }
				elsif(($temp_rec[3] eq "chr".$seg_rec[0] || $temp_rec[3] eq $seg_rec[0]) && abs($temp_rec[4] - $seg_rec[2]) <= $window)
                                {

                                        #print "Edge added 4: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";


                                        #End of CN segment has link
                                        #See if the link goes to another segment

                                        #print "seg_rec[7]: ".$seg_rec[7]." seg_rec[9]: ".$seg_rec[9]."\n";
                                        $other_chr = $temp_rec[0];
                                        $other_pos = $temp_rec[1];

                                        for(my $j = $i+1; $j < scalar(@seg_record); $j++)
                                        {
                                                #Check other CN segments for link

                                                $temp_line = $seg_record[$j];

                                                chomp($temp_line);
                                                @temp_rec_other = split(/\t/, $temp_line);


                                                 #print "seg_line: $seg_line\n";
                                                #print "SEVEN: other_chr: ".$other_chr." temp_rec_other[7] ".$temp_rec_other[7]." other_pos: ".$other_pos." temp_rec_other[8]: ".$temp_rec_other[8]."\n";
                                		 if(($other_chr eq "chr".$temp_rec_other[0] || $other_chr eq $temp_rec_other[0]) && abs($other_pos - $temp_rec_other[1]) <= $window)
                                                {

                                                        #SV link at end of current CN segment goes to start of other CN segment

                                                        # print "Edge added 7: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
                                                        $g->add_edge($seg_record[$i], $seg_record[$j]);

                                                         $explored{$seg_record[$i]} = 0;
                                                        $explored{$seg_record[$j]} = 0;
							$explored_edges{$seg_record[$i]." ".$seg_record[$j]} = 0;
							 $viz->add_edge_once($seg_record[$i], $seg_record[$j]);
                                                        goto A;
                                                }
                                                 #print "seg_line: $seg_line\n";
                                            #    print "EIGHTH: other_chr: ".$other_chr." temp_rec_other[7] ".$temp_rec_other[7]." other_pos: ".$other_pos." temp_rec_other[9]: ".$temp_rec_other[9]."\n";
                                                if($other_chr eq "chr".$temp_rec_other[0] && abs($other_pos - $temp_rec_other[2]) <= $window) #check other end of CN segment
                                                {
                                                        #Adding edges is an idempotent operation, so I removed the "elsif"
                                                        #SV link at end of current CN segment goes to end of other CN segment
                                                 #        print "Edge added 8: \n".$seg_record[$i]."\n".$seg_record[$j]."\n";
                                                        $g->add_edge($seg_record[$i], $seg_record[$j]);
                                                        $explored{$seg_record[$i]} = 0;
                                                        $explored{$seg_record[$j]} = 0;
							$explored_edges{$seg_record[$i]." ".$seg_record[$j]} = 0;
							 $viz->add_edge_once($seg_record[$i], $seg_record[$j]);
                                                        goto A;
                                                }

                                        }

                                }
				

			A:
			}
	#}


}


#Undirected graph built
#DEBUG ABOVE CODE: 2/19/2013

#Now seg_record contains all records representing a copy number gain


my @V = $g->vertices;
#print "seg_record: ".$V[0]."\n";

### HERE: may do DFS for all vertices 
#9/1/14

foreach my $e (@V)
{
	dfs($e);
}

my $u = 0;
my $v = 0;

for($e = 0; $e < scalar(@edges); $e++)
{
	$u = $edges[$e][0];
	$v = $edges[$e][1];

}

#print "SIZE OF EDGES: ".scalar(@edges)."\n";
for($e = 0; $e < scalar(@edges); $e++)
{
        $u = $edges[$e][0];
        $v = $edges[$e][1];

}

close CN;

#print "The undirected graph is $g\n";
#print "The directed graph is $h\n";
#print "Num vertices: ".$g->vertices."\n";

#print "input vertex: ".$seg_record[0]."\n";
my @adj = ();
@adj = $g->neighbors($V[0]);


#my $apsp = $g->SP_Dijkstra();

my $apsp = $h->APSP_Floyd_Warshall(); #All pairs shortest path now calculated

#Now calculate the SCCs
my @scc = $h->strongly_connected_components; #scc contains list of anonymous arrays
#print "number of SCCs is ".scalar(@scc)."\n";

#my @shortest_path = $h->SP_Dijkstra($adj[0], $seg_record[0]);

#print "-------PATH-------\n";
#print @shortest_path;
#print "\n";

#print "adj: ".$adj[0]."\n";
#print "seg_record: ".$V[0]."\n";
#print "2-1 weight: ".$h->get_edge_weight($adj[0], $V[0])."\n";
#print "1-2 weight: ".$h->get_edge_weight($V[0], $adj[0])."\n";

#$apsp->path_length($u, $v);
#print "EDGES PUSHED ON STACK: $count_edges_added\n";

#print "adj: ".$adj[0]."\n";
#print "List of vertices: ".@($g->vertices)."\n";



#print @scc."\n";
my @shortest_path;
my $path_total_weight = INFINITY;
my $temp_weight = INFINITY;


foreach my $e (@scc) #Cycle through all SCCs in the graph to find potential DMs
{

	
	#print "SCC: @$e\n";
	#print "size: ".scalar(@$e)."\n";
	
	@shortest_path = ();
	$temp_weight = INFINITY;
	$path_total_weight = INFINITY;
	foreach my $n (@$e)
	{
		#print "IN FIRST FOR\n";
		@adj = $h->neighbors($n);
		#We now have all the neighbors of the vertex n
	
		foreach my $m (@adj)
		{
			#for each neighbor of n

			#print "IN THIRD FOR\n";
			#if($h->get_edge_weight($n, $m) != INFINITY)
			if((defined $h->get_edge_weight($n, $m)) && $h->get_edge_weight($n, $m) ne '' && (defined $apsp->path_length($m, $n)) && $apsp->path_length($m, $n) ne '')
			{

			#Now return shortest path from m to n
				
				$temp_weight = $h->get_edge_weight($n, $m) + $apsp->path_length($m, $n);			
				#print "THE EDGE WEIGHT OF INTEREST IS :".$h->get_edge_weight($n, $m)."\n";
				#print "temp_weight: $temp_weight\npath_total_weight: $path_total_weight\n";
				if($temp_weight < $path_total_weight)
				{
					$path_total_weight = $temp_weight;			
					@shortest_path = $apsp->path_vertices($m, $n);			
					#print "PATH VERTICES ARE @shortest_path\n";		
				}

			}	



		}



	}	


		for(my $i = 0; $i < scalar(@shortest_path); $i++)
                {

                        chomp($shortest_path[$i]);

                        my @rec_i = split(/\t/, $shortest_path[$i]);
                        my $chr = "chr".$rec_i[0];
                        my $start = $rec_i[1];
                        my $end = $rec_i[2];

                        system("samtools depth -r $chr:$start-$end -Q $min_qual $bam_file > temp1.txt");
                        system("cat temp1.txt | awk 'BEGIN { prev_chr=\"\";prev_pos=0;} { if(\$1==prev_chr && prev_pos+1!=int(\$2)) {for(i=prev_pos+1;i<int(\$2);++i) {printf(\"%s\\t%d\\t0\\n\",\$1,i);}} print; prev_chr=\$1;prev_pos=int(\$2);}' > temp1.zeros.txt");


                        print average("temp1.zeros.txt", $i)."\n";
                        if(average("temp1.zeros.txt", $i) >= $max_cov)
                        {
                                $max_cov = average("temp1.zeros.txt", $i);
                        }
                        if(average("temp1.zeros.txt", $i) <= $min_cov)
                        {
                                $min_cov = average("temp1.zeros.txt", $i);
                        }



                }

	if(scalar(@shortest_path) > $min_cyclic && ($max_cov - $min_cov) < 35) #Once again, this is to ensure we don't get a predicted dmin that has only one amplicon 
        {	
		print "SEGMENTS FOR PREDICTED DOUBLE MINUTE: @shortest_path\n*************************************\n";
		#print REPORT "SEGMENTS FOR PREDICTED DOUBLE MINUTE: @shortest_path\n*************************************\n";
		$dm_count++;
	}

	$min_cov = 100000000;
	$max_cov = 0;
	my $ssize = scalar(@shortest_path);
	$h = $h->delete_path(@shortest_path);

	@shortest_path = ();
	open(VIZ, ">>graph.dot") or die ("Could not create the graph html file!\n");

	if($viz->is_directed())
	{
		#die("GRAPH IS DIRECTED");
	}
	print VIZ $viz->as_graphviz();
	
	close VIZ;
	#goto C;
	B:
}

my @wcc = $h->weakly_connected_components();


foreach my $e (@wcc)
{
  
@shortest_path = @$e;

for(my $i = 0; $i < scalar(@shortest_path); $i++)
{
                        #chomp($shortest_path[$j]);
                        chomp($shortest_path[$i]);

                        my @rec_i = split(/\t/, $shortest_path[$i]);
                        my $chr = "chr".$rec_i[0];
                        my $start = $rec_i[1];
                        my $end = $rec_i[2];

                        system("samtools depth -r $chr:$start-$end -Q 30 $bam_file > temp1.txt");
                        system("cat temp1.txt | awk 'BEGIN { prev_chr=\"\";prev_pos=0;} { if(\$1==prev_chr && prev_pos+1!=int(\$2)) {for(i=prev_pos+1;i<int(\$2);++i) {printf(\"%s\\t%d\\t0\\n\",\$1,i);}} print; prev_chr=\$1;prev_pos=int(\$2);}' > temp1.zeros.txt");


                     #   print average("temp1.zeros.txt", $i)."\n";

			 if(average("temp1.zeros.txt", $i) >= $max_cov)
                        {
                                $max_cov = average("temp1.zeros.txt", $i);
                        }
                        if(average("temp1.zeros.txt", $i) <= $min_cov)
                        {
                                $min_cov = average("temp1.zeros.txt", $i);
                        }

}


	if(scalar(@shortest_path) > $min_non_cyclic && ($max_cov - $min_cov) < 35)
        {
                print "SEGMENTS FOR PREDICTED INCOMPLETE DOUBLE MINUTE: @shortest_path\n*************************************\n";
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
print "There are $dm_count predicted double minutes\n";
#close REPORT;
close SV;

