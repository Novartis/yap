#!/usr/bin/env perl 
=license information
Copyright 2014 Novartis Institutes for Biomedical Research

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=cut

if(@ARGV <4){
    die "Arguments for microrna expression counts  less than 4\n\n";
}

@input_arr = <STDIN>;
$miR_precursor_pos_annot_file = $ARGV[1];
$output_file = $ARGV[3];
print $#input_arr;
    @l = open_file ($miR_precursor_pos_annot_file);
    %hairpin_start = ();
    %hairpin_end = ();
    %hairpin_mature_or_star = ();

    foreach $e (@l){
	@d = split /\t/, $e;
	if(@d != 5){ die "died error not 6  ele in $miR_precursor_pos_annot_file !!!!";}
	$hp = $d[0];
	$matu = $d[1];
	$st = $d[3];
	$ed = $d[4];

	$hairpin_start{$hp} = $st;
	$hairpin_end{$hp} = $ed;
	$hairpin_mature_or_star{$hp} = $matu;
    }
    ##### 3.2 filter the bowtie output to determine whether a hit occ in mature miR or miR* area

    $shift_cutoff = 4;

    %miR_star_2_count_individual = ();
    %miR_mature_2_count_individual = ();




    $line = ();
    $c = 0;
    $pre_read_id = (); # for block detection
    $cur_read_id = (); # for block detection
    foreach $line (@input_arr){ #note the bowtie output is sorted based on the read ids
	if ($line =~ /^@.*/){next;} 
	chomp($line);
	@d = split /\t/, $line;
	#if(@d < 2 ){die "not 9 columns in sorted aligner output file  at line: $c \n";}
	$strand = $d[1];
	if($strand ne '0'){	next;    } # consider + strand only
	$cur_read_id = $d[0];
	$strand = $d[1];
	$start_pos = $d[3] + 1;
	$c ++;     

	if ($c ==1){
	    $pre_read_id = $d[0]; #set pre pointer
	    # check the current line
	    $cur_read_len = length($d[9]);
	    $end_pos = $d[3]+ $cur_read_len;
	    $cur_hp = $d[2];

	    ## ******add block success count	
		    $block_success_miR_or_miR_star_hs = ();

            if(!exists $hairpin_start{$cur_hp}) {die "error, $cur_hp not found in miRbase !!!!\n\n";}

            @mature_or_star_sts = split/,/, $hairpin_start{$cur_hp};
            @mature_or_star_eds = split/,/, $hairpin_end{$cur_hp};
            @mature_or_star_ids = split/,/, $hairpin_mature_or_star{$cur_hp};

            for ($j = 0; $j <= @mature_or_star_sts -1; $j++){
                $shift_st = abs($start_pos - $mature_or_star_sts[$j]);
                $shift_ed = abs($end_pos - $mature_or_star_eds[$j]);
                if ($shift_st <= $shift_cutoff  && $shift_ed  <= $shift_cutoff){
                    $block_success_miR_or_miR_star_hs{$mature_or_star_ids[$j]} = 1;


                }
            }
        }


		   
	else{ #other lines than line 1
	    if($cur_read_id eq $pre_read_id){ # continue in this block

		    $cur_read_len = length($d[9]);
		    $end_pos = $d[3]+ $cur_read_len;
		    $cur_hp = $d[2];
		    #$block_success_miR_or_miR_star_hs{$cur_hp} = 1;
		    #$block_success_miR_or_miR_star_hs{$cur_hp} = 1;
		  if(!exists $hairpin_start{$cur_hp}) {die "error, $cur_hp not found in miRbase !!!!\n\n";}
			  @mature_or_star_sts = split/,/, $hairpin_start{$cur_hp};
                    @mature_or_star_eds = split/,/, $hairpin_end{$cur_hp};
                    @mature_or_star_ids = split/,/, $hairpin_mature_or_star{$cur_hp};

                    for ($j = 0; $j <= @mature_or_star_sts -1; $j++){
                        $shift_st = abs($start_pos - $mature_or_star_sts[$j]);
                        $shift_ed = abs($end_pos - $mature_or_star_eds[$j]);
                        if ($shift_st <= $shift_cutoff  && $shift_ed  <= $shift_cutoff){

                            $block_success_miR_or_miR_star_hs{$mature_or_star_ids[$j]} = 1;


                        }
                    }

            }

 
	    else{ #start new block
		#first process prev block  ******
		$number_of_success_miR_or_miR_star_in_block = scalar(keys(%block_success_miR_or_miR_star_hs));
		if($number_of_success_miR_or_miR_star_in_block > 0 ){
		    $read_redundancy_count = 1;
	       	#	$read_redundancy_count = $d[5]; #modify here
		    foreach $k (keys(%block_success_miR_or_miR_star_hs)){
			if(index($k, '*') >=0){
			    if(!exists($miR_star_2_count_individual{$k})){
			
				$miR_star_2_count_individual{$k} = $read_redundancy_count/$number_of_success_miR_or_miR_star_in_block;
			    }
			    else{
				$miR_star_2_count_individual{$k} += $read_redundancy_count/$number_of_success_miR_or_miR_star_in_block; 
			    }
		    
			}
			else{
			    if (!exists($miR_mature_2_count_individual{$k})){
				$miR_mature_2_count_individual{$k} = $read_redundancy_count/$number_of_success_miR_or_miR_star_in_block;
			    }
			    else{
				$miR_mature_2_count_individual{$k} += $read_redundancy_count/$number_of_success_miR_or_miR_star_in_block;
				#print "detected: $k\t$miR_mature_2_count_individual{$k}\n";
			    }
		       
			}
		    }
	        }

		## *** reset block
		%block_success_miR_or_miR_star_hs = ();
		
		$pre_read_id = $d[0]; #set pre pointer
		# check the current line
		$cur_read_len = length($d[9]);
		$end_pos = $d[3]+ $cur_read_len;
		$cur_hp = $d[2];
   if(!exists $hairpin_start{$cur_hp}) {die "error, $cur_hp not found in miRbase !!!!\n\n";}

                @mature_or_star_sts = split/,/, $hairpin_start{$cur_hp};
                @mature_or_star_eds = split/,/, $hairpin_end{$cur_hp};
                @mature_or_star_ids = split/,/, $hairpin_mature_or_star{$cur_hp};


                for ($j = 0; $j <= @mature_or_star_sts -1; $j++){
                        $shift_st = abs($start_pos - $mature_or_star_sts[$j]);
                        $shift_ed = abs($end_pos - $mature_or_star_eds[$j]);
                        if ($shift_st <= $shift_cutoff  && $shift_ed  <= $shift_cutoff){

                            $block_success_miR_or_miR_star_hs{$mature_or_star_ids[$j]} = 1;

                        }
                }


            }
        }
	    
		
}


    ### *** process the last block
    $number_of_success_miR_or_miR_star_in_block = scalar(keys(%block_success_miR_or_miR_star_hs));
       $read_redundancy_count = 1;
	if($number_of_success_miR_or_miR_star_in_block > 0 ){
	foreach $k (keys(%block_success_miR_or_miR_star_hs)){
	    if(index($k, '*') >=0){
		       if(!exists($miR_star_2_count_individual{$k})){
			   $miR_star_2_count_individual{$k} = $read_redundancy_count/$number_of_success_miR_or_miR_star_in_block;
		       }
		       else{
			   $miR_star_2_count_individual{$k} += $read_redundancy_count/$number_of_success_miR_or_miR_star_in_block; 
		       }
		    
		   }
	    else{
		       if (!exists($miR_mature_2_count_individual{$k})){
			   $miR_mature_2_count_individual{$k} = $read_redundancy_count/$number_of_success_miR_or_miR_star_in_block;
		       }
		       else{
			   #print "$number_of_success_miR_or_miR_star_in_block\n";
			   $miR_mature_2_count_individual{$k} += $read_redundancy_count/$number_of_success_miR_or_miR_star_in_block;
		       }
		       
		   }
	}
    }
    
    open (FILE, ">$output_file") || die "error, cannot open file for reading: $output_file !";
   #### step 4 print out results ########
    $num_of_miR_mature_species = scalar(keys(%miR_mature_2_count_individual));
    $num_of_miR_star_species = scalar(keys(%miR_star_2_count_individual));

    if ($num_of_miR_mature_species > 0 ){
	foreach $k (keys(%miR_mature_2_count_individual)){
	    print FILE "$k\t$miR_mature_2_count_individual{$k}\n";
	}
    }

    if ($num_of_miR_star_species > 0 ){
	foreach $k (keys(%miR_star_2_count_individual)){
	    print FILE "$k\t$miR_star_2_count_individual{$k}\n";
	}
    }
close FILE;
   


   ##### step 5: rm temp files in the input dir(current working directory

   #`rm temp.$in_f.fa   temp.$in_f.fa.bowtie_out`;







sub open_file                                                              
{                                                                          
    my $file = $_[0];                                                  
    open (FILE, "$file") || die "error, cannot open file for reading: $file !";    
    @file = <FILE>;                                                    
    chomp @file;                                                       
    close (FILE) || die "error, cannot close file: $file $!";                       
    return @file;                                                      
}  






