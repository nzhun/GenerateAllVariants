#!/user/bin/perl

use strict;
use warnings;

our %rate=();
our %gene_rate=();
our $cutoff=0.025;
our @nucleotides=("A","T","C","G");

our %trichanges=();
our %amino_tables=();
my $path="/home/local/ARCS/nz2274";
my $file_rate="$path/Resources/mutation_table/3mer_table.txt";
## read amino acid change table
my $file_aminoacid="$path/Resources/mutation_table/amino_acid_codon.txt";


my %complement=(
"A" =>"T",
"T" =>"A",
"C" => "G",
"G" => "C"	
);


sub main {
		#my $dbnsfp="/home/local/ARCS/nz2274/PSAP/Data/exon_transcript.dbnsfp30a.bed.gz";
		## read 3mer
		#my $path="/Users/nazhu/server";
		my $gene="";
		my $range="";
		my $ftranscript=$_[0];
		my $fout=$_[1];
		print "$ftranscript\t$fout\n";
#	exit;	
		my %gene_rate=();
	   open my $OUT, ">$fout";
		print $OUT "#CHROM\tPOS\tREF\tALT\tSTRAND\tCODON\tCODONChange\tAA\tAAChange\tTranscriptID\tExon\tExonicFunc\tMutationRate\t";
	   open IN, "zless $ftranscript|" or die $ftranscript." cannot find!\n";
	   my @lastgenes=();
		my @lastpos=();
		#my @lastintervals=();
		my $line=<IN>;
	    my $strseq="";
		my @exon_starts=();
        my @exon_ends=();
		my $strand="+";
		do{
			if(!$line){last;}
			#print $line."\t";
		#	exit;
			$line=~ s/range=chr//g;
			$line=~ s/strand=//g;
			my @sets=split(/\s+/,$line);
			if(@sets<6){die "uncomplete info in $line ";}
			my @pos_info=split(/:|-/,$sets[1]);
			my @gene_info=split(/_|\./,$sets[0]);
			if(@gene_info==5){$gene_info[3]="";}
			my $chr=$pos_info[0];
		#	print $pos_info[1]."\n";
		 #  print $line."\n";
		#	print "start:".join("\t",@exon_starts)."\n";
		#	print "end:".join("\t",@exon_ends)."\n";
			if($strseq ne "" && $gene_info[2] ne $lastgenes[2]){
				### return the mutation type, it is not gene based, it is position based
				print "start:".$lastgenes[2]."\t".join("\t",@exon_starts)."\n";
				print "end:".$lastgenes[2]."\t".join("\t",@exon_ends)."\n";
				callinfo($strseq,$lastgenes[2],\@lastpos, \@exon_starts,\@exon_ends,$strand,$OUT);		
			
				$strseq="";
				@lastgenes=();
				@lastpos=();
				@exon_starts=();
				@exon_ends=();
				push(@exon_starts,$pos_info[1]);
				push(@exon_ends,$pos_info[2]);
			}else{
			
				push(@exon_starts,$pos_info[1]);
				push(@exon_ends,$pos_info[2]);
				if($strseq ne ""){
					if($gene_info[1] =~ /\d+/ && $gene_info[2] =~ /\d+/ && $gene_info[4] =~ /\d/){
						$pos_info[1]=$pos_info[1] > $lastpos[1] ? $lastpos[1]:$pos_info[1];
						$pos_info[2]=$pos_info[2] < $lastpos[2] ? $lastpos[2]:$pos_info[2];
						$gene_info[4]=$gene_info[4] < $lastgenes[4] ? $lastgenes[4]:$gene_info[4];
			
				    }
			    }
			 }
			 
		
			$strand=$sets[4];   	
			## read next line until start with :\>:
			$line=<IN>;
			my $exonseq="";
			while($line&& !($line =~ /^>/)){
				chomp($line);
				$exonseq=$exonseq."".$line;
 				#print length($line);
				$line=<IN>;
			}
			$strseq.=$exonseq;
			#print @exon_starts."\t".length($exonseq)."\n";
			#print $exonseq."\n"
			@lastgenes=@gene_info;	
			@lastpos=@pos_info;
			if($lastpos[1]>$lastpos[2]){ 
				my $temp=shift(@lastpos); 
				push(@lastpos,$temp); 
			}
		
		#	exit;
		} while($line);
		my $rs=callinfo($strseq,$lastgenes[2],\@lastpos, \@exon_starts,\@exon_ends,$strand,$OUT);
		close IN;
	   close $OUT;
		
	}
	
	
	
	sub load_info {
				
		### transcript 
		
		## load amino acide change
		if (! -e $file_rate) {print "$file_rate cannot find\t";exit;}	
		open IN, $file_rate or die $file_rate." cannot find!\n";
		while(my $line=<IN>){
			if(!($line =~ /^[ACGT]/)){next;}
			chomp($line);
			my @sets=split(/\s+/,$line);
			if(@sets<3){die "incomplete information in $line ";}
			$trichanges{$sets[0].">".$sets[1]}=$sets[2];
		}
		close IN;
		
		open IN, $file_aminoacid or die $file_rate." cannot find!\n";
		while(my $line=<IN>){
			if(!($line =~ /^[ACGT]{3}/)) {next;}
			chomp($line);
			my @sets=split(/\s+/,$line);
			if(@sets<4){die "uncomplete line : $line "}
			$amino_tables{$sets[0]}=$sets[3];
		}
		close IN;
	}	
	
	sub callinfo {
		### $strseq,$lastgenes[2],\@lastpos, \@exon_starts,\@exon_ends
		##############################
		my @sequence=split("",$_[0]);
		my $gene=$_[1];
	  # my $chrom=$_[2];
		my @posinfo=@{$_[2]};
	
		my @exon_starts=@{$_[3]};
		my @exon_ends=@{$_[4]};
		
		my $strand=$_[5];
      my $fileHandle=$_[6];
	#	my $region=$posinfo[0].":".$posinfo[1]."-".$posinfo[2];
	   if($sequence[0] eq "n"){return "";}
		my $pre="";
	    my $cur="";
		my $post="";
	   my $exon_cnt=0;
		my $protein="";
      if(@posinfo <3){
      
			print $gene."\t".join(":",@posinfo)."\n";
		#	exit;
      }
		my $chrom=$posinfo[0];
		
		my $coordinate=$posinfo[1];
	#	exit;
		if($strand eq "-"){
			$coordinate=$posinfo[2];
	#		print "$strand"."\t"."$coordinate"."\t".$posinfo[2]."\t".$posinfo[1]."\n";
		}
		
	#print "$strand"."\t"."$coordinate"."\t".$posinfo[1]."\t".$posinfo[2]."\n";
		
	#print join("::",@exon_starts)."\n";
	#print join("::",@exon_ends)."\n";		
		my @intron=();
		
		my $pointer_exon=0;
	#	if($strand eq "-") 
		
		### check whether the mappability is 1 in this region
		## if return several lines, check whether all of them are 1,
		## if yes.....
		## otherwise: check position wheher in the low mappability region ....
		#print $gene."\t".$chrom."\nstart:".join("\t",@exon_starts)."\nend:".join("\t",@exon_ends)."\t".$strand."\n";

			#continue;
		my $len=0;
		for(my $j=0;$j<@exon_starts;$j++){
				$len=$exon_ends[$j]-$exon_starts[$j]+1 + $len;
		}
	#	print $gene."\t".$len."\t".@sequence."\n";
		
		
		my $BIM=0;
		for(my $i=1;$i<@sequence-3;$i++){
			if($strand eq "-"){
				$coordinate-=1;
				if($coordinate < $exon_starts[$pointer_exon]){
					$pointer_exon+=1;
					if($pointer_exon < @exon_ends ){
						$coordinate=$exon_ends[$pointer_exon];
					}else{
						print $coordinate."\tlast exon: ".$exon_starts[$pointer_exon-1]."\n";
					}
				}
			}else{
			   $coordinate+=1;
			   if($coordinate > $exon_ends[$pointer_exon]){
				   $pointer_exon+=1;
				   #$coordinate=$exon_starts[$pointer_exon];
					
					if($pointer_exon < @exon_ends ){
						$coordinate=$exon_starts[$pointer_exon];
					}else{
						print $coordinate."\tlast exon: ".$exon_ends[$pointer_exon-1]."\n";
					}
			   }	
			}
		#	print "then $gene:coordinate: $pointer_exon:".@exon_starts."\t".$exon_starts[$pointer_exon]."\t".$exon_ends[$pointer_exon]  ."\n";
			### only check the uppercase, as those bases are in exon.
			if($sequence[$i] eq (uc $sequence[$i])  ){
				#print "$coordinate\t".$sequence[$i]."\textron\n";
				$exon_cnt=$exon_cnt+1;
			}
			else{
		#		print "$coordinate\t".$sequence[$i]."\tintron\n";
			### get the bases in the introns, it consists of donor and acceptor of splicing sites.
				next;
			}
		#	print "then2 $gene:coordinate: $pointer_exon:".@exon_starts."\t".$exon_starts[$pointer_exon]."\t".$exon_ends[$pointer_exon]  ."\n";
			if($exon_cnt%3==1){
				$protein="";
				my $k=$i;
				while(length($protein)<3 && $k<@sequence){
					if($sequence[$k] eq (uc $sequence[$k])){
				    	$protein.=$sequence[$k];
					}
					$k=$k+1;
			    }
			}
			$pre=uc $sequence[$i-1];
			$cur=$sequence[$i];
			$post=uc $sequence[$i+1];
			
			
			foreach my $nc(@nucleotides){
		       if($nc eq $sequence[$i]){
		       	   next;
		       }else{
				   my $crate=0;
				   if(exists($trichanges{$pre.$cur.$post.">".$pre.$nc.$post})){
					    $crate=$trichanges{$pre.$cur.$post.">".$pre.$nc.$post};
				   }else{
					   print "Error2: $gene\t".$pre.$cur.$post ." > ". $pre.$nc.$post." does not exist!\n";
				      # return ;
				   }
				#   print $coordinate."\t".$nc." < ".$cur.":\t".$pre.$cur.$post.">".$pre.$nc.$post."\t".$crate."\n";
				   my $cprotein=uc $nc."".(substr $protein,1);
				  
				   if($exon_cnt%3==0){
				   		$cprotein=(substr $protein,0,2).$nc;
				   }elsif($exon_cnt%3==2){
				   		$cprotein=(substr $protein,0,1).$nc.(substr $protein,2);
				   }
				   if(length($protein)<3){next;}
				   if(!exists($amino_tables{$protein})||!exists($amino_tables{$cprotein})){
					   print  "error3 :$gene\t\t$i\t$protein\t $cprotein\n";
					#   return ;
				   }
				  # print  "Amino Acid change \t$protein\t $cprotein\n";
				 #  print $amino_tables{$protein}." > > > ".$amino_tables{$cprotein}."\n";
				 
	 			#print "final $gene:coordinate:$coordinate: $pointer_exon:".@exon_starts."\t".$exon_starts[$pointer_exon]."\t".$exon_ends[$pointer_exon]  ."\n";
			
				   if($amino_tables{$protein} ne $amino_tables{$cprotein}){
				   	   if($amino_tables{$protein} eq "X"){
								print $fileHandle  $chrom."\t".$coordinate."\t".$cur."\t".$nc."\t".$strand."\t".$protein."\t".$cprotein."\t".$amino_tables{$protein}."\t".$amino_tables{$cprotein}."\t".$gene."\t".$pointer_exon."\tstoploss\t".$crate."\n";
			#			   print "stoploss\n";
					   }
				   	   if($amino_tables{$cprotein} eq "X"){
						   ## nonsense
							print $fileHandle  $chrom."\t".$coordinate."\t".$cur."\t".$nc."\t".$strand."\t".$protein."\t".$cprotein."\t".$amino_tables{$protein}."\t".$amino_tables{$cprotein}."\t".$gene."\t".$pointer_exon."\tnonsense\t".$crate."\n";
							
							## chrom, coordinate, ref,alt,$strand,$protein,$cprotein,$amino_tables{$protein},$amino_tables{$cprotein},nonsense,crate;
		#				   print "Nonsense\n";
				   	   }else{
						   if($amino_tables{$protein} ne "X"){
						   ## missense 
						  	 #$p_mis=$p_mis+$crate;
 							print $fileHandle  $chrom."\t".$coordinate."\t".$cur."\t".$nc."\t".$strand."\t".$protein."\t".$cprotein."\t".$amino_tables{$protein}."\t".$amino_tables{$cprotein}."\t".$gene."\t".$pointer_exon."\tmissense\t".$crate."\n";
						
							 	## chrom, coordinate, ref,alt,$strand,$protein,$cprotein,$amino_tables{$protein},$amino_tables{$cprotein},missense,crate;
							 ### check whether the nucletides change in mcap or meta list
				
					      }
				   	   }
				   }else{
						print $fileHandle  $chrom."\t".$coordinate."\t".$cur."\t".$nc."\t".$strand."\t".$protein."\t".$cprotein."\t".$amino_tables{$protein}."\t".$amino_tables{$cprotein}."\t".$gene."\t".$pointer_exon."\tsynonymous\t".$crate."\n";
					
					   ## syn
					#   print "Synonymous: ";
					#print chrom, coordinate, ref,alt,$strand,$protein,$cprotein,$amino_tables{$protein},$amino_tables{$cprotein},synonymous,crate;
	
				   }
		     #  print $crate."\t".$cprotein."\t".$protein."\t".$pre.$cur.$post.">".$pre.$nc.$post."\n";
			   }
			  
			}
		
		
		}

	  
	  
	}
	
	

print "Input: perl variants_generate_forJiaYao.pl  transcript outputFile\n";
#	print "Input perl Mutation_rate_20180814.pl transcript_file outname f\n";
	#print "#filter_with_mappability: 0 or 1, 0 means no mappablity filter, otherwise filter the position with mappablity\n";
my $ftranscript=$ARGV[0]; #"$path/Resources/Transcript/GencodeV19_splice.gz";MY $QMAP=1;
my $fout=$ARGV[1];	
my $QMAP=0;
load_info();

main($ftranscript,$fout);
print "output to $fout\n"
