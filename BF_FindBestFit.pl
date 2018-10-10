## This script is used to compare recombination data simulated using MADpatterns (https://mwhite4.github.io/MADpatterns/) to experimental recombination datasets
## Specifically, the output files of "analyze_events_on_linear_objects" run on recombination data 
## simulated using "crossover_simulation" are compared to output file of "analyze_events_on_linear_objects"
## run on the experimental data, and the Best Fit simulation identified.
##
## The script requires the user to provide
##  1. The directory containing the simulated and experimental analysis files 
##  2. The file name of the experimental analysis file
##  3. Whether the experimental data are from bivalents (0) or gametes (1) e.g. cytological data or backcross data)
##
##
##
##
##


use strict;
use warnings;

my $wd = $ARGV[0]; # directory containing the output files from Matlab program "analyze_events_on_linear_objects"

my $exDATAfile = $ARGV[1]; # file containing output of "analyze_events_on_linear_objects" run on the experimental dataset

my $gamete = 1;  # 1 = experimental data is from gametes; 0 = experimental data is from bivalents

if(exists $ARGV[2]){
	$gamete = $ARGV[2];
}


my @BFdata = `ls $wd/*.csv_line*.csv | grep intervals` ; # creates array of the filenames of the analysis files of the BF simulations
my %BFscore;

##### Get CO/ED/CoC info from BF analysis of experimental data
my @exArray;
my @exArray2;
my %simHash;

open (EXD, "$wd/$exDATAfile"); #open file

my $n = 0;

while(<EXD>){ #creates array of arrays from the experimental input file i.e. and array of the columns of the analysis of the experimental dataset.
	chomp;
	if($n == 0){	
	$n++;
		next;
	} else {
		my @line = split (/[,]/, $_) ;
		my $m = 0;
		foreach (@line){
			chomp;
			if($_ =~ m/NaN/){
			#do nothing
			} else {
				push (@{$exArray[$m]}, $_);
			}
		$m++;	
		}
	}
	$n++;	
}

close EXD;

#### creates smaller array containing experimental CO ED and CoC data;
@{$exArray2[0]} = @{$exArray[6]}; #CO frequency analysis
@{$exArray2[1]} = @{$exArray[4]}; #ED analysis
@{$exArray2[2]} = @{$exArray[8]}; #CoC analysis


##### Get CO/ED/CoC info from BF analysis of simulated data
foreach(@BFdata){
	chomp;
	my @BFarray;
	my @split = split (/\//, $_) ;
	my $BFdataFile = $split[$#split] ; 
	@split = split (/[y.]/, $BFdataFile) ;
	my $sim = "$split[1]".".$split[2]";
	open (BFD, "$_");
	$n =0;
	while(<BFD>){
		chomp;
		if($n == 0){
		$n++;
		next;
		} elsif ($n == 1){
			my @line = split (/[,]/, $_) ;
			my $m = 0;
			foreach (@line){
				chomp;
				if($_ =~ m/NaN/){
					push (@{$BFarray[$m]}, 0);
				} else {
					push (@{$BFarray[$m]}, $_);
				}
			$m++;
			}	
		}elsif ($n == 22){ #This number set the maximum lines read of each file and minimises unnecessary looping through long files. This number should be increased if analysing more than 20 intervals per bivalent.
			last;
		} else {
			my @line = split (/[,]/, $_) ;
			my $m = 0;
			foreach (@line){
				chomp;
				if($_ =~ m/NaN/){
				#do nothing
				} else {
					push (@{$BFarray[$m]}, $_);
				}
			$m++;	
			}
		}
		$n++;	
	}
	close BFD;
	@{$simHash{$sim}} = @BFarray;	
}	

my %COs; #hash containing CO freqs for each sim: key = sim
my %ED; #hash containing ED for each sim: key = sim
my %CoC; #hash containing CoC intervals for each sim: key = sim

my %COsScore; #hash containing score for COs for each sim: key = sim
my %EDScore; #hash containing score for ED for each sim: key = sim
my %CoCScore; #hash containing score for CoC intervals for each sim: key = sim
my %totalScore; #


foreach my $key (keys %simHash) {
	my @array = @{$simHash{$key}};
	print "$key\n";
	if($gamete == 1){ # inputs CO ED and CoC data and tranforms for gamete data: if bivalent data i.e. $gamete = 0 then no transformation performed.
    	
    	my @temparray = @{$array[6]};
    	foreach my $x (@temparray) {
        	$x = $x / 2;
        } # transforms COs for gamete data
    	
    	@{$COs{$key}} = @temparray; #inputs COs
    	
    	@{$CoC{$key}} = @{$array[8]}; #inputs CoC for sim into %CoC (no transform required)
        
##### Inputs ED for sim to %ED and transforms for gamete data #####
    	my $l = scalar @{$array[4]};
    	my $b = 0; #COs bivalent
    	my $g = 0; #COs gamete
    	@{$ED{$key}}=(0) x $l;
    	while($b < $l){
    		$g = 0;
    		while($g <= $b){
    			${$ED{$key}}[$g] = ${$ED{$key}}[$g] + $array[4][$b] * (fact($b)/((fact($g) * (fact($b - $g)))))/(2**$b);
				$g++;
			}
		$b++;
		}
	} else { # inputs CO ED and CoC data with no transformation
		@{$COs{$key}} = @{$array[6]};
		@{$ED{$key}} = @{$array[4]};
		@{$CoC{$key}} = @{$array[8]};
	}
}


############# This section computes the score for correlation between experimental and simulated data
foreach my $key (keys %COs) {
	$COsScore{$key} = 0;
	my $m = 0;
	foreach (@{$COs{$key}}){ # Calculates score for COs
        chomp;
        my $ratio = ($_ + 0.01) / ($exArray2[0][$m] + 0.01); #transforms to avoid illegal division by 0
    	$COsScore{$key} = $COsScore{$key} + abs(log($ratio)/log(2)); #finds absolute value of the log2 of the ratio
    	$m++;
    }
}

foreach my $key (keys %CoC) { # finds sum of squares for two CoC curves
	$CoCScore{$key} = 0;
	my $m = 0;
	foreach (@{$CoC{$key}}){ # Calculates score for CoC
        chomp;
        my $ratio = ($_ + 0.01) / ($exArray2[2][$m] + 0.01); #transforms to avoid illegal division by 0
    	$CoCScore{$key} = $CoCScore{$key} + abs(log($ratio)/log(2)); #finds absolute value of the log2 of the ratio
    	$m++;
    }
}

foreach my $key (keys %ED) { # finds sum of squares of differences for ED
	$EDScore{$key} = 0;
	my $sim = scalar @{$ED{$key}};
	my $exp = scalar @{$exArray2[1]};
	my $m = 0;
	
	if($sim >= $exp){
		foreach (@{$ED{$key}}){ # Calculates score for ED when sim ED array is longer
        	chomp;
        	#print "sim:$_\texp:$exArray2[1][$m]\n";
        	if ($m < $exp){
    			$EDScore{$key} = $EDScore{$key} + (($_ - $exArray2[1][$m])**2);
    		} else {
    			$EDScore{$key} = $EDScore{$key} + ($_**2);
    			#print "sim:$_\texp:0\n";
    		}
    	$m++;
    	}
    } else {
    	foreach (@{$ED{$key}}){ #  # Calculates score for ED when exp ED array is longer
        	chomp;
        	if ($m < $sim){
    			$EDScore{$key} = $EDScore{$key} + ($_ - $exArray2[1][$m])**2;
    			#print "sim:$_\texp:$exArray2[1][$m]\n";
    		} else {
    			$EDScore{$key} = $EDScore{$key} + ($exArray2[1][$m]**2);
    			#print "sim:0\texp:$exArray2[1][$m]\n";
    		}
    	$m++;
    	}
    }
}

$n = 0;

####### This section sums the scores for CO ED and CoC to give totalScore

open OUT, ">$wd/BestFitAnalysis.txt" ; #creating text file to house script
       		
foreach my $key (sort { $COsScore{$a} <=> $COsScore{$b} } keys %COsScore) {
	$totalScore{$key} = $n;
	$n++;
}

$n = 0;

foreach my $key (sort { $CoCScore{$a} <=> $CoCScore{$b} } keys %CoCScore) {
	$totalScore{$key} = $totalScore{$key} + $n;
	$n++;
}		

$n = 0;

foreach my $key (sort { $EDScore{$a} <=> $EDScore{$b} } keys %EDScore) {
	$totalScore{$key} =  $totalScore{$key} + $n;
	$n++;
}


######## Prints sims in order of best (lowest) score

print OUT "\n";
print OUT "Total Score\n";
#print "Total Score\n";

$n = 0;

my $BestSim;

foreach my $key (sort { $totalScore{$a} <=> $totalScore{$b} } keys %totalScore) {
	if ($n == 0){
		$BestSim = $key;
	}
	print OUT "sim_$key\t$totalScore{$key}\n";
	#print "sim_$key\t$totalScore{$key}\n";
	$n++;
}
print OUT "\nCOs Score\n";
#print "\nCOs Score\n";

foreach my $key (sort { $COsScore{$a} <=> $COsScore{$b} } keys %COsScore) {
	#print "SIM:\t$key\t$COsScore{$key}\n";
	print OUT "sim_$key\t$COsScore{$key}\n";
}

print OUT "\nCoC Score\n";
#print "\nCoC Score\n";

foreach my $key (sort { $CoCScore{$a} <=> $CoCScore{$b} } keys %CoCScore) {
	#print "SIM:\t$key\t$CoCScore{$key}\n";
	print OUT "sim_$key\t$CoCScore{$key}\n";
}

print OUT "\nED Score\n";
#print "\nED Score\n";

foreach my $key (sort { $EDScore{$a} <=> $EDScore{$b} } keys %EDScore) {
	print "SIM:\t$key\t$EDScore{$key}\n";
	#print OUT "sim_$key\t$EDScore{$key}\n";
}

close OUT;

sub fact {
    my $n = shift;
    my $f = 1;
    $f *= $n-- while $n > 0;    # Multiply, then decrement
    return $f;
}
