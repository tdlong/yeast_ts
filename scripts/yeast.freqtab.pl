my $Nfounders = 18;
while (my $line = <STDIN>){
	chomp $line;
	my @Line = split(' ',$line);
	my $N = scalar @Line;
	my $outline = "$Line[0]\t$Line[1]";
	# biallelic
	$mygood="OK";
	if (!($Line[2] =~ m/^[ACGT]$/ && $Line[3] =~ m/^[ACGT]$/)){
		$mygood=$mygood."NTBI";
		}
	my $foundcount=1;
	my $sumfreq=0;
	for(my $i=4; $i<$N-1; $i=$i+2){
		my $T = $Line[$i] + $Line[$i+1];
		my $freq = "NA";
		if($T>0){$freq = $Line[$i] / $T;}
		# per sample coverage
		if ($T<15 && $foundcount <= $Nfounders){
			$mygood=$mygood."WTOT";
#			print "$Line[1]\t$foundcount\t$T\n"
			}
		# per sample heterozygosity (founders only)
		my $het=2*$freq*(1-$freq);
		if ($foundcount <= $Nfounders && $het > 0.05){
			$mygood=$mygood."WHET";
#			print "$Line[1]\t$foundcount\t$het\n"
			}
		# sum of frequencies over founders
		if ($foundcount <= $Nfounders){$sumfreq=$sumfreq+$freq;}
		my $temp = sprintf("%.4f", $freq);
		$outline = $outline."\t$temp";
		$foundcount=$foundcount+1;
		}
	if($mygood eq "OK" && $sumfreq > 0.98 && $sumfreq < $Nfounders - 0.98){
		print "$outline\n";
		}
	}

