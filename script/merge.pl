my $sam=$ARGV[0];
my $path=$ARGV[1];
open IN,"$path/$ARGV[0].genomeLength.tab";##GCM10030259.genomeLength.tab
my %hash;
my $strat=1;
my $end;
my $len;
my %sort;
my $n=0;
while(<IN>){
	chomp;
	next if(/SeqID/);
	my @ae=split /\t/,$_;
	$end=$strat+$ae[1]-1;
	$n++;
	$sort{$n}=$ae[0];
	$hash{$ae[0]}->[0]=$strat;
	$hash{$ae[0]}->[1]=$end;
#	print "$ae[0]\t$strat\t$end\n";
	$strat=$end+1;
	$len+=$ae[1];
}
close IN;
system("mkdir -p $path/merge");
open OUT,">$path/merge/$sam.genomeLength.tab";
print OUT"SeqID\tlength\n";
print OUT"Contig1\t$len\n";

open FA,"$ARGV[2]";
my $id;
my %fa;
while(<FA>){
	chomp;
	if(/^>(\S+)/){
		$id=$1;
		
	}else{
		$fa{$id}.=$_;
	}
}
close FA;
open OFA,">$path/merge/$sam.fasta";
print OFA">Contig1\n";
for my $k(sort {$a<=>$b} keys %sort){
	my $kk=$sort{$k};
	if(exists $fa{$kk}){
		my $seq=$fa{$kk};
		print OFA"$seq";
	}
}
print FA"\n";

open NON,"$path/$sam.noncodingGene.tab"; ##GCM10030259.noncodingGene.tab
open ONON,">$path/merge/$sam.noncodingGene.tab";
while(<NON>){
	chomp;
	if(/SeqID/){
		print ONON"$_\n";
	}else{
		my @ae=split /\t/,$_;
		my $st=$hash{$ae[0]}->[0];
		$ae[2]=$ae[2]+$st-1;
		$ae[3]=$ae[3]+$st-1;
		$ae[0]="Contig1";
		my $seq=join("\t",@ae);
		print ONON"$seq\n";
	}
}
close NON;

open CO,"$path/$sam.codingGene.tab";## GCM10030259.codingGene.tab
open OCO,">$path/merge/$sam.codingGene.tab";
while(<CO>){
	chomp;
	if(/SeqID/){
		print OCO"$_\n";
	}else{
		my @ae=split /\t/,$_;
		my $st=$hash{$ae[0]}->[0];
		$ae[2]=$ae[2]+$st-1;
                $ae[3]=$ae[3]+$st-1;
		$ae[0]="Contig1";
                my $seq=join("\t",@ae);	
		print OCO"$seq\n";
	}
}
		
