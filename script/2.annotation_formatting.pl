use strict;
use warnings;
use Getopt::Long;

my $codingGeneFileName = "codingGene.tab";
my $noncodingGeneFileName = "noncodingGene.tab";

my $outputHeader = "";
my $prodigalFile = "";
my $rnammerFile = "";
my $trnascanFile = "";
my $rfamRf = "";
my $rfamFile = "";
my $cogPepID2CogID = "";
my $cogCogID2ClassID = "";
my $ClassID2FirstID = "";
my $cogFile = "";
my $keggFile = "";
my $keggko = "";
my $keggmap ="";
my $keggclass =" ";
my $str;
my $g;
my %codingGene;
my %noncodingGene;
my @ae;
my @subae;
my @subsubae;
my $tempPos;
my $flagOfCOG = 0;
my $flagOfKEGG = 0;
my $flagOfGO = 0;

GetOptions(
    'outputHeader|o=s' => \$outputHeader,
    'prodigalFile=s' => \$prodigalFile,
    'rnammerFile=s' => \$rnammerFile,
    'trnascanFile=s' => \$trnascanFile,
    'rfamRf=s' => \$rfamRf,
    'rfamFile=s' => \$rfamFile,
    'cogPepID2CogID=s' => \$cogPepID2CogID,
    'cogCogID2ClassID=s' => \$cogCogID2ClassID,
    'ClassID2FirstID=s' => \$ClassID2FirstID,
    'cogFile=s' => \$cogFile,
    'keggFile=s' => \$keggFile,
    'keggko=s' => \$keggko,
    'keggmap=s' => \$keggmap,
    'keggclass=s' => \$keggclass,
);

open FIN, $prodigalFile or die "No Prodigal result!\n"; ###gff
#print "$prodigalFile\n";
while (defined($str=<FIN>)){
	chomp $str;
	next if($str=~/^#/);
	#print "$str\n";
	my @ae=split /\s+/,$str;
	if($str=~/=/){
		if($ae[2] eq "CDS"){
			print "$str\n";
			my @be=split /;/,$ae[-1];
			$be[0]=~s/ID=//g;
		        $codingGene{$be[0]}{"seqID"} = $ae[0];
	        	$codingGene{$be[0]}{"source"} = "Prodigal";
	        	$codingGene{$be[0]}{"start"} = $ae[3];
		        $codingGene{$be[0]}{"end"} = $ae[4];
	        	$codingGene{$be[0]}{"strand"} = $ae[6];
	        	$codingGene{$be[0]}{"partial"} = -1;
		}
	}
}
close FIN;

##################Eggnog
my %ko2map;
my %map2class;
my %koclass2first;

if(( -e $keggFile) && (-e $keggko) && ( -e $keggmap) && ( -e $keggclass)){
	$flagOfKEGG=1;
}
print "$flagOfKEGG\n";
if($flagOfKEGG){
	open FIN,"$keggko";
	while(defined($str=<FIN>)){
		chomp $str;
		if($str=~/map/){
			my @arr=split /\t/,$str;
			$arr[0]=~s/path:map//g;
			$arr[1]=~s/ko://g;
			push @{$ko2map{$arr[1]}},($arr[0]);
		#	print "$arr[1]\t$arr[0]\n";
		}
	}
	close FIN;
	open FIN,"$keggmap";
	my $secondary_classify;
	while(defined($str=<FIN>)){
		chomp $str;
		if($str=~/##/){
			$secondary_classify = substr($str,2);
		}else{
			my @ae = split /\t/,$str;
			$map2class{$ae[0]}=$secondary_classify;
		}
	}
	close FIN;
	
	open FIN,"$keggclass";
	while(defined($str=<FIN>)){
		chomp $str;
		my @ae=split /\t/,$str;
		$koclass2first{$ae[0]}=$ae[1];
	}
	close FIN;

	open FIN,"$keggFile";
	while(defined($str=<FIN>)){
		chomp $str;
		next if($str=~/^#/);
		next if($str=~/identity/);
		my @ae=split /\t/,$str;
		my @aa=split /\|/,$ae[2];
		#print "$ae[2]\t$aa[-1]\n";
		if($aa[-1]=~/K/){
			$aa[-1]=~s/\s//g;
			#print "$aa[-1]\n";
			my $st;
			my %tmpst;
			my @be=split /\//,$aa[-1];
			for my $i(@be){
				#print "$i\n";
				if(exists $ko2map{$i}){
					#print "$i\n";
					my @line=@{$ko2map{$i}};
					#print "$i\t@line\n";
					for my $j(@line){
						if(exists $map2class{$j} && exists $koclass2first{$map2class{$j}}){
							$tmpst{$koclass2first{$map2class{$j}}}="";
							#$st.="$koclass2first{$map2class{$j}};";
						}
					}
				}
			}
			for my $k(keys %tmpst){
				$st.="$k;";
			}
			if($st){
				$st=~s/;$//g;
				$codingGene{$ae[0]}{"KO"} ="$aa[-1] | $st";
				$codingGene{$ae[0]}{"KOcolor"} =$st;
			}
		}
	}
	close FIN;
}
my %cogPepCog;
my %cogClass;
my %cogFirst;
if ((-e $cogPepID2CogID) && (-e $cogCogID2ClassID) && (-e $cogFile) && ( -e $ClassID2FirstID)) { $flagOfCOG = 1; }
if ($flagOfCOG) { 
    open FIN, $cogPepID2CogID;
    while (defined($str=<FIN>)){
        chomp $str;
        @ae = split /,/, $str;
        $cogPepCog{$ae[2]} = $ae[6]; ## GENE ID 	COG
    }
    close FIN;
    
    open FIN, $cogCogID2ClassID;
    for (1) { $str = <FIN>; }
    while (defined($str=<FIN>)){
        chomp $str;
        @ae = split /\t/, $str;
        $cogClass{$ae[0]} = $ae[1]; ##COG type
    }
    close FIN;
    
    open FIN,"$ClassID2FirstID";
    while(defined($str=<FIN>)){
	chomp $str;
	@ae = split /\t/,$str;
	$cogFirst{$ae[0]}=$ae[2]; ##type leve1
    }
    close FIN;

    
    open FIN, $cogFile;
    while (defined($str=<FIN>)){
        chomp $str;
	next if($str=~/^#/);
        @ae = split /\t/, $str;
	my $se;
	if(exists $cogPepCog{$ae[1]} && exists $cogClass{$cogPepCog{$ae[1]}}){
		my $letter=$cogClass{$cogPepCog{$ae[1]}};
		my @line=split //,$letter;
		for my $i(@line){
			$se.="$cogFirst{$i};";
		}
	}
	if($se){
		$se=~s/;$//g;
        	$codingGene{$ae[0]}{"COG"} = $ae[1]." | ".$cogPepCog{$ae[1]}." | ".$cogClass{$cogPepCog{$ae[1]}}." | ".$se." | ".$ae[2];
	        $codingGene{$ae[0]}{"COGcolor"} = $se;
	}
    }
	close FIN;
}

open FOUT, ">".$outputHeader.".".$codingGeneFileName;
print FOUT "SeqID\tSource\tStart\tEnd\tStrand\tPartial";
if ($flagOfCOG) { print FOUT "\tCOG\tCOG category"; }
if($flagOfKEGG) { print FOUT"\tKO\tKO category";}
print FOUT "\n";
foreach $g (sort {$codingGene{$b}{"seqID"} cmp $codingGene{$a}{"seqID"} or $codingGene{$a}{"start"}<=>$codingGene{$b}{"start"}} keys %codingGene) {
   print FOUT $codingGene{$g}{"seqID"}, "\t";
    print FOUT $codingGene{$g}{"source"}, "\t";
    print FOUT $codingGene{$g}{"start"}, "\t";
    print FOUT $codingGene{$g}{"end"}, "\t";
    print FOUT $codingGene{$g}{"strand"}, "\t";
    print FOUT $codingGene{$g}{"partial"};
    if ($flagOfCOG) {
        if (exists($codingGene{$g}{"COG"})) {
            print FOUT "\t", $codingGene{$g}{"COG"}, "\t", $codingGene{$g}{"COGcolor"};
        } else {
            print FOUT "\t#\t#";
        }
    }
	if($flagOfKEGG){
		if(exists $codingGene{$g}{"KO"}){
			print FOUT"\t",$codingGene{$g}{"KO"},"\t",$codingGene{$g}{"KOcolor"};
		}else{
			print FOUT"\t#\t#";
		}
	}
    print FOUT "\n";
}
close FOUT;

if (-e $rnammerFile) {
    open FIN, $rnammerFile;
    while (defined($str=<FIN>)){
        chomp $str;
        if ($str!~/^#/) {
            @ae = split /\t/, $str;
            $ae[8] =~ s/s_/S /;
            $g = $ae[0]."_".$ae[3]."_".$ae[4];
            $noncodingGene{$g}{"seqID"} = $ae[0];
            $noncodingGene{$g}{"source"} = "RNAmmer";
            $noncodingGene{$g}{"start"} = $ae[3];
            $noncodingGene{$g}{"end"} = $ae[4];
            $noncodingGene{$g}{"strand"} = $ae[6];
            $noncodingGene{$g}{"type"} = "rRNA";
            $noncodingGene{$g}{"annotation"} = $ae[8];
        }
    }
    close FIN;
}

if (-e $trnascanFile) {
    open FIN, $trnascanFile;
    for (1..3) { $str=<FIN>; }
    while (defined($str=<FIN>)){
        chomp $str;
        @ae = split /\t/, $str;
        foreach (0..scalar(@ae)-1) {
            $ae[$_] =~ s/^\s+|\s+$//g;
        }
        if ($ae[3]>$ae[2]) {
            $g = $ae[0]."_".$ae[2]."_".$ae[3];
            $noncodingGene{$g}{"start"} = $ae[2];
            $noncodingGene{$g}{"end"} = $ae[3];
            $noncodingGene{$g}{"strand"} = "+";
        } else {
            $g = $ae[0]."_".$ae[3]."_".$ae[2];
            $noncodingGene{$g}{"start"} = $ae[3];
            $noncodingGene{$g}{"end"} = $ae[2];
            $noncodingGene{$g}{"strand"} = "-";
        }
        $noncodingGene{$g}{"seqID"} = $ae[0];
        $noncodingGene{$g}{"source"} = "tRNAscan-SE";
        $noncodingGene{$g}{"type"} = "tRNA";
        $noncodingGene{$g}{"annotation"} = $ae[4]." tRNA";
    }
    close FIN;
}

my %rfamClass;
my %rfamAnnotation;
if ((-e $rfamRf) && (-e $rfamFile)) {
    open FIN, $rfamRf;
    while (defined($str=<FIN>)){
        chomp $str;
        @ae = split /\t/, $str;
        $ae[0] =~ s/^\s+|\s+$//g;
        $ae[1] =~ s/^\s+|\s+$//g;
        $ae[-1] =~ s/^\s+|\s+$//g;
        $rfamClass{$ae[0]} = $ae[-1];
        $rfamAnnotation{$ae[0]} = $ae[1];
    }
    close FIN;
    
    open FIN, $rfamFile;
    for (1) { <FIN>; }
    while (defined($str=<FIN>)){
        chomp $str;
        if ($str!~/^#/) {
            @ae = split /\s+/, $str;
            if ($rfamClass{$ae[2]} eq "sRNA" || $rfamClass{$ae[2]} eq "CRISPR") {
                if ($ae[10]>$ae[9]) {
                    $g = $ae[3]."_".$ae[9]."_".$ae[10];
                    $noncodingGene{$g}{"start"} = $ae[9];
                    $noncodingGene{$g}{"end"} = $ae[10];
                } else {
                    $g = $ae[3]."_".$ae[10]."_".$ae[9];
                    $noncodingGene{$g}{"start"} = $ae[10];
                    $noncodingGene{$g}{"end"} = $ae[9];
                }
                $noncodingGene{$g}{"seqID"} = $ae[3];
                $noncodingGene{$g}{"source"} = "RFAM";
                $noncodingGene{$g}{"strand"} = $ae[11];
                $noncodingGene{$g}{"type"} = $rfamClass{$ae[2]};
                $noncodingGene{$g}{"annotation"} = $rfamAnnotation{$ae[2]};
            }
        }
    }
}

open FOUT, ">".$outputHeader.".".$noncodingGeneFileName;
print FOUT "SeqID\tSource\tStart\tEnd\tStrand\tType\tAnnotation";
print FOUT "\n";
foreach $g (sort {$noncodingGene{$b}{"seqID"} cmp $noncodingGene{$a}{"seqID"} or $noncodingGene{$a}{"start"}<=>$noncodingGene{$b}{"start"}} keys %noncodingGene) {
    print FOUT $noncodingGene{$g}{"seqID"}, "\t";
    print FOUT $noncodingGene{$g}{"source"}, "\t";
    print FOUT $noncodingGene{$g}{"start"}, "\t";
    print FOUT $noncodingGene{$g}{"end"}, "\t";
    print FOUT $noncodingGene{$g}{"strand"}, "\t";
    print FOUT $noncodingGene{$g}{"type"}, "\t";
    if (exists($noncodingGene{$g}{"annotation"})) { print FOUT $noncodingGene{$g}{"annotation"}, "\n"; } else { print FOUT "#\n"; }
}
close FOUT;
