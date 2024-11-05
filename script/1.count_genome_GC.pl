use strict;
use warnings;
use Getopt::Long;

my $contigGCFileName = "genomeGC.tab";
my $contigLengthFileName = "genomeLength.tab";

my $inputFile = "";
my $outputHeader;
my $windowSize = 1000;
my $stepSize = 200;
my @ae;
my @fileList;
my $i = 0;

GetOptions(
    'inputFile|i=s' => \$inputFile,
    'outputHeader|o=s' => \$outputHeader,
    'windowSize|w=i' => \$windowSize,
    'stepSize|s=i' => \$stepSize,
);

my $flagForScan;
my $str;
my ($stpos, $enpos);
my $id = "";
my $seq;
my $subseq;
my $file;
my $t;
my %countOfBase;
my %length;
open FOUT, ">".$outputHeader.".".$contigGCFileName || die "Can not open file !\n";
print FOUT "SeqID\tstartPosition\tendPosition\tGCcontent\tGCskew\n";
open FIN, $inputFile || warn "There is no sequence file $file!\n";
do {
    $flagForScan = defined($str=<FIN>);
    if ($flagForScan) {
        chomp $str;
        $stpos = index ($str, ">");
    } else {
        $stpos = -1;
    }
    if ($stpos!=-1 || !$flagForScan) {
        if ($id ne "") {
            $length{$id} = length $seq;
            $seq = uc $seq;
            $t = 0;
            while ($t+$windowSize-1<length($seq)) {
                $subseq = substr $seq, $t, $windowSize-1;
                @ae = split //, $subseq;
                %countOfBase = ("A"=>0, "T"=>0, "C"=>0, "G"=>0);
                foreach (@ae) {
                    if ($_ eq "A" || $_ eq "T" || $_ eq "C" || $_ eq "G") { $countOfBase{$_}++; }
                    elsif ($_ eq "R") { $countOfBase{"A"}+=0.5; $countOfBase{"G"}+=0.5; }
                    elsif ($_ eq "Y") { $countOfBase{"C"}+=0.5; $countOfBase{"T"}+=0.5; }
                    elsif ($_ eq "M") { $countOfBase{"A"}+=0.5; $countOfBase{"C"}+=0.5; }
                    elsif ($_ eq "K") { $countOfBase{"G"}+=0.5; $countOfBase{"T"}+=0.5; }
                    elsif ($_ eq "S") { $countOfBase{"G"}+=0.5; $countOfBase{"C"}+=0.5; }
                    elsif ($_ eq "W") { $countOfBase{"A"}+=0.5; $countOfBase{"T"}+=0.5; }
                    elsif ($_ eq "H") { $countOfBase{"A"}+=1/3; $countOfBase{"T"}+=1/3; $countOfBase{"C"}+=1/3; }
                    elsif ($_ eq "B") { $countOfBase{"G"}+=1/3; $countOfBase{"T"}+=1/3; $countOfBase{"C"}+=1/3; }
                    elsif ($_ eq "V") { $countOfBase{"G"}+=1/3; $countOfBase{"A"}+=1/3; $countOfBase{"C"}+=1/3; }
                    elsif ($_ eq "D") { $countOfBase{"G"}+=1/3; $countOfBase{"A"}+=1/3; $countOfBase{"T"}+=1/3; }
                    elsif ($_ eq "N") { $countOfBase{"A"}+=1/4; $countOfBase{"T"}+=1/4; $countOfBase{"C"}+=1/4; $countOfBase{"G"}+=1/4; }
                }
                print FOUT $id, "\t", $t, "\t", $t+$windowSize, "\t";
		if ($countOfBase{"G"}+$countOfBase{"C"} >0){
	                printf FOUT "%.4f\t%.4f\n", (($countOfBase{"G"}+$countOfBase{"C"})/length($subseq), ($countOfBase{"G"}-$countOfBase{"C"})/($countOfBase{"G"}+$countOfBase{"C"}));
		}else{
			printf FOUT "%.4f\t%.4f\n",(0,0);
		}
                $t += $stepSize;
            }
        }
        if ($flagForScan) {
            $enpos = index($str, " ");
            if ($enpos!=-1) {
                $id = substr($str, $stpos+1, $enpos-$stpos-1);
            } else {
                $id = substr($str, $stpos+1);
            }
            $seq = "";
        }
    } else {
        chomp $str;
        $seq = $seq.$str;
    }
} while ($flagForScan);
close FOUT;

open FOUT, ">".$outputHeader.".".$contigLengthFileName || die "Can not open file!\n";
print FOUT "SeqID\tlength\n";
foreach (keys %length) {
    print FOUT $_, "\t", $length{$_}, "\n";
}
close FOUT;
