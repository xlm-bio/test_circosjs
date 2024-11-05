use strict;
use warnings;
use Getopt::Long;

my $contigGCFileName = "genomeGC.tab";
my $contigLengthFileName = "genomeLength.tab";
my $codingGeneFileName = "codingGene.tab";
my $noncodingGeneFileName = "noncodingGene.tab";

my $inputHeader = "";
my $outputXML = "";
my $CGViewFile = "";
my $COGColorFile = "";
my $contigName = "";
my $legendPosition = "";
my $defaultcolor = "gray";
my $length = 0;
my $radius = 300;
my $height = 800;
my $width = 800;
my $str;
my $flagOfGC = 1;
my $flagOfCOG = 1;
my $flagOfCDS = 1;
my $flagOfNoncodingGene = 1;
my %COGcolor = ();
my %legendColor = ();
my %header = ();
my @CGViewMain = ();
my @textUperCenter = ();
my @textLowerCenter = ();
my @COGList = ();
my @dataGCContent = ();
my @dataGCskew = ();
my @dataCDS = ();
my @dataNoncoding = ();
my @ae;
my @subae;

GetOptions(
    'inputHeader|i=s' => \$inputHeader,
    'outputXML|o=s' => \$outputXML,
    'CGViewFile=s' => \$CGViewFile,
    'COGColorFile=s' => \$COGColorFile,
    'contigName=s' => \$contigName,
);

if ($inputHeader eq "") { die "Argument -i is empty!\n"; }

open FIN, $CGViewFile || die "There is no CGView config file!\n";
while (defined($str=<FIN>)){
    chomp $str;
    @ae = split /=/, $str;
    if (scalar(@ae)==2) {
        $ae[1] =~s/"//g;
        if ($ae[0] eq "legendPosition") { $legendPosition = $ae[1]; }
        elsif ($ae[0] eq "textUperCenter") { push @textUperCenter, $ae[1]; }
        elsif ($ae[0] eq "textLowerCenter") { push @textLowerCenter, $ae[1];; }
        elsif ($ae[0] eq "GCskew+") { $legendColor{"GCskew+"} = $ae[1]; }
        elsif ($ae[0] eq "GCskew-") { $legendColor{"GCskew-"} = $ae[1]; }
        elsif ($ae[0] eq "GCcontent") { $legendColor{"GCcontent"} = $ae[1]; }
        elsif ($ae[0] eq "rRNA") { $legendColor{"rRNA"} = $ae[1]; }
        elsif ($ae[0] eq "tRNA") { $legendColor{"tRNA"} = $ae[1]; }
        elsif ($ae[0] eq "sRNA") { $legendColor{"sRNA"} = $ae[1]; }
        elsif ($ae[0] eq "CRISPR") { $legendColor{"CRISPR"} = $ae[1]; }
        elsif ($ae[0] eq "CDS") { $legendColor{"CDS"} = $ae[1]; }
        else {
            if ($ae[0] eq "backboneRadius") { $radius = $ae[1]; }
            if ($ae[0] eq "height") { $height = $ae[1]; }
            if ($ae[0] eq "width") { $width = $ae[1]; }
            push @CGViewMain, $str;
        }
    }
}
close FIN;

if (-e $COGColorFile) {
    open FIN, $COGColorFile;
    while (defined($str=<FIN>)){
        chomp $str;
        @ae = split /\t/, $str;
        push @COGList, $ae[0];
        $COGcolor{$ae[0]} = $ae[2];
    }
    close FIN;
} else {
    $flagOfCOG = 0;
}

if (!(-e $inputHeader.".".$contigLengthFileName)) { die "There is no genome length information!\n"; } else {
    open FIN, $inputHeader.".".$contigLengthFileName;
    %header = ();
    while (defined($str=<FIN>)){
        chomp $str;
        @ae = split /\t/, $str;
        if (scalar(keys %header)==0) {
            for (0..scalar(@ae)-1) { $header{$ae[$_]} = $_; }
            if (!exists($header{"SeqID"})) { die "There is no header 'SeqID' in file $contigLengthFileName\n"; }
            if (!exists($header{"length"})) { die "There is no header 'length' in file $contigLengthFileName\n"; }    
        } else {
            if ($ae[$header{"SeqID"}] eq $contigName) { $length = $ae[$header{"length"}]; }
        }
    }
    close FIN;
}
if (!(-e $inputHeader.".".$contigGCFileName)) { $flagOfGC = 0; } else {
    open FIN, $inputHeader.".".$contigGCFileName;
    %header = ();
    while (defined($str=<FIN>)){
        chomp $str;
        @ae = split /\t/, $str;
        if (scalar(keys %header)==0) {
            for (0..scalar(@ae)-1) { $header{$ae[$_]} = $_; }
            if (!exists($header{"SeqID"})) { die "There is no header 'SeqID' in file $contigGCFileName\n"; }
            if (!exists($header{"startPosition"})) { die "There is no header 'startPosition' in file $contigGCFileName\n"; } 
            if (!exists($header{"endPosition"})) { die "There is no header 'endPosition' in file $contigGCFileName\n"; } 
            if (!exists($header{"GCcontent"})) { die "There is no header 'GCcontent' in file $contigGCFileName\n"; }
            if (!exists($header{"GCskew"})) { die "There is no header 'GCskew' in file $contigGCFileName\n"; }
        } else {
            if ($ae[$header{"SeqID"}] eq $contigName) {
                push @dataGCContent, [$ae[$header{"startPosition"}], $ae[$header{"endPosition"}], $ae[$header{"GCcontent"}]];
                push @dataGCskew, [$ae[$header{"startPosition"}], $ae[$header{"endPosition"}], $ae[$header{"GCskew"}]];
            }
        }
    }
    close FIN;
}
if (!(-e $inputHeader.".".$codingGeneFileName)) { $flagOfCDS = 0; } else {
    open FIN, $inputHeader.".".$codingGeneFileName;
    %header = ();
    while (defined($str=<FIN>)){
        chomp $str;
        @ae = split /\t/, $str;
        if (scalar(keys %header)==0) {
            for (0..scalar(@ae)-1) { $header{$ae[$_]} = $_; }
            if (!exists($header{"SeqID"})) { die "There is no header 'SeqID' in file $codingGeneFileName\n"; }
            if (!exists($header{"Start"})) { die "There is no header 'Start' in file $codingGeneFileName\n"; } 
            if (!exists($header{"End"})) { die "There is no header 'End' in file $codingGeneFileName\n"; } 
            if (!exists($header{"Strand"})) { die "There is no header 'Strand' in file $codingGeneFileName\n"; }
            if (!exists($header{"COG"})) { $flagOfCOG = 0; }
        } else {
            if ($ae[$header{"SeqID"}] eq $contigName) {
                if ($flagOfCOG==0) {
                    push @dataCDS, [$ae[$header{"Start"}], $ae[$header{"End"}], $ae[$header{"Strand"}]];
                } else {
                    if ($ae[$header{"COG"}] eq "#") {
                        push @dataCDS, [$ae[$header{"Start"}], $ae[$header{"End"}], $ae[$header{"Strand"}]];
                    } else {
                        @subae = split /\s\|\s/, $ae[$header{"COG"}];
                        #if (index($subae[1],"[")!=-1) {
                        #    $subae[1] = substr $subae[1], 0, index($subae[1],"[")-1;
                        #}
                        push @dataCDS, [$ae[$header{"Start"}], $ae[$header{"End"}], $ae[$header{"Strand"}], @subae];
                    }  
                }
            }
        }
    }
    close FIN;
}
if (!(-e $inputHeader.".".$noncodingGeneFileName)) { $flagOfNoncodingGene = 0; } else {
    open FIN, $inputHeader.".".$noncodingGeneFileName;
    %header = ();
    while (defined($str=<FIN>)){
        chomp $str;
        @ae = split /\t/, $str;		
		if (scalar(keys %header)==0) {
			for (0..scalar(@ae)-1) { $header{$ae[$_]} = $_; }
			if (!exists($header{"SeqID"})) { die "There is no header 'SeqID' in file $noncodingGeneFileName\n"; }
			if (!exists($header{"Start"})) { die "There is no header 'Start' in file $noncodingGeneFileName\n"; } 
			if (!exists($header{"End"})) { die "There is no header 'End' in file $noncodingGeneFileName\n"; } 
			if (!exists($header{"Strand"})) { die "There is no header 'Strand' in file $noncodingGeneFileName\n"; }
			if (!exists($header{"Type"})) { die "There is no header 'Type' in file $noncodingGeneFileName\n"; }
			if (!exists($header{"Annotation"})) { die "There is no header 'Annotation' in file $noncodingGeneFileName\n"; }
		} else {
			if ($ae[$header{"SeqID"}] eq $contigName) {
				push @dataNoncoding, [$ae[$header{"Start"}], $ae[$header{"End"}], $ae[$header{"Strand"}], $ae[$header{"Type"}], $ae[$header{"Annotation"}]];
			}
		}
    }
    close FIN;
}

open FOUT, ">".$outputXML || die "Can not open output file $outputXML!\n";
print FOUT '<?xml version="1.0" encoding="ISO-8859-1" ?>', "\n";
print FOUT '<cgview'; foreach (@CGViewMain) { print FOUT " ", $_; } print FOUT ' sequenceLength="', $length, '">', "\n";
if (scalar(@textUperCenter)) {
    print FOUT '<legend position="upper-center" backgroundOpacity="0.8">', "\n";
    foreach (@textUperCenter) { print FOUT "\t", '<legendItem textAlignment="center" text="', $_, '" />', "\n"; }
    print FOUT "</legend>\n";
}
if (scalar(@textLowerCenter)) {
    print FOUT '<legend position="lower-center" backgroundOpacity="0.8">', "\n";
    foreach (@textLowerCenter) { print FOUT "\t", '<legendItem textAlignment="center" text="', $_, '" />', "\n"; }
    print FOUT "</legend>\n";
}
my $meanOfGCskew = 0;
my $SDOfGCskew = 0;
my $meanOfGCContent = 0;
my $SDOfGCContent = 0;
if ($flagOfGC) {    
    print FOUT '<legend position="lower-right" backgroundOpacity="0.8">', "\n";
    foreach (@dataGCskew) {
        $meanOfGCskew += $_->[2];
        $SDOfGCskew += $_->[2]**2;
    }
    $meanOfGCskew = $meanOfGCskew/scalar(@dataGCskew);
    $SDOfGCskew = sqrt($SDOfGCskew/scalar(@dataGCskew)-$meanOfGCskew**2);
    foreach (@dataGCContent) {
        $meanOfGCContent += $_->[2];
        $SDOfGCContent += $_->[2]**2;
    }
    $meanOfGCContent = $meanOfGCContent/scalar(@dataGCskew);
    $SDOfGCContent = sqrt($SDOfGCContent/scalar(@dataGCskew)-$meanOfGCContent**2);
    $meanOfGCskew = sprintf "%.2f", $meanOfGCskew;
    $SDOfGCskew = sprintf "%.2f", $SDOfGCskew;
    $meanOfGCContent = sprintf "%.2f", $meanOfGCContent;
    $SDOfGCContent = sprintf "%.2f", $SDOfGCContent;
    print FOUT "\t", '<legendItem textAlignment="left" text="GC-skew = ', $meanOfGCskew."(".$SDOfGCskew.")" , '" />', "\n";
    print FOUT "\t", '<legendItem textAlignment="left" text="GC-content = ', $meanOfGCContent."(".$SDOfGCContent.")" , '" />', "\n";
    print FOUT "</legend>\n";
}
print FOUT '<legend position="', $legendPosition, '" textAlignment="left" backgroundOpacity="0.8">', "\n";
if ($flagOfCDS && !$flagOfCOG) {
    if (exists ($legendColor{"CDS"})) { print FOUT "\t", '<legendItem text="CDS" drawSwatch="true" swatchColor="', $legendColor{"CDS"}, '" />', "\n"; }
    else { print FOUT '<legendItem text="CDS" drawSwatch="true" swatchColor="', $defaultcolor ,'" />', "\n"; }
}
if ($flagOfNoncodingGene) {
    if (exists ($legendColor{"rRNA"})) { print FOUT "\t", '<legendItem text="rRNA" drawSwatch="true" swatchColor="', $legendColor{"rRNA"}, '" />', "\n"; }
    else { print FOUT '<legendItem text="rRNA" drawSwatch="true" swatchColor="', $defaultcolor ,'" />', "\n"; }
    if (exists ($legendColor{"tRNA"})) { print FOUT "\t", '<legendItem text="tRNA" drawSwatch="true" swatchColor="', $legendColor{"tRNA"}, '" />', "\n"; }
    else { print FOUT '<legendItem text="tRNA" drawSwatch="true" swatchColor="', $defaultcolor ,'" />', "\n"; }
    if (exists ($legendColor{"sRNA"})) { print FOUT "\t", '<legendItem text="sRNA" drawSwatch="true" swatchColor="', $legendColor{"sRNA"}, '" />', "\n"; }
    else { print FOUT '<legendItem text="sRNA" drawSwatch="true" swatchColor="', $defaultcolor ,'" />', "\n"; }
    if (exists ($legendColor{"CRISPR"})) { print FOUT "\t", '<legendItem text="CRISPR" drawSwatch="true" swatchColor="', $legendColor{"CRISPR"}, '" />', "\n"; }
    else { print FOUT '<legendItem text="CRISPR" drawSwatch="true" swatchColor="', $defaultcolor ,'" />', "\n"; }
}
if ($flagOfGC) {
    if (exists ($legendColor{"GCskew+"})) { print FOUT "\t", '<legendItem text="GCskew+" drawSwatch="true" swatchColor="', $legendColor{"GCskew+"}, '" />', "\n"; }
    else { print FOUT '<legendItem text="GCskew+" drawSwatch="true" swatchColor="', $defaultcolor ,'" />', "\n"; }
    if (exists ($legendColor{"GCskew-"})) { print FOUT "\t", '<legendItem text="GCskew-" drawSwatch="true" swatchColor="', $legendColor{"GCskew-"}, '" />', "\n"; }
    else { print FOUT '<legendItem text="GCskew-" drawSwatch="true" swatchColor="', $defaultcolor ,'" />', "\n"; }
    if (exists ($legendColor{"GCcontent"})) { print FOUT "\t", '<legendItem text="GC%" drawSwatch="true" swatchColor="', $legendColor{"GCcontent"}, '" />', "\n"; }
    else { print FOUT '<legendItem text="GC%" drawSwatch="true" swatchColor="', $defaultcolor ,'" />', "\n"; }
}
if ($flagOfCDS && $flagOfCOG) {
    foreach (@COGList) {
        if (exists ($COGcolor{$_})) { print FOUT "\t", '<legendItem text="COG ', $_, '" drawSwatch="true" swatchColor="', $COGcolor{$_}, '" />', "\n"; }
        else { print FOUT '<legendItem text="COG ', $_, '" drawSwatch="true" swatchColor="', $defaultcolor ,'" />', "\n"; }
    }
}
print FOUT "</legend>\n";

my $mouseover;
# direct two circle
if ($flagOfCDS) {
    if ($flagOfCOG) {
        print FOUT '<featureSlot showShading="false" strand="direct">', "\n";
        foreach (reverse @dataCDS) {
            if ($_->[2] eq "+") {
                if (scalar(@{$_})>3) {
			$_->[5]=~s/\s//g;
			my @ae=split /;/,$_->[6];
			$ae[0]=~s/^\s//g;
			print "$_->[6]\t$ae[0]\n";
                   $mouseover = $_->[5]."; ".$_->[6]."; ".$_->[4]."; ".$_->[0]." to ".$_->[1];	
                 	
        #            print FOUT "\t", '<feature color="', $COGcolor{substr($_->[6],0,1)}, '" decoration="clockwise-arrow" label="', $_->[5], '" mouseover="', $mouseover, '">', '<featureRange start="', $_->[0], '" stop="', $_->[1], '" /></feature>', "\n";
        		print FOUT"\t",'<feature color="',$COGcolor{$ae[0]},'" decoration="clockwise-arrow" label="',$_->[5],'" mouseover="',$mouseover, '">', '<featureRange start="', $_->[0], '" stop="', $_->[1], '" /></feature>', "\n";
                } else {
                    print FOUT "\t", '<feature color="', $legendColor{"CDS"}, '" decoration="clockwise-arrow" >', '<featureRange start="', $_->[0], '" stop="', $_->[1], '" /></feature>', "\n";
                }
            }
        }
        print FOUT '</featureSlot>', "\n";
    } else {
        foreach (reverse @dataCDS) {
            if ($_->[2] eq "+") {
                print FOUT "\t", '<feature color="', $legendColor{"CDS"}, '" decoration="clockwise-arrow" >', '<featureRange start="', $_->[0], '" stop="', $_->[1], '" /></feature>', "\n";
            }
        }
        print FOUT '</featureSlot>', "\n";
    }
    print FOUT '<featureSlot featureThickness="0.5" showShading="true" strand="direct"><feature color="rgb(102,102,102)" decoration="arc"><featureRange start="1" stop="', $length, '"></featureRange></feature></featureSlot>', "\n";
}
if ($flagOfNoncodingGene) {
    print FOUT '<featureSlot showShading="false" strand="direct">', "\n";
    foreach (reverse @dataNoncoding) {
        if ($_->[2] eq "+") {
            $mouseover = $_->[4]."; ".$_->[0]." to ".$_->[1];
            print FOUT "\t", '<feature color="', $legendColor{$_->[3]}, '" decoration="clockwise-arrow" label="', $_->[4], '" mouseover="', $mouseover, '">', '<featureRange start="', $_->[0], '" stop="', $_->[1], '" /></feature>', "\n";
        }
    }
    print FOUT '</featureSlot>', "\n";
    print FOUT '<featureSlot featureThickness="0.5" showShading="true" strand="direct"><feature color="rgb(102,102,102)" decoration="arc"><featureRange start="1" stop="', $length, '"></featureRange></feature></featureSlot>', "\n";
}
# reverse two circle
if ($flagOfCDS) {
    if ($flagOfCOG) {
        print FOUT '<featureSlot showShading="false" strand="reverse">', "\n";
        foreach (@dataCDS) {
            if ($_->[2] eq "-") {
                if (scalar(@{$_})>3) {
		    my @ae=split /;/,$_->[6];
                    $mouseover = $_->[5]."; ".$_->[6]."; ".$_->[4]."; ".$_->[0]." to ".$_->[1];
                    print FOUT "\t", '<feature color="', $COGcolor{$ae[0]}, '" decoration="counterclockwise-arrow" label="', $_->[5], '" mouseover="', $mouseover, '">', '<featureRange start="', $_->[0], '" stop="', $_->[1], '" /></feature>', "\n";
                } else {
                    print FOUT "\t", '<feature color="', $legendColor{"CDS"}, '" decoration="counterclockwise-arrow" >', '<featureRange start="', $_->[0], '" stop="', $_->[1], '" /></feature>', "\n";
                }
            }
        }
        print FOUT '</featureSlot>', "\n";
    } else {
        foreach (@dataCDS) {
            if ($_->[2] eq "-") {
                print FOUT "\t", '<feature color="', $legendColor{"CDS"}, '" decoration="counterclockwise-arrow" >', '<featureRange start="', $_->[0], '" stop="', $_->[1], '" /></feature>', "\n";
            }
        }
        print FOUT '</featureSlot>', "\n";
    }
    print FOUT '<featureSlot featureThickness="0.5" showShading="true" strand="reverse"><feature color="rgb(102,102,102)" decoration="arc"><featureRange start="1" stop="', $length, '"></featureRange></feature></featureSlot>', "\n";
}
if ($flagOfNoncodingGene) {
    print FOUT '<featureSlot showShading="false" strand="reverse">', "\n";
    foreach (@dataNoncoding) {
        if ($_->[2] eq "-") {
            $mouseover = $_->[4]."; ".$_->[0]." to ".$_->[1];
            print FOUT "\t", '<feature color="', $legendColor{$_->[3]}, '" decoration="counterclockwise-arrow" label="', $_->[4], '" mouseover="', $mouseover, '">', '<featureRange start="', $_->[0], '" stop="', $_->[1], '" /></feature>', "\n";
        }
    }
    print FOUT '</featureSlot>', "\n";
    print FOUT '<featureSlot featureThickness="0.5" showShading="true" strand="reverse"><feature color="rgb(102,102,102)" decoration="arc"><featureRange start="1" stop="', $length, '"></featureRange></feature></featureSlot>', "\n";
}

if ($flagOfGC) {
    my $proportion;
    print FOUT '<featureSlot showShading="false" strand="reverse" featureThickness="', $radius/5, '">', "\n";
    print FOUT "\t", '<feature decoration="arc">', "\n";
    foreach (@dataGCContent) {
        $proportion = sprintf "%.3f", ($_->[2]-$meanOfGCContent)/$SDOfGCContent/10;
        if ($proportion<-0.5) { $proportion = -0.5; }
        if ($proportion>0.5) { $proportion = 0.5; }
        if ($proportion<0) {
            print FOUT "\t\t", '<featureRange color="', $legendColor{"GCcontent"}, '" start="', $_->[0]+1, '" stop="', $_->[1], '" proportionOfThickness="', -$proportion, '" radiusAdjustment="', 0.5+$proportion/2, '" />', "\n";
        } else {
            print FOUT "\t\t", '<featureRange color="', $legendColor{"GCcontent"}, '" start="', $_->[0]+1, '" stop="', $_->[1], '" proportionOfThickness="', $proportion, '" radiusAdjustment="', 0.5+$proportion/2, '" />', "\n";
        }
    }
    print FOUT "</feature>\n</featureSlot>\n";
    print FOUT '<featureSlot featureThickness="0.5" showShading="true" strand="reverse"><feature color="rgb(102,102,102)" decoration="arc"><featureRange start="1" stop="', $length, '"></featureRange></feature></featureSlot>', "\n";

    print FOUT '<featureSlot showShading="false" strand="reverse" featureThickness="', $radius/5, '">', "\n";
    print FOUT "\t", '<feature decoration="arc">', "\n";
    foreach (@dataGCskew) {
        if ($SDOfGCskew==0) { $SDOfGCskew = 0.01; }
        $proportion = sprintf "%.3f", $_->[2]/$SDOfGCskew/10;
        if ($proportion<-0.5) { $proportion = -0.5; }
        if ($proportion>0.5) { $proportion = 0.5; }
        if ($proportion<0) {
            print FOUT "\t\t", '<featureRange color="', $legendColor{"GCskew-"}, '" start="', $_->[0]+1, '" stop="', $_->[1], '" proportionOfThickness="', -$proportion, '" radiusAdjustment="', 0.5+$proportion/2, '" />', "\n";
        } else {
            print FOUT "\t\t", '<featureRange color="', $legendColor{"GCskew+"}, '" start="', $_->[0]+1, '" stop="', $_->[1], '" proportionOfThickness="', $proportion, '" radiusAdjustment="', 0.5+$proportion/2, '" />', "\n";
        }
    }
    print FOUT "</feature>\n</featureSlot>\n";
}

print FOUT "</cgview>\n";
close FOUT;
