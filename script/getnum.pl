#!/usr/bin/perl -w
use strict;

my $file=$ARGV[0];

open FA,"$file" or die $!;
my $length;
while(<FA>){
	chomp;
	next if(/^>/);
	$length+=length($_);
}
close FA;

my $num=sprintf "%.4f",($length*2/1000000);
print "$num\n";
