#!/usr/bin/perl
die "perl $0 *fq > *fa\n" if(@ARGV != 1);
my $fq=shift;
open(IN,$fq) || die;
while(my $line=<IN>){
	my $seq=<IN>; <IN>; <IN>;
	print ">$line";
	print $seq;
}
	
