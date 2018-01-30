#!/usr/bin/perl
my $inputfile = "PairList.csv";

open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(',',$line);
    my $pidl = $array[1];
    my $pidr = $array[2];
    my $file1 = "./Compare/Compare_".$pidl."_".$pidr.".tif";
    system("cp $file1 .");
    print $pidl,",",$pidr,"\n";
}

