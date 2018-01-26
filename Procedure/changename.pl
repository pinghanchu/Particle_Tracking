#!/usr/bin/perl
my $inputfile = "test.txt";

open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(".pdf",$line);
   
    my $file1 = $array[0].".tif";
    my $file2 = $line;
    system("mv $file2 $file1");
    print $file1,",",$file2,"\n";
}

