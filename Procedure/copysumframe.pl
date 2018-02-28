#!/usr/bin/perl
my $index = 232;
my $inputpdfR = "./Data/Shot".$index."/Clean_Data_Shot".$index."_Cam_18158/sumTrackR.tif";
my $inputpdfL = "./Data/Shot".$index."/Clean_Data_Shot".$index."_Cam_18333/sumTrackL.tif";
my $outputpdfR = "./sumTrackR".$index.".tif";
my $outputpdfL = "./sumTrackL".$index.".tif";                                                                                                                                   
system("cp $inputpdfR $outputpdfR");
system("cp $inputpdfL $outputpdfL");

