#!/usr/bin/perl
$sum = 4745;
for($i = 0;$i<$sum;$i++){
    #$file1 = "./frame/frame_sum0_".$i.".tif";
    $file2 = "./frame/frame_".$i.".tif";
    $j = $sum-$i-1;
    $file3 = "./invframe/frame_".$j.".tif";
    #system("mv $file1 $file2");
    system("cp $file2 $file3");
}
