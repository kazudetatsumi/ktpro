#!/usr/bin/perl
#submit.pl 
# this script read states of Sun grid engine working queues and 
# write .machines of WIEN2k file for parallel job. 

$listfile="tmplist";
$machinefile=".machines";
system "rm .machines";
system "rm tmplist";
system "qstat > $listfile";
open(LS,"$listfile");
chomp($line=<LS>);
chomp($line=<LS>);
$i=0;
while(!eof(LS)){
  chomp($line=<LS>);@values=split(/\ +/,$line);
  $_=$values[-2];
  if(/mach/){
  @dateset=split(/\//,$values[-4]);
  @timeset=split(/\:/,$values[-3]);
  $i++;
  $mm[$i]=$dateset[0];$dd[$i]=$dateset[1];$yy[$i]=$dateset[2];
  $t[$i]=$timeset[0];$m[$i]=$timeset[1];$s[$i]=$timeset[2];
  $totalsecond[$i]=$s[$i]+60*$m[$i]+60*60*$t[$i]+60*60*24*($dd[$i]+$mm[$i]*31+$yy[$i]*365);
  print "$dateset[0] $dateset[1] $dateset[2] $dateset[3]\n";
  print "$timeset[0] $timeset[1] $timeset[2] $timeset[3] $totalsecond[$i]\n";
  @values2=split(/\@/,$values[-2]);
  $queue[$i]=$values2[-1];$master[$i]=$values[-1];
  $ncore[$i]=$values[$#arr];
  print "ncore, $ncore[$i]\n";
  }
}
$j=0;
$numqueue=$i;
$lid=1;
for ($i=2;$i<=$numqueue;$i++){
     if($totalsecond[$i] >= $totalsecond[$lid]){$lid=$i;}
}
for($i=1;$i<=$numqueue;$i++){
  if($totalsecond[$i] == $totalsecond[$lid]){
    $wid=$i;
  }
}
open(MAC,">$machinefile");
for ($i=1;$i<=$ncore[$wid];$i++){
  print MAC "1:$queue[$wid] \n";
}     
