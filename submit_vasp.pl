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
  $queue[$i]=$values[-2];$master[$i]=$values[-1];
  chop($queue[$i]);chop($queue[$i]);
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
    $j++;
    $wid[$j]=$i;
  }
}
$numw=$j;
#Grouping for dual core!!
for($i=1;$i<=$numw;$i++){
$test[$i]=0;
$deg[$i]=1;
}
$ne=0;
for($i=1;$i<=$numw;$i++){
    
    print "$test[$i]\n";
    if ($test[$i] == 0){
    $ne++;
    $qname[$ne]=$queue[$wid[$i]];
    for($j=$i+1;$j<=$numw;$j++){
       if($queue[$wid[$j]] eq $queue[$wid[$i]]){
          $test[$j] = 1;
	  $deg[$ne]++;}
    }
  }
}
#######

$nne=$ne;  
$lapw0="lapw0:";
for ($ne=1;$ne<=$nne;$ne++){
  $lapw0="$lapw0"."$qname[$ne]:$deg[$ne] ";
}
print "$lapw0 \n";
open(MAC,">$machinefile");
#  print MAC "$lapw0\n";
for ($i=1;$i<=$numw;$i++){
#  print MAC "1:$queue[$wid[$i]] \n";
  print MAC "$queue[$wid[$i]] \n";
}     
