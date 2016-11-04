#!/usr/bin/perl
$ciffile=$ARGV[0];
$tmpfile="tmp";
@values=split(/\./,$ciffile);
$structfile=$values[0].".struct";
open(CIF,"$ciffile");
open(TMP,">$tmpfile");
$OCCUPFound=0;

#print "please enter symbol of Brave lattice\n";
#$lattice=<STDIN>;
#chop($lattice);

while($OCCUPFound==0)
{
 chomp($line=<CIF>);$_=$line;
 $endFound=0;
 if(/_symmetry_space_group_name_H-M/){@values=split(/\'/,$line);$values[1]=~s/M/m/go;$values[1]=~s/N/n/go;$symgp=$values[1];}
 if(/_symmetry_Int_Tables_number/){@values=split(/\ +/,$line);$symno=$values[1];}
 
 if(/_symmetry_equiv_pos_as_xyz/)
 {
   for($i=1;$i<=1000;$i++)
   {
     chomp($line=<CIF>);$_=$line;
     if(/_cell_length_a/){@values=split(/\ +/,$line);$numop=$i-1;$A=$values[1]/0.529177;last;}
     else{$op[$i]=$line;print "$op[$i]\n";}
   }
 } 
 if(/cell_length_b/) {@values=split(/\ +/,$line);$B=$values[1]/0.529177;}
 if(/cell_length_c/) {@values=split(/\ +/,$line);$C=$values[1]/0.529177;}
 if(/cell_angle_alpha/) {@values=split(/\ +/,$line);$alpha=$values[1];}
 if(/cell_angle_beta/) {@values=split(/\ +/,$line);$beta=$values[1];}
 if(/cell_angle_gamma/) {@values=split(/\ +/,$line);$gamma=$values[1];}
 if(/atom_site_occupancy/){$OCCUPFound=1;}
}
$numele=0;
$i=0;
while(!eof(CIF))
{
 chomp($line=<CIF>);@values=split(/\ +/,$line);
 if($values[0] eq "loop_"){last;}
 $i++;$x[$i]=$values[2]+1;$y[$i]=$values[3]+1;$z[$i]=$values[4]+1;

 if($numele == 0){$numele++;$ele[$numele]=$values[1];}
 if($numele >= 1)
 {
   $match=0;
   foreach $ele(@ele)
   {
    if($values[1] eq $ele){$match=1;}
   }
   if($match == 0){$numele++;$ele[$numele]=$values[1];}
 }
 for($j=1;$j<=$numele;$j++)
 {
   if($values[1] eq $ele[$j]){$ne[$i]=$j;}
 }
}
$numatom=$i;

close(CIF);

for($i=1;$i<=$numatom;$i++)
{
 if(($x[$i] <= 0.0)||($x[$i] >= 1.0)){$x[$i]=$x[$i]-int(eval($x[$i]));}
 if(($y[$i] <= 0.0)||($y[$i] >= 1.0)){$y[$i]=$y[$i]-int(eval($y[$i]));}
 if(($z[$i] <= 0.0)||($z[$i] >= 1.0)){$z[$i]=$z[$i]-int(eval($z[$i]));}
 $tmp=abs($y[$i]-0.08333);
 print "$tmp\n"; 
 if(abs($x[$i]-0.11111) <= 0.0001){$x[$i]=1/9;}
 if(abs($x[$i]-0.22222) <= 0.0001){$x[$i]=2/9;}
 if(abs($x[$i]-0.44444) <= 0.0001){$x[$i]=4/9;}
 if(abs($x[$i]-0.55556) <= 0.0001){$x[$i]=5/9;}
 if(abs($x[$i]-0.77778) <= 0.0001){$x[$i]=7/9;}
 if(abs($x[$i]-0.88889) <= 0.0001){$x[$i]=8/9;}
 if(abs($x[$i]-0.33333) <= 0.0001){$x[$i]=1/3;}
 if(abs($x[$i]-0.66667) <= 0.0001){$x[$i]=2/3;}
 if(abs($x[$i]-0.16667) <= 0.0001){$x[$i]=1/6;}
 if(abs($x[$i]-0.83333) <= 0.0001){$x[$i]=5/6;}
 if(abs($x[$i]-0.91667) <= 0.0001){$x[$i]=11/12;}
 if(abs($x[$i]-0.58333) <= 0.0001){$x[$i]=7/12;}
 if(abs($x[$i]-0.41667) <= 0.0001){$x[$i]=5/12;}
 if(abs($x[$i]-0.08333) <= 0.0001){$x[$i]=1/12;}
 if(abs($y[$i]-0.11111) <= 0.0001){$y[$i]=1/9;}
 if(abs($y[$i]-0.22222) <= 0.0001){$y[$i]=2/9;}
 if(abs($y[$i]-0.44444) <= 0.0001){$y[$i]=4/9;}
 if(abs($y[$i]-0.55556) <= 0.0001){$y[$i]=5/9;}
 if(abs($y[$i]-0.77778) <= 0.0001){$y[$i]=7/9;}
 if(abs($y[$i]-0.88889) <= 0.0001){$y[$i]=8/9;}
 if(abs($y[$i]-0.33333) <= 0.0001){$y[$i]=1/3;}
 if(abs($y[$i]-0.66667) <= 0.0001){$y[$i]=2/3;}
 if(abs($y[$i]-0.16667) <= 0.0001){$y[$i]=1/6;}
 if(abs($y[$i]-0.83333) <= 0.0001){$y[$i]=5/6;}
 if(abs($y[$i]-0.91667) <= 0.0001){$y[$i]=11/12;}
 if(abs($y[$i]-0.58333) <= 0.0001){$y[$i]=7/12;}
 if(abs($y[$i]-0.41667) <= 0.0001){$y[$i]=5/12;}
 if(abs($y[$i]-0.08333) <= 0.0001){$y[$i]=1/12;}
 if(abs($z[$i]-0.11111) <= 0.0001){$z[$i]=1/9;}
 if(abs($z[$i]-0.22222) <= 0.0001){$z[$i]=2/9;}
 if(abs($z[$i]-0.44444) <= 0.0001){$z[$i]=4/9;}
 if(abs($z[$i]-0.55556) <= 0.0001){$z[$i]=5/9;}
 if(abs($z[$i]-0.77778) <= 0.0001){$z[$i]=7/9;}
 if(abs($z[$i]-0.88889) <= 0.0001){$z[$i]=8/9;}
 if(abs($z[$i]-0.33333) <= 0.0001){$z[$i]=1/3;}
 if(abs($z[$i]-0.66667) <= 0.0001){$z[$i]=2/3;}
 if(abs($z[$i]-0.16667) <= 0.0001){$z[$i]=1/6;}
 if(abs($z[$i]-0.83333) <= 0.0001){$z[$i]=5/6;}
 if(abs($z[$i]-0.91667) <= 0.0001){$z[$i]=11/12;}
 if(abs($z[$i]-0.58333) <= 0.0001){$z[$i]=7/12;}
 if(abs($z[$i]-0.41667) <= 0.0001){$z[$i]=5/12;}
 if(abs($z[$i]-0.08333) <= 0.0001){$z[$i]=1/12;}
}
# symmetry operation

$l=0;
for($i=1;$i<=$numatom;$i++)
{
 
 # do symmetry operation
 for($j=1;$j<=$numop;$j++)
 {
  @ops=split(/\,/,$op[$j]);
  
  $ops[0]=~s/x/$x[$i]/go;$ops[0]=~s/y/$y[$i]/go;$ops[0]=~s/z/$z[$i]/go;$rx[$j]=eval($ops[0]);
  $ops[1]=~s/x/$x[$i]/go;$ops[1]=~s/y/$y[$i]/go;$ops[1]=~s/z/$z[$i]/go;$ry[$j]=eval($ops[1]);
  $ops[2]=~s/x/$x[$i]/go;$ops[2]=~s/y/$y[$i]/go;$ops[2]=~s/z/$z[$i]/go;$rz[$j]=eval($ops[2]);
 if($rx[$j] > 0.0){$rx[$j]=$rx[$j]-int(eval($rx[$j]));} else{$rx[$j]=$rx[$j]-int(eval($rx[$j]))+1;}  
 if($ry[$j] > 0.0){$ry[$j]=$ry[$j]-int(eval($ry[$j]));} else{$ry[$j]=$ry[$j]-int(eval($ry[$j]))+1;}  
 if($rz[$j] > 0.0){$rz[$j]=$rz[$j]-int(eval($rz[$j]));} else{$rz[$j]=$rz[$j]-int(eval($rz[$j]))+1;}  
# $rx[$j]=eval($rx[$j]);
# $ry[$j]=eval($ry[$j]);
# $rz[$j]=eval($rz[$j]);
 
 }

 # eliminate the same position & determine the element number
 for($j=1;$j<=$numop;$j++)
 {
  $double=0;
  for($o=-1;$o<=1;$o++){                       #eliminate atoms duplicated by lattice translation
  for($p=-1;$p<=1;$p++){
  for($q=-1;$q<=1;$q++){
  
  for($k=1;$k<=$j-1;$k++)
  {
    $dx=abs($rx[$j]-$rx[$k]+$o);$dy=abs($ry[$j]-$ry[$k]+$p);$dz=abs($rz[$j]-$rz[$k]+$q);
    if(($dx <= 0.001)&&($dy <= 0.001)){if(($dz<= 0.001)){$double=1;}}
  }
  
  }
  }
  }
  
  if($double == 0){$l++;$rrx[$l]=$rx[$j];$rry[$l]=$ry[$j];$rrz[$l]=$rz[$j];$rrne[$l]=$ne[$i];}
   
 }
 $sortmulti[$i]=$l;  
 print TMP "$rrx[13]\n";
}
  
$totatom=$l;









# Define RMT value for each element

for($i=1;$i<=$numele;$i++)
{
  print "Enter the RMT value for $ele[$i]\n";
  $RMT[$i]=<STDIN>;
  chop($RMT[$i]);
}

for($i=1;$i<=$numatom;$i++)
{
# print "$ele[$ne[$i]],$RMT[$ne[$i]],$x[$i],$y[$i],$z[$i]\n";
}

# output

open(STR,">$structfile");
@values=split(/\./,$structfile);print STR "$values[0]\n";
#printf (STR "%s   LATTICE,NONEQUIV.ATOMS:%3d%d_%s\n",$lattice,$numatom,$symno,$symgp);
printf (STR "P   LATTICE,NONEQUIV.ATOMS:%3d  1_P1\n",$numatom);
#printf (STR "%s   LATTICE,NONEQUIV.ATOMS:%3d 1 P1\n",$lattice,$numatom);
print STR "MODE OF CALC=RELA unit=bohr\n";
printf (STR "%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n",$A,$B,$C,$alpha,$beta,$gamma);

$numFe=0;
$numCo=0;
$numH=0;
$numC=0;
$numO=0;
$numMg=0;


for($i=1;$i<=$numatom;$i++)
{
 $multi=$sortmulti[$i]-$sortmulti[$i-1];
 $k=-1*$i;
 if($ele[$ne[$i]] eq "Mn"){$znum=25;$numFe++;$sortnum=$numFe;}
 if($ele[$ne[$i]] eq "Fe"){$znum=26;$numFe++;$sortnum=$numFe;}
 if($ele[$ne[$i]] eq "Co"){$znum=27;$numCo++;$sortnum=$numCo;}
 if($ele[$ne[$i]] eq "H"){$znum=1;$numH++;$sortnum=$numH;}
 if($ele[$ne[$i]] eq "He"){$znum=2;$numHe++;$sortnum=$numHe;}
 if($ele[$ne[$i]] eq "Li"){$znum=3;$numLi++;$sortnum=$numLi;}
 if($ele[$ne[$i]] eq "Ni"){$znum=28;$numNi++;$sortnum=$numNi;}
 if($ele[$ne[$i]] eq "Cu"){$znum=29;$numCu++;$sortnum=$numCu;}
 if($ele[$ne[$i]] eq "C"){$znum=6;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "Zr"){$znum=40;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "Ce"){$znum=58;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "Pr"){$znum=59;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "O"){$znum=8;$numO++;$sortnum=$numO;}
 if($ele[$ne[$i]] eq "F"){$znum=9;$numO++;$sortnum=$numF;}
 if($ele[$ne[$i]] eq "Ba"){$znum=56;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "Ca"){$znum=20;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "P"){$znum=15;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "As"){$znum=33;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "Mg"){$znum=12;$numMg++;$sortnum=$numMg;}
 if($ele[$ne[$i]] eq "Zn"){$znum=30;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "Sr"){$znum=38;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "B"){$znum=5;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "C"){$znum=6;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "N"){$znum=7;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "Y"){$znum=39;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "Al"){$znum=13;$numC++;$sortnum=$numC;}
 if($ele[$ne[$i]] eq "Rh"){$znum=45;$numC++;$sortnum=$numC;}

 
# if($x[$i] >=0.0){$x[$i]=$x[$i]-int($x[$i]);} else{$x[$i]=$x[$i]-int($x[$i]-1);}  
# if($y[$i] >=0.0){$y[$i]=$y[$i]-int($y[$i]);} else{$y[$i]=$y[$i]-int($y[$i]-1);}  
# if($z[$i] >=0.0){$z[$i]=$z[$i]-int($z[$i]);} else{$z[$i]=$z[$i]-int($z[$i]-1);}  
 printf (STR "ATOM%4d: X=%10.8f Y=%10.8f Z=%10.8f\n",$k,$x[$i],$y[$i],$z[$i]);
 printf (STR "          MULT=%2d          ISPLIT= 8\n",$multi);
 for($l=$sortmulti[$i-1]+2;$l<=$sortmulti[$i];$l++)
 {
 if($rrx[$l] >=0.0){$rrx[$l]=$rrx[$l]-int($rrx[$l]);} else{$rrx[$l]=$rrx[$l]-int($rrx[$l]-1);}  
 if($rry[$l] >=0.0){$rry[$l]=$rry[$l]-int($rry[$l]);} else{$rry[$l]=$rry[$l]-int($rry[$l]-1);}  
 if($rrz[$l] >=0.0){$rrz[$l]=$rrz[$l]-int($rrz[$l]);} else{$rrz[$l]=$rrz[$l]-int($rrz[$l]-1);}  

  printf (STR "ATOM%4d: X=%10.8f Y=%10.8f Z=%10.8f\n",$k,$rrx[$l],$rry[$l],$rrz[$l]);
 }
# printf (STR "%-2s%-2d       NPT=  781  R0=0.00010000 RMT=%10.4f   Z:%5.1f\n",$ele[$ne[$i]],$sortnum,$RMT[$ne[$i]],$znum);
if($ele[$ne[$i]] eq "H")
 {printf (STR "%-2s         NPT=  781  R0=0.00010000 RMT=%10.4f   Z:%5.1f\n",$ele[$ne[$i]],$RMT[$ne[$i]],$znum);}
else
 {printf (STR "%-2s         NPT=  781  R0=0.00005000 RMT=%10.4f   Z:%5.1f\n",$ele[$ne[$i]],$RMT[$ne[$i]],$znum);}
 print STR "LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000\n";
 print STR "                     0.0000000 1.0000000 0.0000000\n";
 print STR "                     0.0000000 0.0000000 1.0000000\n";
}
print STR "   0      NUMBER OF SYMMETRY OPERATIONS\n";
close(STR);



