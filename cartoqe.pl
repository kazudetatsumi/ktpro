#!/usr/bin/perl
# Description: This script converts atomic coordinations of a crystal stored in a .car 
# file to those in a Quantum Espresso input .in.
# Kazuyoshi TATSUMI 2019 12

$infile="$ARGV[0]";
@values=split(/\./,$infile);
$outfile=$values[0].".pw.in";
open(IN,"$infile")||die "cannot open $infile\n";
$PBCFound=0;
while($PBCFound!=1)
{
  chomp($line=<IN>);@values=split(/\ +/,$line);if($values[0] eq "PBC"){$PBCFound=1;}
}
@values=split(/\ +/,$line);$A=$values[1];$B=$values[2];$C=$values[3];
$alpha=$values[4]/180*3.141592653589790;$beta=$values[5]/180*3.141592653589790;
$gamma=$values[6]/180*3.141592653589790;
#Converting A,B,C,alpha,beta,gamma to a11,a12,a13.....

#print "Rhombohedral? [y/n] \n";
#$yesno=<STDIN>;
#if($yesno eq "y\n"){
#$at=sqrt(2/3*(1-cos($alpha))*$A**2);
#$ct=sqrt(1/3*$A**2+2/3*cos($alpha)*$A**2);
#$a1=$at;$a2=0;$a3=$ct;
#$b1=-1/2*$at;$b2=sqrt(3)/2*$at;$b3=$ct;
#$c1=-1/2*$at;$c2=-sqrt(3)/2*$at;$c3=$ct;
#print "cartitian is set for rhombohedral lattice\n";
#}
#A axis // X,  B axis on XY
#print "A axis along x and B axis on xy? [y/n] \n";
#$yesno=<STDIN>;
#print $yesno;
#if($yesno eq "y\n") {
$a1=$A;$a2=0;$a3=0;
$b1=$B*cos($gamma);$b2=$B*sin($gamma);$b3=0;
$c1=$C*cos($beta);$c2=($C*cos($alpha)-$c1*cos($gamma))/sin($gamma);
$c3=sqrt($C**2-$c1**2-$c2**2);
print "cartitian is set as A axis //x, B axis on xy\n";
#}

$det=$a1*($b2*$c3-$b3*$c2)-$a2*($b1*$c3-$b3*$c1)+$a3*($b1*$c2-$b2*$c1);
$b11=($b2*$c3-$b3*$c2)/$det;$b12=-($b1*$c3-$b3*$c1)/$det;$b13=($b1*$c2-$b2*$c1)/$det;
$b21=-($a2*$c3-$a3*$c2)/$det;$b22=($a1*$c3-$a3*$c1)/$det;$b23=-($a1*$c2-$a2*$c1)/$det;
$b31=($a2*$b3-$a3*$b2)/$det;$b32=-($a1*$b3-$a3*$b1)/$det;$b33=($a1*$b2-$a2*$b1)/$det;

$tst11=$b11*$a11+$b12*$a21+$b13*$a31;
$tst12=$b11*$a12+$b12*$a22+$b13*$a32;
$tst13=$b11*$a13+$b12*$a23+$b13*$a33;
print "volume: $det \n";
print "please check lattice vectors\n";
print " $a1  $b1  $c1  \n";
print " $a2  $b2  $c2  \n";
print " $a3  $b3  $c3  \n";
#print " $b11  $b12  $b13  \n";
#print " $b21  $b22  $b23  \n";
#print " $b31  $b32  $b33  \n";

# For atoms
print "please enter number of elements \n";
$NumEle=<STDIN>;
for($i=1;$i<=$NumEle;$i++)
 {
   print "please enter the $i th element name\n";
   $Elename[$i]=<STDIN>;
   $ELename[$i]=chop($Elename[$i]);
   print "$Elename[$i]\n";
   $Num[$i]=0;
   if ($Elename[$i] eq "Na"){$AM[$i]=22.98976928;$pseu[$i]="Na.pbe-spn-kjpaw_psl.0.2.UPF";}
   elsif ($Elename[$i] eq "Cl") {$AM[$i]=35.453;$pseu[$i]="Cl.pbe-n-kjpaw_psl.0.1.UPF";}
 }
print "please enter the TOTAL number of atoms \n";
$NumAtom=<STDIN>;
close(IN);

$k=0;
for($i=1;$i<=$NumEle;$i++)
 {
   open(IN,"$infile");
   $PBCFound=0;
   while($PBCFound!=1)
    {
       chomp($line=<IN>);@values=split(/ +/,$line);
       if($values[0] eq "PBC"){$PBCFound=1;print "$PBCFound\n";}
    } 
   $Comp[$i]=0;
   for($j=1;$j<=$NumAtom;$j++)
    {
      chomp($line=<IN>);@values=split(/ +/,$line);
#      print "$values[1] $values[2] $values[7]\n";
      if($values[7] eq "$Elename[$i]"){$k++;$x[$k]=$values[1];$y[$k]=$values[2];$z[$k]=$values[3];$Ele[$k]=$Elename[$i];$Comp[$i]++;
                                      print "$values[7]\n";
				        }
    }
   
   close(IN);
 }
print "$k\n";

for($i=1;$i<=$NumAtom;$i++)
 {
   $u[$i]=$x[$i]*$b11+$y[$i]*$b12+$z[$i]*$b13;
   $v[$i]=$x[$i]*$b21+$y[$i]*$b22+$z[$i]*$b23;
   $w[$i]=$x[$i]*$b31+$y[$i]*$b32+$z[$i]*$b33;
   
 }
   
$particles=$Comp[1];
for($i=2;$i<=$NumEle;$i++)
 {
   $particles=$particles." ".$Comp[$i];
 }
 
print "$particles \n";

#sorting
$preb=0;
for($i=1;$i<=$NumEle;$i++)
{
  if($i>=2){$preb=$Comp[$i-1]+$preb;}
  $exchange=1;
  while($exchange!=0){
  $exchange=0;
  for($j=$preb+1;$j<=$preb+$Comp[$i]-1;$j++)
  {
    if($w[$j+1] > $w[$j]){$tmpu=$u[$j+1];$tmpv=$v[$j+1];$tmpw=$w[$j+1];$u[$j+1]=$u[$j];$v[$j+1]=$v[$j];$w[$j+1]=$w[$j];
			  $u[$j]=$tmpu;$v[$j]=$tmpv;$w[$j]=$tmpw;$exchange=1;}
			 }
   }
}

#OUTPUT
open (OUT, ">$outfile");
  print OUT "&control\n";
  print OUT "    calculation = 'scf'\n";
  print OUT "    tprnfor = .true.\n";
  print OUT "    tstress = .true.\n";
  print OUT "    pseudo_dir = '/home/kazu/code/q-e-qe-6.4.1/pseudo/'\n";
  print OUT "/\n";
  print OUT "&system\n";
  print OUT "    ibrav = 0\n";
  print OUT "    nat = 8\n";
  printf (OUT "    ntyp = %d \n", $NumEle);
  print OUT "    ecutwfc = 70.0\n";
  print OUT "/\n";
  print OUT "&electrons\n";
  print OUT "    diagonalization = 'david'\n";
  print OUT "    conv_thr = 1.0d-9\n";
  print OUT "/\n";   
  print OUT "ATOMIC_SPECIES\n";
  for($i=1;$i<=$NumEle;$i++)
   { 
     printf (OUT " %s %12.8f  %s \n", $Elename[$i],$AM[$i],$pseu[$i]);
  }
  print OUT "ATOMIC_POSITIONS crystal\n";
  for($i=1;$i<=$NumAtom;$i++)
   { 
     printf (OUT " %s %22.16f%22.16f%22.16f  \n", $Ele[$i], $u[$i],$v[$i],$w[$i]);
  }
  print OUT "CELL_PARAMETERS angstrom \n";
  printf (OUT "%22.16f%22.16f%22.16f \n",$a1,$a2,$a3);
  printf (OUT "%22.16f%22.16f%22.16f \n",$b1,$b2,$b3);
  printf (OUT "%22.16f%22.16f%22.16f \n",$c1,$c2,$c3);
  print OUT "K_POINTS automatic\n";
  print OUT " 8 8 8 1 1 1\n";
close(OUT);
