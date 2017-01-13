#!/usr/bin/perl
# .cif file is converted as a POSCAR format file by this script
# the crystal structure in the .cif file must be discribed with P1 symmetry.

$infile=$ARGV[0];
@values=split(/\./,$infile);
$posfile=$values[0].".poscar";
#Reading the cif file
open(IN,"$infile") || die "Cannot open $infile ";
$OccupFound=0;
while($OccupFound!=1)
 {
   chomp($line=<IN>);@values=split(/ +/,$line);
   if($values[0] eq "_atom_site_occupancy"){$OccupFound=1;}
   if($values[0] eq "_cell_length_a"){$A=$values[1];}
   if($values[0] eq "_cell_length_b"){$B=$values[1];}
   if($values[0] eq "_cell_length_c"){$C=$values[1];}
   if($values[0] eq "_cell_angle_alpha"){$alpha=$values[1]/180*3.141592653589790;}
   if($values[0] eq "_cell_angle_beta"){$beta=$values[1]/180*3.141592653589790;}
   if($values[0] eq "_cell_angle_gamma"){$gamma=$values[1]/180*3.141592653589790;}
 }
#Converting A,B,C,alpha,beta,gamma to ax,ay,az.....

$cx=0;$cy=0;$cz=$C;
$bx=0; $by=$B*sin($alpha); $bz=$B*cos($alpha);
$az=$A*cos($beta);
$ay=($A*cos($gamma)-$az*cos($alpha))/sin($alpha);
$ax=sqrt($A**2-$ay**2-$az**2); 

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
 }
print "please enter the TOTAL number of atoms \n";
$NumAtom=<STDIN>;
close(IN);

$k=0;
for($i=1;$i<=$NumEle;$i++)
 {
   open(IN,"$infile");
   $OccupFound=0;
   while($OccupFound!=1)
    {
       chomp($line=<IN>);@values=split(/ +/,$line);
       if($values[0] eq "_atom_site_occupancy"){$OccupFound=1;print "$OccupFound\n";}
    } 
   $Comp[$i]=0;
   for($j=1;$j<=$NumAtom;$j++)
    {
      chomp($line=<IN>);@values=split(/ +/,$line);
      print "$values[1]\n";
      if($values[1] eq "$Elename[$i]"){$k++;$x[$k]=$values[2];$y[$k]=$values[3];$z[$k]=$values[4];$Ele[$k]=$Elename[$i];$Comp[$i]++;}
    }
   
   close(IN);
 }
print "$k\n";
if($NumAtom != $k){die "Number of atoms is unmatched !";}

for($i=1;$i<=$NumAtom;$i++)
 {
   print "$Ele[$i],$x[$i],$y[$i],$z[$i]\n";
 }
   
$particles=$Comp[1];
for($i=2;$i<=$NumEle;$i++)
 {
   $particles=$particles." ".$Comp[$i];
 }
 
print "$particles \n";


#OUTPUT
open (OUT, ">$posfile");
  print OUT "$infile\n";
  print OUT "   1.0000000\n";
  printf (OUT "%22.16f%22.16f%22.16f \n",$ax,$ay,$az);
  printf (OUT "%22.16f%22.16f%22.16f \n",$bx,$by,$bz);
  printf (OUT "%22.16f%22.16f%22.16f \n",$cx,$cy,$cz);
  printf (OUT " %s \n",$particles);
  print OUT "Selective Dynamics\n";
  print OUT "Direct\n";
  for($i=1;$i<=$NumAtom;$i++)
   { $t=abs($z[$i]);
     if($t<=0.1700){
     printf (OUT "%11.5f%11.5f%11.5f  F  F  F  %s \n", $x[$i],$y[$i],$z[$i],$Ele[$i])
     }
     elsif(abs($t-0.500)<=0.0001){
     printf (OUT "%11.5f%11.5f%11.5f  F  F  F  %s \n", $x[$i],$y[$i],$z[$i],$Ele[$i])
     }
     else{
     printf (OUT "%11.5f%11.5f%11.5f  T  T  T  %s \n", $x[$i],$y[$i],$z[$i],$Ele[$i]);
     }
  }
close(OUT);
   
