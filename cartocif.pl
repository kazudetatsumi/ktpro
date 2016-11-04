#!/usr/bin/perl
$incar="INCAR";
$posfile="CONTCAR";
$elefile="ele";
$ciffile=$ARGV[0].".cif";

# getting element name
system "grep VRHF POTCAR > $elefile";
open(ELE,"$elefile");
while(!eof(ELE))
{
  chomp($line=<ELE>);@values=split(/\=/,$line);@values2=split(/\:/,$values[1]);
  $j++;
  $values2[0]=~ s/ //g;  
  $ele[$j]=$values2[0];
}
$numele=$j;
close(ELE);

# getting the number of Oxygen and Hydrogen atoms
open(POS,"$posfile");
chomp($line=<POS>);
chomp($line=<POS>);@values=split(/\ +/,$line);$amplitude=$values[1];
chomp($line=<POS>);@values=split(/\ +/,$line);
  $ax=$values[1]*$amplitude;$ay=$values[2]*$amplitude;$az=$values[3]*$amplitude;
  $A=sqrt($ax**2+$ay**2+$az**2);
  print "$A\n";
chomp($line=<POS>);@values=split(/\ +/,$line);
  $bx=$values[1]*$amplitude;$by=$values[2]*$amplitude;$bz=$values[3]*$amplitude;
  $B=sqrt($bx**2+$by**2+$bz**2);
  print "$B\n";
chomp($line=<POS>);@values=split(/\ +/,$line);
  $cx=$values[1]*$amplitude;$cy=$values[2]*$amplitude;$cz=$values[3]*$amplitude;
  $C=sqrt($cx**2+$cy**2+$cz**2);
  print "$C\n";
  $tmp=($bx*$cx+$by*$cy+$bz*$cz)/($B*$C);
  $alpha=atan2((sqrt(1-$tmp**2)),$tmp)/3.141592654*180;
  $tmp=($ax*$cx+$ay*$cy+$az*$cz)/($A*$C);
  $beta=atan2((sqrt(1-$tmp**2)),$tmp)/3.141592654*180;
  $tmp=($ax*$bx+$ay*$by+$az*$bz)/($A*$B);
  $gamma=atan2((sqrt(1-$tmp**2)),$tmp)/3.141592654*180;
chomp($line=<POS>);
chomp($line=<POS>);@values=split(/\ +/,$line);
  $numatom=0;
for($j=1;$j<=$numele;$j++){$nume[$j]=$values[$j];$numatom=$numatom+$nume[$j];}

chomp($line=<POS>);$_=$line;if(/elective/){chomp($line=<POS>);}
for($j=1;$j<=$numatom;$j++)
{
  chomp($line=<POS>);@values=split(/\ +/,$line);
  $x[$j]=$values[1];$y[$j]=$values[2];$z[$j]=$values[3];
}
close(POS);

open(CIF,">$ciffile");
print CIF "data_$ARGV[0]\n";
print CIF "_symmetry_space_group_name_H-M    'P1'\n";
print CIF "_symmetry_Int_Tables_number       1\n";
print CIF "_symmetry_cell_setting            triclinic\n";
print CIF "loop_\n";
print CIF "_symmetry_equiv_pos_as_xyz\n";
print CIF "  x,y,z \n";

printf (CIF "_cell_length_a                   %7.4f\n",$A);
printf (CIF "_cell_length_b                   %7.4f\n",$B);
printf (CIF "_cell_length_c                   %7.4f\n",$C);
printf (CIF "_cell_angle_alpha                %8.4f\n",$alpha);
printf (CIF "_cell_angle_beta                %8.4f\n",$beta);
printf (CIF "_cell_angle_gamma                %8.4f\n",$gamma);
print CIF "loop_\n";
print CIF "_atom_site_label\n";
print CIF "_atom_site_type_symbol\n";
print CIF "_atom_site_fract_x\n";
print CIF "_atom_site_fract_y\n";
print CIF "_atom_site_fract_z\n";
print CIF "_atom_site_U_iso_or_equiv\n";
print CIF "_atom_site_adp_type\n";
print CIF "_atom_site_occupancy\n";
$k=0;
for($j=1;$j<=$numele;$j++)
{
  for($i=1;$i<=$nume[$j];$i++)
  { 
    $k++;
    printf (CIF "%-7s%-3s%10.5f%10.5f%10.5f   0.00000  Uiso   1.00\n",$ele[$j].$i,$ele[$j],$x[$k],$y[$k],$z[$k]);
  }
}

close(CIF);
