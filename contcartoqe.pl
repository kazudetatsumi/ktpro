#!/usr/bin/perl
# usage: contcartoqe.pl outputfilename [-vc] 
$incar="INCAR";
$posfile="CONTCAR";
$elefile="ele";

# getting element name
system "grep VRHF POTCAR > $elefile";
open(ELE,"$elefile");
while(!eof(ELE))
{
  chomp($line=<ELE>);@values=split(/\=/,$line);@values2=split(/\:/,$values[1]);
  $j++;
  $values2[0]=~ s/ //g;  
  $ele[$j]=$values2[0];
  print $ele[$j];
  if ($ele[$j] eq "Na"){$AM[$j]=22.98976928;$pseu[$j]="Na.pbe-spn-kjpaw_psl.0.2.UPF";}
  elsif ($ele[$j] eq "Cl") {$AM[$j]=35.453;$pseu[$j]="Cl.pbe-n-kjpaw_psl.0.1.UPF";}  
  elsif ($ele[$j] eq "Ni") {$AM[$j]=58.6934;$pseu[$j]="Ni.pbe-n-kjpaw_psl.0.1.UPF";}  
  elsif ($ele[$j] eq "La") {$AM[$j]=138.90547;$pseu[$j]="La.GGA-PBE-paw-v1.0.UPF";}
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
print $ARGV[1];
#output
if($ARGV[1] eq "-vc"){
$outfile="pw.vc-relax.".$ARGV[0].".in";
open (OUT, ">$outfile");
  print OUT "&control\n";
  print OUT "    calculation = 'vc-relax'\n";
  print OUT "    tprnfor = .true.\n";
  print OUT "    tstress = .true.\n";
  print OUT "    pseudo_dir = '/home/kazu/code/q-e-qe-6.4.1/pseudo/'\n";
  print OUT "    etot_conv_thr = 1.0d-9\n";
  print OUT "/\n";
  print OUT "&system\n";
  print OUT "    ibrav = 0\n";
  printf (OUT "    nat = %d\n", $numatom);
  printf (OUT "    ntyp = %d \n", $numele);
  print OUT "    ecutwfc = 70.0\n";
  print OUT "/\n";
  print OUT "&electrons\n";
  print OUT "    diagonalization = 'david'\n";
  print OUT "    conv_thr = 1.0d-9\n";
  print OUT "/\n";
  print OUT "&ions\n";
  print OUT "/\n";
  print OUT "&cell\n";
  print OUT "/\n";
  print OUT "ATOMIC_SPECIES\n";
  for($i=1;$i<=$numele;$i++)
   { 
     printf (OUT " %s %12.8f  %s \n", $ele[$i],$AM[$i],$pseu[$i]);
  }
  print OUT "ATOMIC_POSITIONS crystal\n";
  $k=0;
  for($j=1;$j<=$numele;$j++)
  {
     for($i=1;$i<=$nume[$j];$i++)
     { 
       $k++;
       printf (OUT "%s %22.16f%22.16f%22.16f  \n", $ele[$j],$x[$k],$y[$k],$z[$k]);
     }
  }
  print OUT "CELL_PARAMETERS angstrom \n";
  printf (OUT "%22.16f%22.16f%22.16f \n",$ax,$ay,$az);
  printf (OUT "%22.16f%22.16f%22.16f \n",$bx,$by,$bz);
  printf (OUT "%22.16f%22.16f%22.16f \n",$cx,$cy,$cz);
  print OUT "K_POINTS automatic\n";
  print OUT " 8 8 8 1 1 1\n";
}
else{
$outfile="pw.scf.".$ARGV[0].".in";
open (OUT, ">$outfile");
  print OUT "&control\n";
  print OUT "    calculation = 'scf'\n";
  print OUT "    tprnfor = .true.\n";
  print OUT "    tstress = .true.\n";
  print OUT "    pseudo_dir = '/home/kazu/code/q-e-qe-6.4.1/pseudo/'\n";
  print OUT "/\n";
  print OUT "&system\n";
  print OUT "    ibrav = 0\n";
  printf (OUT "    nat = %d\n", $numatom);
  printf (OUT "    ntyp = %d\n", $numele);
  print OUT "    ecutwfc = 70.0\n";
  print OUT "/\n";
  print OUT "&electrons\n";
  print OUT "    diagonalization = 'david'\n";
  print OUT "    conv_thr = 1.0d-9\n";
  print OUT "/\n";   
  print OUT "ATOMIC_SPECIES\n";
  for($i=1;$i<=$numele;$i++)
   { 
     printf (OUT " %s %12.8f  %s \n", $ele[$i],$AM[$i],$pseu[$i]);
  }
  print OUT "ATOMIC_POSITIONS crystal\n";
  $k=0;
  for($j=1;$j<=$numele;$j++)
  {
     for($i=1;$i<=$nume[$j];$i++)
     { 
       $k++;
       printf (OUT "%s %22.16f%22.16f%22.16f  \n", $ele[$j],$x[$k],$y[$k],$z[$k]);
     }
  }
  print OUT "CELL_PARAMETERS angstrom \n";
  printf (OUT "%22.16f%22.16f%22.16f \n",$ax,$ay,$az);
  printf (OUT "%22.16f%22.16f%22.16f \n",$bx,$by,$bz);
  printf (OUT "%22.16f%22.16f%22.16f \n",$cx,$cy,$cz);
  print OUT "K_POINTS automatic\n";
  print OUT " 8 8 8 1 1 1\n";
close(OUT);
}
