#!/usr/bin/perl
$out="OUTCAR";
$pos="POSCAR";
open(OUT,"$out")||die "Cannot open $out file";
open(POS,"$pos")||die "Cannot open $pos file";
$incar="INCAR";
$posfile="POSCAR";
$elefile="ele";
$ciffile=$ARGV[0].".cif";

# getting element name
system "grep VRHF POTCAR > $elefile";
open(ELE,"$elefile");
while(!eof(ELE))
{
  chomp($line=<ELE>);@values=split(/\=/,$line);@values2=split(/\:/,$values[1]);
  $j++;
  $ele[$j]=$values2[0];
}
$numele=$j;
close(ELE);

# getting the number of Oxygen and Hydrogen atoms
open(POS,"$posfile");
chomp($line=<POS>);
chomp($line=<POS>);@values=split(/\ +/,$line);
foreach  (@values){
	if(/[0-9]/)	{ $amplitude=$_;}
}
chomp($line=<POS>);@values=split(/\ +/,$line);
  $ax=$values[1]*$amplitude;$ay=$values[2]*$amplitude;$az=$values[3]*$amplitude;
  $A=sqrt($ax**2+$ay**2+$az**2);
#  print "$A\n";
chomp($line=<POS>);@values=split(/\ +/,$line);
  $bx=$values[1]*$amplitude;$by=$values[2]*$amplitude;$bz=$values[3]*$amplitude;
  $B=sqrt($bx**2+$by**2+$bz**2);
#  print "$B\n";
chomp($line=<POS>);@values=split(/\ +/,$line);
  $cx=$values[1]*$amplitude;$cy=$values[2]*$amplitude;$cz=$values[3]*$amplitude;
  $C=sqrt($cx**2+$cy**2+$cz**2);
#  print "$C\n";
  $tmp=($bx*$cx+$by*$cy+$bz*$cz)/($B*$C);
  $alpha=atan2((sqrt(1-$tmp**2)),$tmp)/3.141592654*180;
  $tmp=($ax*$cx+$ay*$cy+$az*$cz)/($A*$C);
  $beta=atan2((sqrt(1-$tmp**2)),$tmp)/3.141592654*180;
  $tmp=($ax*$bx+$ay*$by+$az*$bz)/($A*$B);
  $gamma=atan2((sqrt(1-$tmp**2)),$tmp)/3.141592654*180;
chomp($line=<POS>);
$_=$line;
if(/[a-zA-Z]/){chomp($line=<POS>);}
@values=split(/\ +/,$line);
print "$line\n";
  $numatom=0;
for($j=1;$j<=$numele;$j++){$nume[$j]=$values[$j];$numatom=$numatom+$nume[$j];}
print "numatom= $numatom \n";

chomp($line=<POS>);$_=$line;if(/elective/){chomp($line=<POS>);}
for($j=1;$j<=$numatom;$j++)
{
  chomp($line=<POS>);@values=split(/\ +/,$line);
  $x[$j]=$values[1];$y[$j]=$values[2];$z[$j]=$values[3];
  $_=$line;
  if(/F/){$fix[$j]=1;}
  else{$fix[$j]=0;}
}
close(POS);
$i=0;
while(!eof(OUT)){
  chomp($line=<OUT>);$_=$line;
  $fmax=0.00000;
  if(/sigma->0/){$ene=$line;}
  if(/TOTAL-FORCE/){
     $i++;
     chomp($line=<OUT>);
     for($j=1;$j<=$numatom;$j++){
        
	chomp($line=<OUT>);@values=split(/\ +/,$line);
	$fx[$j]=$values[4];
	$fy[$j]=$values[5];
	$fz[$j]=$values[6];
	$f[$j]=sqrt($values[4]**2+$values[5]**2+$values[6]**2);
	if($fix[$j] eq 0){
	if($f[$j] >= $fmax){$fmax=$f[$j];}
	}

    }
    printf ("%3d th FMAX= %7.4f   $ene \n",$i,$fmax);
    }
    }
    close(OUT);
    
	
