#!/usr/bin/perl
#This script read the qe output xml file and create a POSCAR file with the calculated atomic structure
#The generated input file is used to obtain further accurate atomic/electronic structure calculation.
#Kazuyoshi TATSUMI 2019 12 25
#latmat was transposed to work correctly.
use XML::TreePP;
use PDL;
use PDL::Matrix;
$atoangs = 0.529177210903;
$tpp = XML::TreePP->new();
$output_str = $tpp->parsefile( "./pwscf.save/data-file-schema.xml" )->{'qes:espresso'}->{output}->{atomic_structure};
$posfile="outfile.poscar";

$species = $tpp->parsefile( "./pwscf.save/data-file-schema.xml" )->{'qes:espresso'}->{output}->{atomic_species}{species};
for ($cnt = 0; $cnt < @$species; $cnt++) {
    $Elename[$cnt] = $species->[$cnt]{"-name"};
	$AM[$cnt] = $species->[$cnt]->{mass};
	$pseu[$cnt] = $species->[$cnt]->{pseudo_file};
	$Comp[$cnt] = 0;
}
$NumEle = $cnt;

$cell = $output_str->{cell};
@values1 = split(/\ +/, $cell->{a1}); 
@values2 = split(/\ +/, $cell->{a2}); 
@values3 = split(/\ +/, $cell->{a3}); 
$latmat = transpose mpdl [$values1[0], $values1[1], $values1[2]],[$values2[0], $values2[1], $values2[2]], [$values3[0], $values3[1], $values3[2]];
$inv_latmat =  matinv $latmat;
$latmat_angs = $latmat * $atoangs;
print $latmat_angs;

$atoms = $output_str->{atomic_positions}{atom};
for ($cnt = 0; $cnt < @$atoms; $cnt++) {
	$Ele[$cnt] = $atoms->[$cnt]{"-name"};
	@values = split(/\ +/, $atoms->[$cnt]{"#text"});
	$vecs[$cnt] = $inv_latmat x vpdl [$values[0], $values[1], $values[2]]; 
	for ($i=0;$i<=$NumEle;$i++)
	 {
	   if($Ele[$cnt] eq $Elename[$i]){$Comp[$i]++;}
	 }
}
$NumAtom = $cnt;
$particles=$Comp[0];
for ($i=1;$i<=$NumEle;$i++)
 {
   $particles=$particles." ".$Comp[$i];
 }
#OUTPUT
open (OUT, ">$posfile");
  print OUT "atom str from ./pwscf.save/data-file-schema.xml \n";
  print OUT "   1.0000000\n";
  #printf (OUT "%22.16f%22.16f%22.16f \n",$ax,$ay,$az);
  #printf (OUT "%22.16f%22.16f%22.16f \n",$bx,$by,$bz);
  #printf (OUT "%22.16f%22.16f%22.16f \n",$cx,$cy,$cz);
  printf (OUT "%22.16f%22.16f%22.16f \n",$latmat_angs->at(0,0),$latmat_angs->at(1,0),$latmat_angs->at(2,0));
  printf (OUT "%22.16f%22.16f%22.16f \n",$latmat_angs->at(0,1),$latmat_angs->at(1,1),$latmat_angs->at(2,1));
  printf (OUT "%22.16f%22.16f%22.16f \n",$latmat_angs->at(0,2),$latmat_angs->at(1,2),$latmat_angs->at(2,2));
  printf (OUT " %s \n",$particles);
  print OUT "Selective Dynamics\n";
  print OUT "Direct\n";
  for($i=0;$i<$NumAtom;$i++){
     printf (OUT "%14.8f%14.8f%14.8f  T  T  T  %s \n", $vecs[$i]->at(0,0),$vecs[$i]->at(1,0),$vecs[$i]->at(2,0)),$Ele[$i];
  }
close(OUT);

