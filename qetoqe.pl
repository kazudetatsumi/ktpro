#!/usr/bin/perl
#This script read the qe output xml file and create a qe new input file with the calculated atomic structure
#The generated input file is used to obtain further accurate atomic/electronic structure calculation.
#Kazuyoshi TATSUMI 2019 12 18
#latmat was transposed to work correctly.
use XML::TreePP;
use PDL;
use PDL::Matrix;
$atoangs = 0.529177210903;
$tpp = XML::TreePP->new();
$output_str = $tpp->parsefile( "./pwscf.save/data-file-schema.xml" )->{'qes:espresso'}->{output}->{atomic_structure};
$outfile="outfile.in";

$species = $tpp->parsefile( "./pwscf.save/data-file-schema.xml" )->{'qes:espresso'}->{output}->{atomic_species}{species};
for ($cnt = 0; $cnt < @$species; $cnt++) {
    $Elename[$cnt] = $species->[$cnt]{"-name"};
	$AM[$cnt] = $species->[$cnt]->{mass};
	$pseu[$cnt] = $species->[$cnt]->{pseudo_file};
}
$NumEle = $cnt;

$cell = $output_str->{cell};
@values1 = split(/\ +/, $cell->{a1}); 
@values2 = split(/\ +/, $cell->{a2}); 
@values3 = split(/\ +/, $cell->{a3}); 
$latmat = transpose mpdl [$values1[0], $values1[1], $values1[2]],[$values2[0], $values2[1], $values2[2]], [$values3[0], $values3[1], $values3[2]];
#$latmat = mpdl [$values1[0], $values1[1], $values1[2]],[$values2[0], $values2[1], $values2[2]], [$values3[0], $values3[1], $values3[2]];
$inv_latmat =  matinv $latmat;
$latmat_angs = $latmat * $atoangs;
print $latmat_angs;

$atoms = $output_str->{atomic_positions}{atom};
for ($cnt = 0; $cnt < @$atoms; $cnt++) {
	$Ele[$cnt] = $atoms->[$cnt]{"-name"};
	@values = split(/\ +/, $atoms->[$cnt]{"#text"});
	$vecs[$cnt] = $inv_latmat x vpdl [$values[0], $values[1], $values[2]]; 
}
$NumAtom = $cnt;


#OUTPUT
if($ARGV[0] eq "-vc"){
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
  printf (OUT "    nat = %d\n", $NumAtom);
  printf (OUT "    ntyp = %d \n", $NumEle);
  print OUT "    ecutwfc = 70.0\n";
  print OUT "    occupations='smearing',\n";
  print OUT "    smearing='mv',\n";
  print OUT "    degauss=0.01,\n";
  print OUT "/\n";
  print OUT "&electrons\n";
  print OUT "    diagonalization = 'david'\n";
  print OUT "    conv_thr = 1.0d-9\n";
  print OUT "    mixing_beta = 0.4\n";
  print OUT "    mixing_ndim = 45\n";
  print OUT "/\n";
  print OUT "&ions\n";
  print OUT "/\n";
  print OUT "&cell\n";
  print OUT "/\n";
  print OUT "ATOMIC_SPECIES\n";
  for($i=0;$i<$NumEle;$i++)
   {
     printf (OUT " %s %12.8f  %s\n", $Elename[$i], $AM[$i], $pseu[$i]);
  }
  print OUT "ATOMIC_POSITIONS crystal\n";
  for($i=0;$i<$NumAtom;$i++)
   {
     printf (OUT " %s %12.8f%12.8f%12.8f  \n", $Ele[$i], $vecs[$i]->at(0,0),$vecs[$i]->at(1,0),$vecs[$i]->at(2,0));
  }
  print OUT "CELL_PARAMETERS angstrom \n";
  printf (OUT "%22.16f%22.16f%22.16f \n",$latmat_angs->at(0,0),$latmat_angs->at(1,0),$latmat_angs->at(2,0));
  printf (OUT "%22.16f%22.16f%22.16f \n",$latmat_angs->at(0,1),$latmat_angs->at(1,1),$latmat_angs->at(2,1));
  printf (OUT "%22.16f%22.16f%22.16f \n",$latmat_angs->at(0,2),$latmat_angs->at(1,2),$latmat_angs->at(2,2));
  print OUT "K_POINTS automatic\n";
  print OUT " 5 5 3 0 0 0\n";
close(OUT);
}
else{
open (OUT, ">$outfile");
  print OUT "&control\n";
  print OUT "    calculation = 'scf'\n";
  print OUT "    tprnfor = .true.\n";
  print OUT "    tstress = .true.\n";
  print OUT "    pseudo_dir = '/home/kazu/code/q-e-qe-6.4.1/pseudo/'\n";
  print OUT "/\n";
  print OUT "&system\n";
  print OUT "    ibrav = 0\n";
  printf (OUT "    nat = %d\n", $NumAtom);
  printf (OUT "    ntyp = %d \n", $NumEle);
  print OUT "    ecutwfc = 70.0\n";
  print OUT "/\n";
  print OUT "&electrons\n";
  print OUT "    diagonalization = 'david'\n";
  print OUT "    conv_thr = 1.0d-9\n";
  print OUT "/\n";
  print OUT "ATOMIC_SPECIES\n";
  for($i=0;$i<$NumEle;$i++)
   {
     printf (OUT " %s %12.8f  %s\n", $Elename[$i], $AM[$i], $pseu[$i]);
  }
  print OUT "ATOMIC_POSITIONS crystal\n";
  for($i=0;$i<$NumAtom;$i++)
   {
     printf (OUT " %s %12.8f%12.8f%12.8f  \n", $Ele[$i], $vecs[$i]->at(0,0),$vecs[$i]->at(1,0),$vecs[$i]->at(2,0));
  }
  print OUT "CELL_PARAMETERS angstrom \n";
  printf (OUT "%22.16f%22.16f%22.16f \n",$latmat_angs->at(0,0),$latmat_angs->at(1,0),$latmat_angs->at(2,0));
  printf (OUT "%22.16f%22.16f%22.16f \n",$latmat_angs->at(0,1),$latmat_angs->at(1,1),$latmat_angs->at(2,1));
  printf (OUT "%22.16f%22.16f%22.16f \n",$latmat_angs->at(0,2),$latmat_angs->at(1,2),$latmat_angs->at(2,2));
  print OUT "K_POINTS automatic\n";
  print OUT " 8 8 8 1 1 1\n";
close(OUT);
}
