#!/usr/bin/perl
$infile="POSCAR";
$bin=$ARGV[0];
open(IN,"$infile")||die "Cannot open $infile ";
 chomp($line=<IN>); @values=split(/ +/, $line);
 chomp($line=<IN>); @values=split(/ +/, $line);
    $fac=$values[1];
 chomp($line=<IN>); @values=split(/ +/, $line);
	$ax=$values[1]*$fac;$ay=$values[2]*$fac;$az=$values[3]*$fac;
 chomp($line=<IN>); @values=split(/ +/, $line);
	$bx=$values[1]*$fac;$by=$values[2]*$fac;$bz=$values[3]*$fac;
 chomp($line=<IN>); @values=split(/ +/, $line);
	$cx=$values[1]*$fac;$cy=$values[2]*$fac;$cz=$values[3]*$fac;
 $volume=$ax*($by*$cz-$bz*$cy)+$ay*($bz*$cx-$bx*$cz)+$az*($bx*$cy-$by*$cx);
 print "[volume]: $volume\n";
 
 $iax=($by*$cz-$bz*$cy)/$volume;$iay=($bz*$cx-$bx*$cz)/$volume;$iaz=($bx*$cy-$by*$cx)/$volume;
 $ibx=($cy*$az-$cz*$ay)/$volume;$iby=($cz*$ax-$cx*$az)/$volume;$ibz=($cx*$ay-$cy*$ax)/$volume;
 $icx=($ay*$bz-$az*$by)/$volume;$icy=($az*$bx-$ax*$bz)/$volume;$icz=($ax*$by-$ay*$bx)/$volume;
 printf ("%12.4f %12.4f %12.4f\n",$iax,$iay,$iaz);
 printf ("%12.4f %12.4f %12.4f\n",$ibx,$iby,$ibz);
 printf ("%12.4f %12.4f %12.4f\n",$icx,$icy,$icz);
 
 $ka=int(sqrt($iax**2+$iay**2+$iaz**2)/$bin);
 $kb=int(sqrt($ibx**2+$iby**2+$ibz**2)/$bin);
 $kc=int(sqrt($icx**2+$icy**2+$icz**2)/$bin);

 $verifier=$ka*$kb*$kc*$volume;
 print "[kmesh]: $ka x $kb x $kc [density]: $verifier Angs^3 \n";
 
 close(IN);
 
 
 
