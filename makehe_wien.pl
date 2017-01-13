#!/usr/bin/perl
# usage makehe_wien.pl case [-nc]
#                             -nc: non-complex calc

if($ARGV[1] ne "-nc"){
$klist=$ARGV[0].".klist";
$in1c=$ARGV[0].".in1c";

$klisthe=$klist."h";
$in1che=$in1c."h";

open(IN1C,"$in1c") || die "cannot open $in1c\n";
open(IN1CHE,">$in1che");
$KVECTORSFound=0;
while($KVECTORSFound==0)
{
 chomp($line=<IN1C>);
 $_=$line;
 if(/K-VECTORS/){$KVECTORSFound=1;}
 else{print IN1CHE "$line\n";}
}
@values=split(/\ +/,$line);
$values[4]=3.5;
printf (IN1CHE "%s %s %s%7.1f%10.1f      %s %s\n",$values[0],$values[1],$values[2],$values[3],$values[4],$values[5],$values[6]);
close(IN1C);
close(IN1CHE);

open(KLIST,"$klist") || die "cannot open $klist\n";
open(KLISTHE,">$klisthe");

chomp($line=<KLIST>);@values=split(/\ +/,$line);
$values[8]=3.5;$lastterm="        $values[9] $values[10] $values[11] $values[12]  $values[13]  $values[14]  $values[15]";
printf (KLISTHE "%10d%10d%10d%10d%10d%5.1f%5.1f%5.1f%s \n", $values[1],$values[2],$values[3],$values[4],$values[5],$values[6],$values[7],$values[8],$lastterm);
while(!eof(KLIST))
{
 chomp($line=<KLIST>);
 print KLISTHE "$line\n";
}
close(KLIST);
close(KLISTHE);

}

if($ARGV[1] eq "-nc"){
$klist=$ARGV[0].".klist";
$in1c=$ARGV[0].".in1";

$klisthe=$klist."h";
$in1che=$in1c."h";

open(IN1C,"$in1c") || die "cannot open $in1c\n";
open(IN1CHE,">$in1che");
$KVECTORSFound=0;
while($KVECTORSFound==0)
{
 chomp($line=<IN1C>);
 $_=$line;
 if(/K-VECTORS/){$KVECTORSFound=1;}
 else{print IN1CHE "$line\n";}
}
@values=split(/\ +/,$line);
$values[4]=3.5;
printf (IN1CHE "%s %s %s%7.1f%10.1f      %s %s\n",$values[0],$values[1],$values[2],$values[3],$values[4],$values[5],$values[6]);
close(IN1C);
close(IN1CHE);

open(KLIST,"$klist") || die "cannot open $klist\n";
open(KLISTHE,">$klisthe");

chomp($line=<KLIST>);@values=split(/\ +/,$line);
$values[8]=3.5;$lastterm="        $values[9] $values[10] $values[11] $values[12]  $values[13]  $values[14]  $values[15]";
#printf (KLISTHE "%10d%5d%5d%5d%5d%5.1f%5.1f%5.1f%s \n", $values[1],$values[2],$values[3],$values[4],$values[5],$values[6],$values[7],$values[8],$lastterm);
printf (KLISTHE "%10d%10d%10d%10d%10d%5.1f%5.1f%5.1f%s \n", $values[1],$values[2],$values[3],$values[4],$values[5],$values[6],$values[7],$values[8],$lastterm);
while(!eof(KLIST))
{
 chomp($line=<KLIST>);
 print KLISTHE "$line\n";
}
close(KLIST);
close(KLISTHE);
}
 
 
