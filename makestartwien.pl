#!/usr/bin/perl
# usage: makestartwien.pl case [-nc]
#                         -nc: non-complex calc. 
$outfile="$ARGV[0].csh";


if ($ARGV[1] ne "-nc"){

open(OUT,">$outfile");
print OUT "#!/bin/csh \n";
print OUT "#\$ -cwd  \n";
print OUT "submit.pl \n";
print OUT "rsh machx \"cd \$PWD;x dstart -c >> log1\" & \n";
print OUT "rsh machx \"cd \$PWD;x dstart -c -up>> log2\" & \n";
print OUT "rsh machx \"cd \$PWD;x dstart -c -dn>> log3\" & \n";
print OUT "wait \n";
print OUT "cp $ARGV[0].clmsum $ARGV[0]_grnd/$ARGV[0]_grnd.clmsum \n";
print OUT "cp $ARGV[0].clmup $ARGV[0]_grnd/$ARGV[0]_grnd.clmup \n";
print OUT "cp $ARGV[0].clmdn $ARGV[0]_grnd/$ARGV[0]_grnd.clmdn \n";
print OUT "runsp_lapw -p -ec 0.001 -i 60\n";
print OUT "save_lapw save \n";
print OUT "cp $ARGV[0].in1ch $ARGV[0].in1c \n";
print OUT "cp $ARGV[0].klisth $ARGV[0].klist \n";
print OUT "x lapw1 -c -p -up  \n";
print OUT "x lapw1 -c -p -dn  \n";
print OUT "x lapw2 -c -p -qtl -up \n";
print OUT "x lapw2 -c -p -qtl -dn \n";
close(OUT);
system "chmod 744 $ARGV[0].csh";
system "makehe_wien.pl $ARGV[0]";

}

if($ARGV[1] eq "-nc"){

open(OUT,">$outfile");
print OUT "#!/bin/csh \n";
print OUT "#\$ -cwd  \n";
print OUT "submit.pl \n";
print OUT "rsh machx \"cd \$PWD;x dstart  >> log1\" & \n";
print OUT "rsh machx \"cd \$PWD;x dstart  -up>> log2\" & \n";
print OUT "rsh machx \"cd \$PWD;x dstart  -dn>> log3\" & \n";
print OUT "wait \n";
print OUT "cp $ARGV[0].clmsum $ARGV[0]_grnd/$ARGV[0]_grnd.clmsum \n";
print OUT "cp $ARGV[0].clmup $ARGV[0]_grnd/$ARGV[0]_grnd.clmup \n";
print OUT "cp $ARGV[0].clmdn $ARGV[0]_grnd/$ARGV[0]_grnd.clmdn \n";
print OUT "runsp_lapw -p -ec 0.001 -i 60\n";
print OUT "save_lapw save \n";
print OUT "cp $ARGV[0].in1h $ARGV[0].in1 \n";
print OUT "cp $ARGV[0].klisth $ARGV[0].klist \n";
print OUT "x lapw1  -p -up  \n";
print OUT "x lapw1  -p -dn  \n";
print OUT "x lapw2  -p -qtl -up \n";
print OUT "x lapw2  -p -qtl -dn \n";
close(OUT);
system "chmod 744 $ARGV[0].csh";
system "makehe_wien.pl $ARGV[0] -nc";

}
