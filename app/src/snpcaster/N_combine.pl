#!/usr/bin/perl


if($ARGV[2] eq "L"){
	open(OUT2,"> $ARGV[0]"."_n_combine_posi.tsv");
}

if(!$ARGV[1]){
	if($ARGV[1] eq '0'){
		$amount=0;
	}else{
		$amount=100;	
	}
}else{
	$amount=$ARGV[1];	
}

open(OUT,"> $ARGV[0]"."_N".$amount.".fasta");

$n="";
for($i=0;$i<$amount;$i++){
	$n.="N";	
}


$i=0;
$b=0;
@name=();
open(FILE,"$ARGV[0]");
print $ARGV[0];
while(<FILE>){
    chomp;
    if(/>([\S|\s]+)/){
	if($i==0){                     #dereat fasta type
    		#print ">".$ARGV[0]."\n";
		print OUT ">".$ARGV[0]."\n";	
		$i++;
	}else{
		$com="@";
		push @nuc,$com;
	}
     	$name[$b]=$1;
	#print $_."\n";
	$b++;
     }else{
	if(/[ATCGNX]/g || /[atgcnx]/g){
	    @nu=split(//,$_);
	    push @nuc,@nu;                 #input nuc_data 
	}
	
    }
    

}
close(FILE);

print "$trash\n";
#print OUT "$trash";     #output fasta_type

$linecount=0;

$nuc_count=$#nuc;
$a=0;
$fla=0;

for($i=$x;$i<$nuc_count+1;$i++){
	
	    if($nuc[$i] eq "@"){
		
		if($fla==1){
			$end=$i+$amount*$a-$a;
			$fla=0;

			print $amount."\n";
		}	
		print OUT $n;
	    	print $name[$a]."\t".$start."\t".$end."\n";
		print OUT2 $name[$a]."\t".$start."\t".$end."\n";	
		$a++;
	     }else{
		if($fla==0){
			$start=$i+$a*$amount-$a+1;
			$fla=1;
		}   
		print OUT $nuc[$i];                  #output nuc_type
		
	     }
}

$end=$i+$amount*$a-$a;

print $name[$a]."\t".$start."\t".$end."\n";
print OUT2 $name[$a]."\t".$start."\t".$end."\n";

close(OUT);
close(OUT2);
