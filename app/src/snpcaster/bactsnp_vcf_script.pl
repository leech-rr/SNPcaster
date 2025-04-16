#!/usr/bin/perl
# リファレンスとlist内の配列を1塩基ずつ比較し、異なる場合に出力する
# 
# $ARGV[0]: 参照配列
# $ARGV[1]: 比較する配列
use strict;
use warnings;
use Bio::SeqIO;

# variables
my $ref_file = $ARGV[0];
my $fas_file = $ARGV[1];
my $out_file = "result_vcf.txt";
my %Sequences = ();
my %Count = ();


# 処理開始
&main;

# メイン処理
sub main {
    # Ref配列を取得
    my $ref_id = &get_sequences("Ref", $ref_file);
    
    # 配列を取得
    &get_sequences("Fas", $fas_file);
    
    # check & 出力
    &output_data($ref_id, $out_file);
    
    exit;
}

# 配列を取得
sub get_sequences {
    my ($strain, $file) = @_;
    my $in = Bio::SeqIO->new(-file=>$file, -format=>"fasta", -alphabet=>'dna');

    my $id;
    my $sequence_count = 0;
    while (my $seq = $in->next_seq()) {        
        if ($sequence_count > 1) {
            die "Error: Multiple sequences found in [$file]. Only a single sequence is allowed.\n";
        }
        $id = $seq->id;
        # 1文字ずつハッシュ代入
        my $count = 1;
        foreach my $cha (split //, $seq->seq) {
            $Sequences{$strain}{$count} = $cha;
            $count++;
        }
        $sequence_count++;
    }
    return $id;
}

# check & 出力
sub output_data {
    my ($ref_id, $out_file) = @_;
    open(my $out, '>', $out_file) or die "cannot open > '$out_file': $!\n";
    
    # count数
    my $max = scalar( keys(%{$Sequences{"Ref"}}) );
    
    for (my $count = 1; $count <= $max; $count++ ) {
        
        my $cha_Ref = $Sequences{"Ref"}{$count};
        my $cha_Fas = $Sequences{"Fas"}{$count};
        
        next if ($cha_Ref eq $cha_Fas);
        next if ($cha_Ref eq "-" or $cha_Ref eq "N");
        next if ($cha_Fas eq "-" or $cha_Fas eq "N");
        
        print $out "$ref_id\t$count\t" . "." . "\t$cha_Ref\t$cha_Fas\n";
    }
    
    close $out;
}
