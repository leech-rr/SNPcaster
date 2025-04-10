#!/usr/bin/perl
# "core.full.aln"ファイルから“-”や混合塩基の箇所を省き、"core.full.fas"ファイルに出力する
# 出力1: core.full.fas     --> "core.full.aln"から“-”や混合塩基の箇所を取り除いた配列
# 出力2: core_region.tsv  --> "core.full.aln"から“-”や混合塩基を除いた箇所の位置
# 出力3: matching_list.txt --> core.full.aln(original)とcore.full.fas(extract)の位置の対応表
# $ARGV[0]:"core.full.aln"ファイル名 （引数が無い場合は"core.full.aln"で実行）
use strict;
use warnings;
use File::Basename;

# ファイル名
my $input_core  = "core.full.aln";
my $output_core = "core.full.fas";
my $output_sum  = "core_region.tsv";
my $output_match = "matching_list.txt";

# 混合塩基のリスト
my @mixed_bases = qw (- n b d h k m r s v w y N B D H K M R S V W Y);

# ハッシュ名
my %Bases;
my %N_Bases;
my %Sequences;

# 処理開始
&main;

# メイン処理
sub main {
    # 引数チェック
    &check_parameters;
    
    # fileからseqsequenceを取得する
    my @list = get_sequence($input_core);
    
    # 取得したseqsequenceから混合塩基部分を取り除く
    &get_data($output_sum, $output_match, @list);
    
    # fileに出力する
    &output_data($output_core, @list);
    
    exit;
}

# 引数チェック
sub check_parameters {
    if (defined $ARGV[0]) {
        
        # 変数に代入
        $input_core = $ARGV[0];
        
        my ($name) = basename($input_core) =~ /^(.*?)\.[^.]+$/;
        $output_core = $name . ".fas";
    }
    
    # 既存ファイルを削除
    unlink $output_core if ( -e $output_core );
    unlink $output_sum if ( -e $output_sum );
}

# fileからデータを取得する
sub get_data {
    my ($file_sum, $file_match, @list) = @_;
    open(my $out_1, '>', $file_sum) or die "cannot open > '$file_sum': $!\n";
    open(my $out_2, '>', $file_match) or die "cannot open > '$file_match': $!\n";
    
    my $length = 0;
    
    # 混合塩基の位置を割り出す
    foreach my $strain (@list) {
        
        $length = length $Sequences{$strain};
#        print "$strain : $length\n";
        
        my $counter = 1;
        foreach my $base (split //, $Sequences{$strain}) {
            $N_Bases{$counter}++ if (grep {$_ eq $base} @mixed_bases);
            $counter++;
        }
        $N_Bases{$counter}++;
    }
    
    # 混合塩基以外の位置を出力する
    my $count = 1;
    my ($last_base, $onset, $flag) = (0) x 3;
    
    print $out_1 "Start\tEnd\n";
    print $out_2 "extract\toriginal\n";
    
    for (my $i = 1; $i <= $length + 1; $i++) {
        
        if ($N_Bases{$i}) {
            
            if ($flag) {
                print $out_1 "$onset\t$last_base\n";
                $Bases{$onset} = $last_base;
            }
            $flag = 0;
            
        } else {
            
            print $out_2 "$count\t$i\n";
            $count++;
            
            $onset = $i unless ($flag);
            $flag = 1;
        }
        
        $last_base = $i;
    }
    
    # ハッシュの初期化
    %N_Bases = ();
    
    close $out_1;
    close $out_2;
}

# データを出力する
sub output_data {
    my ($file, @list) = @_;
    open(my $out, '>>', $file) or die "cannot open > '$file': $!\n";
    
    foreach my $strain (@list) {
        
        print $out ">$strain\n";
        
        my $string = "";
        
        # ハッシュ
        for my $start (sort {$a <=> $b} keys %Bases) {
            
            my $length = $Bases{$start} - $start + 1;
            $start -= 1;
            my $part = substr $Sequences{$strain}, $start, $length;
            $string .= $part;
        }
        
        # 80字ずつ区切る
        $string =~ s|(\w{80})|$1\n|g;
        print $out "$string\n";
    }
    
    close $out;
}

# ファイルからsequenceを取得する
sub get_sequence {
    my ($file) = @_;
    open(my $in, '<', $file) or die "cannot open < '$file': $!\n";
    
    my @list = (  );
    my ($name, $sequence) = ("") x 2;
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;
        
        # blank line
        if ($line =~ m/^\s*$/) {
            next;
            
        # comment line
        } elsif ($line =~ m/^\s*#/) {
            next;
            
        # fasta header line
        } elsif ($line =~ m/^>/) {
            
            if ($name) {
                # ハッシュに代入
                $sequence =~ s/\s//g;
                $Sequences{$name} = $sequence;
                
                $sequence = '';
            }
            
            if ($line =~ /\s+/) {
                ($name) = ($line =~ m/^>(.+?)\s+/);
            } else {
                ($name) = ($line =~ m/^>(.+)/);
            }
            
            push (@list, $name);
            
        # keep line
        } else {
            $line =~ s/\s//g;
            $sequence .= $line;
        }
    }
    
    # 最終行分
    $sequence =~ s/\s//g;
    $Sequences{$name} = $sequence;
    
    close $in;
    return @list;
}