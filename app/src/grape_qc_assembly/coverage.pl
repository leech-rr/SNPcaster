#!/usr/bin/perl
# coverageを算出し、ファイル出力する
# $ARGV[0]:菌株リスト名
# $ARGV[1]:genomeサイズ (Mbp)
# $ARGV[2]:拡張子("_12s.fastq")
use strict;
use warnings;

# ファイル名
my $list_file;
my $out_file = "coverage.txt";

# 変数
my $size;
my $extension;


# 処理開始
&main;

# メイン処理
sub main {
    # 引数を代入
    $list_file = $ARGV[0];
    $size = $ARGV[1] * 1000000;
    $extension = $ARGV[2];

    print "[Coverage.pl]\n";
    print "Genome size for coverage calculation: $ARGV[1] Mbp\n";

    my @list = get_list($list_file);
    
    open(my $out, '>', $out_file) or die "cannot open > '$out_file': $!\n";
    print $out "strain\tno. of reads\ttotal read length\tcoverage\n";
    
    foreach my $strain (@list) {
        
        my $file = $strain.$extension;
        
        # ファイルからデータを取得
        my ($no, $lenght, $cov) = get_data($file, $size);
        
        # 出力
        print $out "$strain\t$no\t$lenght\t$cov\n";
    }
    
    close $out;
    
    # 処理の終了
    print "done!\n";
    exit;
}

# listファイルからファイル名を取得
sub get_list {
    my ($list_file) = shift;
    open(my $in, '<', $list_file) or die "cannot open < '$list_file': $!\n";
    
    my @list;
    
    while (my $line = <$in>) {
        # 末尾の改行を削除
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;
        next if ($line eq "");
        # タブで分割
        my $strain = (split /\t/, $line)[0];
        print $strain;
        push (@list, $strain);
    }
    
    return @list;
    close $in;
}

sub get_data {
    my ($in_file, $size) = @_;
    open(my $in, '<', $in_file) or die "cannot open < '$in_file': $!\n";
    
    my $no = 0;
    my $length = 0;
    
    while (my $line = <$in>) {
        
        next if ($. % 4 != 2);
        
        # 末尾の改行を削除
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;
        next if ($line eq "");
        
        $no++;
        $length += length($line);
    }
    
    my $cov = $length / $size;
    
    return ($no, $length, $cov);
    close $in;
}
