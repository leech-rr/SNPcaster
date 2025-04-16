#!/usr/bin/perl
# coverageが基準値以上の菌株リストを作成
# $ARGV[0]: 出力ファイル カバレッジが基準値以上の菌株リスト
# $ARGV[1]: 基準とするカバレッジの値 (デフォルトは40)
use strict;
use warnings;

# ファイル名
my $in_file = "coverage.txt";
my $out_file = $ARGV[0];

# 変数
my $min_coverage = 40;
$min_coverage = $ARGV[1] if ($ARGV[1]);


# 処理開始
&main;

# メイン処理
sub main {
    # data取得、出力
    &get_data($in_file, $out_file);
    
    # 処理の終了
    print "done!\n";
    exit;
}

# dataを取得
sub get_data {
    my ($in_file, $out_file) = @_;
    open(my $in, '<', $in_file) or die "cannot open < '$in_file': $!\n";
    open(my $out, '>', $out_file) or die "cannot open > '$out_file': $!\n";
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;
        next if ($line eq "");
        
        # 1行目はスキップ
        next if ($. == 1);
        
        # タブで分割
        my ($strain, $coverage) = (split /\t/, $line) [0, 3];
        
        $coverage += 0;
        
        if ($coverage >= $min_coverage) {
            print $out "$strain\n";
        }
    }
    
    close $in;
    close $out;
}
