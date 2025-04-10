#!/usr/bin/perl
# core_extract.pl を実行し作成された core_region.tsv から、recombination.tsv とmask.tsvの範囲を取り除く
# 
# $ARGV[0]: core_region.tsv か 同じ形式のファイル (startとendをタブ区切り)
# $ARGV[1]: recombination.tsv か mask.tsv
# $ARGV[2]: タグ (recombination.tsv, mask.tsv の1列目)（指定が無い場合は全て）
use strict;
use warnings;
use File::Basename;

# ファイル名
my $input_sum;
my $input_recombination;
my $input_tag;
my $output_file;

# ハッシュ名
my %Sums = ();
my %Rcms = ();
my %Data = ();

# 処理開始
&main;

# メイン処理
sub main {
    # 引数チェック
    &check_parameters;
    
    # 出力ファイル名
    my ($a, $b) = (basename($input_sum) =~ /^(.*)(\..*)$/);
    $output_file = $a . "_without_recombination" . $b;
    
    # データ取得
    &get_data_sum($input_sum);
    &get_data_rcm($input_recombination);
    
    &check_overlap;
    
    # 出力
    &output_data ($output_file);
    
    print "done!\n";
    exit;
}

# 引数チェック
sub check_parameters {
#    print '$ARGV[1]\n';
#    print '@ARGV\n';
    if (@ARGV < 2) {
        print 'Usage: perl core_rcm/core_remove_rep_rcm.pl.pl core_region.tsv recombination.tsv [tag]', "\n";
        print '  $1  core_region.tsv', "\n";
        print '  $2  recombination.tsv or mask.tsv', "\n";
        print '  $3  Tag from recombination.tsv or mask.tsv (All tags by default)', "\n";
        exit;
    }
    
    # 変数代入
    $input_sum = $ARGV[0];
    $input_recombination = $ARGV[1];
    $input_tag = $ARGV[2] if (defined $ARGV[2]);
}

# データを取得する
sub get_data_sum {
    my ($file) = @_;
    open(my $in, '<', $file) or die "cannot open < '$file': $!\n";
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;
        next if ($. == 1);
        next if ($line eq "");
        
        # タブで分割
        my ($start, $end) = split (/\t/, $line, 2);
        $Sums{$start} = $end;
        
        my $len = $end - $start;
        for (my $i = 0; $i <= $len; $i++) {
            $Data{$start + $i}++;
        }
    }
    close $in;
}

# データを取得する
sub get_data_rcm {
    my ($file) = @_;
    open(my $in, '<', $file) or die "cannot open < '$file': $!\n";
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;
        next if ($. == 1);
        next if ($line eq "");
        
        # タブで分割
        my ($tag, $start, $end) = split (/\t/, $line, 3);
        next if ($input_tag && $tag ne $input_tag);
        $Rcms{$start} = $end;
    }
    close $in;
}

# 重複チェック
sub check_overlap {
    # Sumsハッシュ
    foreach my $S1 (sort {$a <=> $b} keys %Sums) {
        my $S2 = $Sums{$S1};
        
        # Rcmsハッシュ
        foreach my $start (sort {$a <=> $b} keys %Rcms) {
            my $end = $Rcms{$start};
            
            next if ($start < $S1 and $end < $S1);
            next if ($S2 < $start and $S2 < $end);
            my $len = $end - $start;
            
            # 重複要素を配列に追加
            for (my $i = 0; $i <= $len; $i++) {
                my $val = $start + $i;
                delete $Data{$val} if (exists($Data{$val}));
            }
        }
    }
}

# データ出力
sub output_data {
    my ($out_file) = @_;
    open(my $out, '>', $out_file) or die "cannot open > '$out_file': $!\n";
    print $out "Start\tEnd\n";
    
    # ループ内で使う変数
    my ($onset, $max, $count) = (0) x 3;
    
    foreach my $no (sort {$a <=> $b} keys %Data) {
        # 1つ目
        if ($count == 0) {
            $max = $no;
            $onset = $no;
            $count++;
            next;
        }
        
        if ($no == $max + 1) {
            $max++;
            
        } else {
            print $out "$onset\t$max\n";
            
            $max = $no;
            $onset = $no;
        }
    }
    
    # 最終行を出力
    print $out "$onset\t$max\n";
    close $out;
}
