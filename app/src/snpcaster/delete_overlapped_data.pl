#!/usr/bin/perl
# 引数で渡されたファイルから、重複データを除去し出力する
# $ARGV[0]:処理対象もととなるのファイル名
use strict;
use warnings;

# ファイル名
my $input_file;
my $output_file;

# ハッシュ名
my %Values = ();

# 変数
my $dir = './';


# 処理開始
&main;

# メイン処理
sub main {
    # 引数チェック
    &check_parameters;
    
    # 出力ファイル名
    my ($a, $b) = ($input_file =~ /^(.*)(\..*)$/);
    $output_file = $a . "_without_overlap" . $b;
    
    # ファイルからデータ取得
    &collect_data($input_file);
    
    # ファイル出力
    &output_file($output_file);
    
    # 処理の終了
    print "done!\n";
    exit;
}

# 引数チェック
sub check_parameters {
    if (@ARGV < 1) {
        print 'Usage: perl delete_overlapped_data.pl file.txt', "\n";
        print '  $1  File Name', "\n";
        exit;
    }
    
    # 引数の代入
    $input_file = $ARGV[0];
}

# ファイルからデータを取得する
sub collect_data {
    my ($file) = @_;
    
    open(my $in, '<', $file) or die "cannot open < '$file': $!\n";
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/^\s+//;          # 行頭のスペースを削除
        next if ($line =~ /^\D+/);  # 数字で始まらない行は飛ばす
        
        # タブで分割
        my ($start, $end) = (split /\t/, $line) [0, 1];
        
        # Start値の方が小さい時は入れかえる
        if ($start > $end) {
            my $a = $start;
            $start = $end;
            $end = $a;
        }
        
        # ハッシュに代入
        if (exists $Values{$start}) {
            next if ($Values{$start} > $end);
            $Values{$start} = $end;
            
        } else {
            $Values{$start} = $end;
        }
    }
    close $in;
}

# 処理: ファイルに出力する
sub output_file {
    my ($file) = @_;
    open(my $out, '>', $file) or die "cannot open > '$file': $!\n";
    
    # タイトル行を出力
    print $out "Start\tEnd\n";
    
    # ループ内で使う変数
    my ($onset, $maxEnd, $flag) = (0) x 5;
    my $count = 0;
    
    # ハッシュ
    for my $start (sort {$a <=> $b} keys %Values) {
        my $end = $Values{$start};
        
        # 1行目
        if ($count == 0) {
            $onset = $start;
            $maxEnd = $end;
            $count++;
            next;
        }
        
        # チェック
        if ($start > $maxEnd + 1) {
            # MaxEndを出力
            print $out "$onset\t$maxEnd\n";
            
            # リセット
            $onset = $start;
            $maxEnd = $end;
            
        } else {
            $maxEnd = $end if ($maxEnd < $end);
        }
    }
    
    # 最終行を出力
    print $out "$onset\t$maxEnd\n";
    
    close $out;
}

