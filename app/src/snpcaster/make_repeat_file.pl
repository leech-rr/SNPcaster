#!/usr/bin/perl
# tandems、seq.repeats、seq_seq.coordsのファイルを元に、mask.tsvファイルを作成する
# $ARGV[0]:参照配列名(repeats.tsvの1列目として使う)
use strict;
use warnings;

# ファイル名
my $tandem_file = 'tandems';
my $exact_file = 'seq.repeats';
my $repeat_file = 'seq_seq.coords';
my $mask_file = 'mask.tsv';
my $phage_file = 'phage.tsv';

# ハッシュ名
my %Repeat = ();

# 変数
my $dir = './';
my $tag;


# 処理開始
&main;

# メイン処理
sub main {
    # 引数チェック
    &check_parameters;
    
    # 'tandems'ファイルからデータを取得
    &collect_tandems($tandem_file);
    
    # 'seq.repeats'ファイルからデータを取得
    &collect_exact($exact_file);
    
    # 'seq_seq.coords'ファイルからデータ取得
    &collect_repeat;
    
    # 'phage.tsv'ファイルからデータ取得
    my @phage_files = get_file_names($phage_file);
    if (@phage_files) {
        &collect_phage($phage_file);
    }
    
    # ファイル出力
    &output_file;
    
    # 処理の終了
    print "done!\n";
    exit;
}

# 引数チェック
sub check_parameters {
    unless (@ARGV == 1) {
        print "Set parameter!\n";
        exit;
    }
    
    # 引数の代入
    $tag = $ARGV[0];
}

# 処理: 'tandems'ファイルからデータを取得する
sub collect_tandems {
    my ($file) = @_;
    
    open(my $in, '<', $file) or die "cannot open < '$file': $!\n";
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/^\s+//;          # 行頭のスペースを削除
        next if ($line =~ /^\D+/);  # 数字で始まらない行は飛ばす
        
        # 複数スペースで分割
        my ($start, $extent) = (split /\s+/, $line) [0, 1];
        
        # Endを求める:End = Start + Extent - 1
        my $end = $start + $extent - 1;
        
        # Start値の方が小さい時は入れかえる
        if ($start > $end) {
            my $a = $start;
            $start = $end;
            $end = $a;
        }
        
        # ハッシュに代入
        if (exists $Repeat{$start}) {
            next if ($Repeat{$start} > $end);
            $Repeat{$start} = $end;
            
        } else {
            $Repeat{$start} = $end;
        }
    }
    close $in;
}

# 処理: 'seq.repeats'ファイルからデータを取得する
sub collect_exact {
    my ($file) = @_;
    
    open(my $in, '<', $file) or die "cannot open < '$file': $!\n";
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/^\s+//;          # 行頭のスペースを削除
        next if ($line =~ /^\D+/);  # 数字で始まらない行は飛ばす
        
        # 複数スペースで分割
        my ($start1, $start2, $length) = (split /\s+/, $line) [0, 1, 2];
        
        # Endを求める
        # End1: Start1 + Length - 1
        my $end1 = $start1 + $length - 1;
        
        my $end2;
        if ($start2 =~ /r/) {
            # End2 r有: Start2 - Length + 1
            $start2 =~ s/r//;   # rを削除
            $end2 = $start2 - $length + 1;
            
        } else {
            # End2 r無: Start2 + Length - 1
            $end2 = $start2 + $length - 1;
        }
        
        # Start値の方が小さい時は入れかえる
        if ($start1 > $end1) {
            my $a = $start1;
            $start1 = $end1;
            $end1 = $a;
        }
        if ($start2 > $end2) {
            my $a = $start2;
            $start2 = $end2;
            $end2 = $a;
        }
        
        # ハッシュに代入
        if (exists $Repeat{$start1}) {
            $Repeat{$start1} = $end1 if ($end1 > $Repeat{$start1});
        } else {
            $Repeat{$start1} = $end1;
        }
        
        if (exists $Repeat{$start2}) {
            next if ($Repeat{$start2} > $end2);
            $Repeat{$start2} = $end2;
        } else {
            $Repeat{$start2} = $end2;
        }
    }
    close $in;
}

# 処理: 'seq_seq.coords'ファイルからデータを取得する
sub collect_repeat {
    open(my $in, '<', $repeat_file) or die "cannot open < '$repeat_file': $!\n";
    
    OUTER: while (my $line = <$in>) {
        chomp $line;
        $line =~ s/^\s+//;          # 行頭のスペースを削除
        next if ($line =~ /^\D+/);  # 数字で始まらない行は飛ばす
        next if ($line eq "");      # 空行は飛ばす
        
        # 複数スペースで分割
        my ($start, $end, $e2, $len2) = (split /\s+/, $line) [0, 1, 4, 7];
        
        # E1とE2とLEN2の値が同じ行はスキップ
        if ($e2 && $len2) {
            next if ($end == $e2 && $end == $len2);
        }
        
        # Start値の方が小さい時は入れかえる
        if ($start > $end) {
            my $a = $start;
            $start = $end;
            $end = $a;
        }
        
        # ハッシュに代入
        if (exists $Repeat{$start}) {
            next if ($Repeat{$start} > $end);
            $Repeat{$start} = $end;
            
        } else {
            $Repeat{$start} = $end;
        }
    }
    close $in;
}

# ディレクトリ内にあるファイル名を取得する
sub get_file_names {
    my $file_name = shift;
    opendir(my $dh, $dir) || die "can't opendir $dir: $!";
    
    # 特定の名前のファイルだけをリストに保存
    my @files = sort grep { /$file_name$/i } readdir($dh);
    
    closedir $dh;
    
    return @files;
}

# 処理: 'mask.tsv'ファイルからデータを取得する
sub collect_phage {
    my ($file) = @_;
    
    open(my $in, '<', $file) or die "cannot open < '$file': $!\n";
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;   # 末尾の改行を削除
        
        # タブで分割
        my ($chr, $start, $end) = (split /\t/, $line) [0..2];
        
        # $startが数字以外の場合はスキップ
        next if ($start =~ m/\D+/ || $end =~ m/\D+/);
        
        # ハッシュに代入
        if (exists $Repeat{$start}) {
            next if ($Repeat{$start} > $end);
            $Repeat{$start} = $end;
            
        } else {
            $Repeat{$start} = $end;
        }
    }
    close $in;
}


# 処理: ファイルに出力する
sub output_file {
    open(my $out, '>', $mask_file) or die "cannot open > '$mask_file': $!\n";
    
    # タイトル行を出力
    print $out "Tag\tStart\tEnd\n";
    
    # ループ内で使う変数
    my ($preStart, $preEnd, $onset, $maxEnd, $flag) = (0) x 5;
    my $count = 0;
    
    # ハッシュ
    for my $start (sort {$a <=> $b} keys %Repeat) {
        my $end = $Repeat{$start};
        
        # 1行目
        if ($count == 0) {
            $preStart = $start;
            $onset = $start;
            $preEnd = $end;
            $maxEnd = $end;
            $count++;
            next;
        }
        
        # チェック
        if ($start > $maxEnd + 1) {
            # MaxEndを出力
            print $out "$tag\t$onset\t$maxEnd\n";
            
            # リセット
            $preStart = $start;
            $onset = $start;
            $preEnd = $end;
            $maxEnd = $end;
            
        } else {
            $maxEnd = $end if ($maxEnd < $end);
        }
    }
    
    # 最終行を出力
    print $out "$tag\t$onset\t$maxEnd\n";
    
    close $out;
}

