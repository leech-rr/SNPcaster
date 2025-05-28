#!/usr/bin/perl
# 隣接SNPを削除し、mask.tsvのデータと被らないデータだけを並び替え出力する
# 'mask.tsv'と'snp_position.csv'が必要
# 
# $ARGV[0]: 隣接SNPの間隔
# $ARGV[1]: Repeat ファイル名
# $ARGV[2]: indelが含まれているファイルを出力する場合は、1を入力する
use strict;
use warnings;
use File::Copy;
use List::Util qw(any);
use Getopt::Long;

# Command args
my ($help, $mask_region_file, $gap, $snp_file);    # placeholder for command line options
$gap = 0;

# Output files
# processing order is:
# 1. Create snp_position_sample_only.csv ($output_1)
# 2. Remove cluster SNP ($output_2-4) but will be skipped if $gap is 0
# 3. Remove masking region ($output_5-7)
# 4. Create final_snp.fas ($output_8)
my $output_1 = 'snp_position_sample_only.csv';   # $snp_fileのうち、サンプル間で差異があるSNPのみを出力したファイル
my $output_2 = 'snp_position_without_clusterSNP.csv';   # 隣接SNPを除いたファイル
my $output_3 = 'snp_position_without_clusterSNP_sample_only.csv';   # 隣接SNPを除いたファイル
my $output_4 = 'removed_clusterSNP.csv';               # 隣接SNPを含んだファイル
my $output_5 = 'snp_position_after_masking.csv';       # maskなど指定領域を除いたファイル
my $output_6 = 'snp_position_after_masking_sample_only.csv';       # snp_position_without_specified_region.csvのうち、サンプル間で差異があるSNPのみを出力したファイル
my $output_7 = 'masked_region.csv';       # 除かれた領域のSNPファイル
my $output_8 = 'snp_position_final.csv';
my $output_9 = 'final_snp.fas';

# 変数
my $title;

# ハッシュ名
my %Snps = ();
my %Mask = ();
my %Delete = ();


# 処理開始
&main;

# メイン処理
sub main {
    &check_parameters;
    
    &process_sequence(
        $snp_file,
        $mask_region_file,
        $gap,
        $output_1,
        $output_2, 
        $output_3, 
        $output_4, 
        $output_5, 
        $output_6,
        $output_7,
        $output_8,
    );
    # 縦横変換
    &make_transpose($output_8, $output_9);
    
    exit 0;
}

sub print_usage {
    print "Usage: $0 [-h] [-d mask_file -n mask_name] [-g gap] snp_file.csv\n";
    print '  $1  Input SNP CSV file', "\n";
    print '  -c  Cluster SNP max distance(for 0, no cluster SNP will be removed. Default: 0)', "\n";
    print '  -d  Mask region tsv file path"', "\n";
}

# 引数チェック
sub check_parameters {
    GetOptions(
        'h'   => \$help,
        'd=s' => \$mask_region_file,
        'c=i' => \$gap,
    ) or do {
        print "Error: Invalid option.\n";
        print_usage();
        exit 1;
    };
    if ($help) {
        print_usage();
        exit 0;
    }

    if (@ARGV < 1) {
        print "Error: Please specify the input SNP csv file.\n";
        print_usage();
        exit 1;
    }
    $snp_file = $ARGV[0];

    # Print args
    print "Input SNP file: $snp_file\n";
    print "Mask region file: ", defined $mask_region_file ? $mask_region_file : '-', "\n";
    print "Cluster SNP max distance: $gap\n";
}

# snpファイルからデータ抽出
sub get_snp {
    my ($in_file) = @_;
    open(my $in, '<', $in_file) or die "cannot open < '$in_file': $!\n";
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;   # 末尾の改行を削除
        next if ($line eq "");
        
        # 1行目、タイトルの設定
        if ($. == 1) {
            $title = $line;
            next;
        }
        
        my @cols = split /,/, $line;
        my ($chrom, $position, $ref, @samples) = (@cols[0..2], @cols[3..$#cols]);
        $Snps{$chrom} = [] unless exists $Snps{$chrom};
        push @{$Snps{$chrom}}, {line => $line, chrom => $chrom, position => $position, ref => $ref, samples => \@samples};
    }
    close $in;

    # acsending order by position
    foreach my $chrom (keys %Snps) {
        @{$Snps{$chrom}} = sort { $a->{position} <=> $b->{position} } @{$Snps{$chrom}};
    }
}

sub is_sample_snp_different {
    my ($snp_hash_ref) = @_;
    my $samples = $snp_hash_ref->{samples};
    my $first = $samples->[0];
    return any { $_ ne $first } @$samples;
}

sub is_cluster_snp {
    my ($cur_snp_hash_ref, $next_snp_hash_ref) = @_;
    my $cur_pos = $cur_snp_hash_ref->{position};
    my $next_pos = $next_snp_hash_ref->{position};
    if ($next_pos - $cur_pos > $gap) {
        return 0;
    }
    my $cur_ref = $cur_snp_hash_ref->{ref};
    my $next_ref = $next_snp_hash_ref->{ref};
    my $cur_samples = $cur_snp_hash_ref->{samples};
    my $next_samples = $next_snp_hash_ref->{samples};
    for my $i (0 .. $#{ $cur_samples }) {
        # Check if the sample base for both of the positions are SNP.
        # SNP means the sample base is different from the reference.
        if ($cur_samples->[$i] ne $cur_ref
            && $next_samples->[$i] ne $next_ref) {
            return 1;
        }
    }
    return 0;
}

sub is_in_mask_region {
    my ($chrom, $snp_hash_ref) = @_;
    my $position = $snp_hash_ref->{position};
    foreach my $region (@{$Mask{$chrom}}) {
        my $start = $region->{start};
        my $end = $region->{end};
        if ($position < $start) {
            return 0;
        }
        if ($position <= $end) {
            return 1;
        }
    }
    return 0;
}

sub process_sequence {
    my $arg_count = scalar(@_);
    unless ($arg_count == 11) {
        die "Error: Invalid number of arguments ($arg_count). Expected 11.\n";
    }

    # Essential arguments
    my (
        $snp_file,
        $mask_region_file,
        $gap,
        $no_mask_sample_only_file, 
        $wo_cluster_snp_file, 
        $wo_cluster_snp_sample_only_file, 
        $removed_cluster_file,
        $w_mask_file,
        $w_mask_sample_only_file,
        $masked_snp_file,
        $final_snp_file
    ) = @_;

    &get_snp($snp_file);
    open(my $out_no_mask_sample_only, '>', $no_mask_sample_only_file) or die "cannot open > '$no_mask_sample_only_file': $!\n";
    open(my $out_final_snp, '>', $final_snp_file) or die "cannot open > '$final_snp_file': $!\n";
    print $out_no_mask_sample_only "$title\n";
    print $out_final_snp "$title\n";

    my $remove_cluster_snp = $gap > 0 ? 1 : 0;
    my $out_wo_cluster_snp = undef;
    my $out_wo_cluster_snp_sample_only = undef;
    my $out_removed_cluster = undef;
    if($remove_cluster_snp){
        open($out_wo_cluster_snp, '>', $wo_cluster_snp_file) or die "cannot open > '$wo_cluster_snp_file': $!\n";
        open($out_wo_cluster_snp_sample_only, '>', $wo_cluster_snp_sample_only_file) or die "cannot open > '$wo_cluster_snp_sample_only_file': $!\n";
        open($out_removed_cluster, '>', $removed_cluster_file) or die "cannot open > '$removed_cluster_file': $!\n";
        print $out_wo_cluster_snp "$title\n";
        print $out_wo_cluster_snp_sample_only "$title\n";
        print $out_removed_cluster "$title\n";
    }

    my $do_masking = $mask_region_file ? 1 : 0;
    my $out_w_mask = undef;
    my $out_w_mask_sample_only = undef;
    my $out_masked_snp = undef;
    if ($do_masking) {
        &get_rep($mask_region_file);
        open($out_w_mask, '>', $w_mask_file) or die "cannot open > '$w_mask_file': $!\n";
        open($out_w_mask_sample_only, '>', $w_mask_sample_only_file) or die "cannot open > '$w_mask_sample_only_file': $!\n";
        open($out_masked_snp, '>', $masked_snp_file) or die "cannot open > '$masked_snp_file': $!\n";
        print $out_w_mask "$title\n";
        print $out_w_mask_sample_only "$title\n";
        print $out_masked_snp "$title\n";
    }

    my @fields = split /,/, $title;
    foreach my $chrom (keys %Snps) {
        my $clustered_with_prev_snp = 0;
        my $clustered_with_next_snp = 0;
        for my $i (0 .. $#{ $Snps{$chrom} }) {
            my $snp_hash_ref = $Snps{$chrom}->[$i];
            my $line = $snp_hash_ref->{line};
            my $sample_snp_different = is_sample_snp_different($snp_hash_ref);

            if ($sample_snp_different) {
                print $out_no_mask_sample_only "$line\n";
            }

            if ($remove_cluster_snp){
                my $next_snp_hash_ref = $Snps{$chrom}->[$i + 1];
                if($next_snp_hash_ref
                    && is_cluster_snp($snp_hash_ref, $next_snp_hash_ref)
                ){
                    $clustered_with_next_snp = 1;
                }
                if ($clustered_with_prev_snp
                    || $clustered_with_next_snp
                ) {
                    print $out_removed_cluster "$line\n";
                    $clustered_with_prev_snp = $clustered_with_next_snp;
                    $clustered_with_next_snp = 0;
                    next;
                }
                # not clustered SNP
                print $out_wo_cluster_snp "$line\n";
                if ($sample_snp_different) {
                    print $out_wo_cluster_snp_sample_only "$line\n";
                }
            }

            # remove masked region
            if ($do_masking) {
                if(is_in_mask_region($chrom, $snp_hash_ref)){
                    print $out_masked_snp "$line\n";
                    next;
                }
                print $out_w_mask "$line\n";
                if ($sample_snp_different) {
                    print $out_w_mask_sample_only "$line\n";
                }
            }
            print $out_final_snp "$line\n";
        }
    }
    close $out_no_mask_sample_only;
    close $out_final_snp;
    if ($do_masking) {
        close $out_w_mask;
        close $out_w_mask_sample_only;
        close $out_masked_snp;
    }
    if ($remove_cluster_snp){
        close $out_wo_cluster_snp;
        close $out_wo_cluster_snp_sample_only;
        close $out_removed_cluster;
    }
}


sub get_rep {
    my ($in_file) = shift;
    open(my $in, '<', $in_file) or die "cannot open < '$in_file': $!\n";
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;   # 末尾の改行を削除
        next if ($. == 1);          # 1行目はスキップ
        next if ($line eq "");
        
        # タブで分割
        my ($chrom, $start, $end) = (split /\t/, $line) [0..2];
        
        $Mask{$chrom} = [] unless exists $Mask{$chrom};
        push @{$Mask{$chrom}}, {start => $start, end => $end};
    }
    close $in;

    # acsending order by start position
    foreach my $chrom (keys %Mask) {
        @{$Mask{$chrom}} = sort { $a->{start} <=> $b->{start} } @{$Mask{$chrom}};
    }
}

# 縦横の並べ替え
sub make_transpose {
    my ($in_file, $out_file) = @_;
    open(my $in, '<', $in_file) or die "cannot open < '$in_file': $!\n";
    open(my $out, '>', $out_file) or die "cannot open > '$out_file': $!\n";
    
    # ファイルデータを二次元配列にして@transに移す
    my $maxrow = 0;
    my @trans;
    
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\x0D?\x0A?$//;
        $line =~ tr|\t|,|;
        next if ($line eq "");
        
        my @row = split(/,/, $line);
        # chromとposition要素を削除
        splice (@row, 0, 2);
        push(@trans, \@row);
        $maxrow = @row if ($maxrow < @row);
    }
    close $in;
    
    # 縦横変換したものをファイルに出力
    for (my $r = 0; $r < $maxrow; $r++) {
        print $out ">";
        
        my $i = 0;
        foreach my $rowref (@trans) {
            print $out "\n" if ($i == 1);
            print $out $rowref->[$r] if ($rowref->[$r]);
            $i++;
        }
        print $out "\n";
    }
    close $out;
}
