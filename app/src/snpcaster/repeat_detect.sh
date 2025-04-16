#!/bin/bash
# 参照ファイルから、tandems、 seq.repeats、 seq_seq.coordsを作成する
# $1:参照配列のファイル名

CMDNAME=`basename $0`

# parameter check
if [ "$#" -lt 1 ]; then
  echo "Usage: bash $CMDNAME ref.fas"
  echo '  $1  Reference Sequence File'
  exit 1
fi

# 参照配列FASTAファイル
REF=$1

# 参照ファイル名（拡張子なし）
filename=${REF%.*}

# FASTAファイルの1行目を書き換える
sed -i '1d' $REF
sed -i "1i >${filename}" $REF

# mummer4-4.0.0rc1
source activate mummer4
exact-tandems $REF 50 > tandems
repeat-match -n 50 $REF > seq.repeats

# large inexact repeats
nucmer --maxmatch --nosimplify --prefix=seq_seq $REF $REF
show-coords -r seq_seq.delta > seq_seq.coords
conda deactivate

# 不要ファイルの削除
rm seq_seq.delta

# tandems、seq.repeats、seq_seq.coordsのファイルを元に、repeats.tsvファイルを作成
source activate perl-bioperl-core
perl ${BAC1}/snpcaster/make_repeat_file.pl $filename
conda deactivate
