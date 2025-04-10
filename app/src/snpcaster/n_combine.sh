#!/usr/bin
# ドラフトゲノム等の配列間を100個のNでつなぎ、1本の配列にする。

#_usage
function _usage(){
	cat <<__EOF__

  Usage: $(basename "$0") [INPUT]
  INPUT input fasta file  (Required!)
__EOF__
}	

# 引数代入
INPUT="$1"
if [ -z $INPUT ]; then
  	_usage
  	exit 0
fi

DIR0=$BAC1/snpcaster
source activate perl-bioperl-core
perl $DIR0/N_combine.pl $INPUT
conda deactivate
