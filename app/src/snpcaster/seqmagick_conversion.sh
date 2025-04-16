#!/usr/bin
# final_snp.fas 作成後実行する。 Ref無しバージョンを作成し、 nexus変換も作成する
# final_snp.fas が必要。別ファイルを使用する場合は引数で指定できる
# final_snp.fas --> final_snp_woRef.fas, final_snp.nex, final_snp_woRef.nex

# 引数代入
if [[ -n "$1" ]]; then
  FILE_FULL=$1
else
  FILE_FULL=final_snp.fas
  [ ! -e ${FILE_FULL} ] && echo "$FILE_FULL not found." && exit 1
fi

FILE=${FILE_FULL%.*}
#echo $FILE.fas

# final_snp.fasを別名コピー、 先頭2行を削除
cp ${FILE}.fas ${FILE}_woRef.fas
sed -i -e '1,2d' ${FILE}_woRef.fas
# nexus conversion
source activate seqmagick
seqmagick convert --output-format nexus --alphabet dna ${FILE}.fas ${FILE}.nex
seqmagick convert --output-format nexus --alphabet dna ${FILE}_woRef.fas ${FILE}_woRef.nex
conda deactivate
