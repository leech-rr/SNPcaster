# calculate core_genome_size
# 入力ファイル：args[1]
# 出力ファイル：core_genome.txt（合計の数値）


import pandas as pd
import sys

#for arg in sys.argv:
#  print(arg)

args=sys.argv
input_file_name=args[1] 
df=pd.read_table(input_file_name,sep="[\t]",engine='python')
df['D']=df['End']-df['Start']+1
total_value=df['D'].sum()
f=open('core_genome.txt','w')
f.write(str(total_value)+'\n')
f.close
