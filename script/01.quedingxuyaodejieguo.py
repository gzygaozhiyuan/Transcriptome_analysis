import sys
import argparse
import numpy as np
import pandas as pd
import os

arg = argparse.ArgumentParser()
arg.add_argument('-GL', '--geneloc', required=True,type=str,help='Gene location file[Required]')
arg.add_argument('-IP', '--inputpath', required=True,type=str,help='eGWAS result path[Required]')
arg.add_argument('-OP', '--outputpath', required=True,type=str,help='Output path[Required]')
arg.add_argument('-D', '--distance', required=True,type=int,help='Cis distance[Required]')

args = arg.parse_args()

geneloc = args.geneloc
inputpath = args.inputpath
outputpath = args.outputpath
distance = args.distance

if inputpath[-1] != '/':
    inputpath += '/'
if outputpath[-1] != '/':
    outputpath += '/'

f1 = open(geneloc,'r')

for i in f1:
    if i.startswith('geneid'):
        pass
    else:
        ls = i.strip().split('\t')
        gene_name = ls[0]
        chr = int(ls[1][3:])
        start = int(ls[2])
        end = int(ls[3])
        gene_range = [start-distance, end+distance]
        # 判断文件是否存在
        if os.path.exists(inputpath + gene_name + '.txt_ck_result.assoc.txt') == False:
            print('Warning: ' + gene_name + ' is not exists!')
            pass
        else:
            # 读取对应的eGWAS结果
            egwas = pd.read_csv(inputpath + gene_name + '.txt_ck_result.assoc.txt',
                                sep='\t', header=0)
            # egwas添加一列type，值根据ps列进行判断，如果ps在gene_range范围内，则为cis，否则为trans
            egwas['type'] = egwas.apply(lambda x: 'cis' if gene_range[0] <= x['ps'] <= gene_range[1] and x['chr'] == chr else 'trans', axis=1)
            # 将p值小于0.05的结果保存
            egwas = egwas[egwas['p_wald'] < 1/500000]
            # 保存结果
            # egwas.to_csv(outputpath + gene_name + '_sa_result.assoc.txt', sep='\t', index=False, header=True)
            # 保存cis结果
            egwas[egwas['type']=='cis'].to_csv(outputpath + gene_name + '_sa_result.cis.assoc.txt', sep='\t', index=False, header=True)
            # 保存trans结果
            egwas[egwas['type']=='trans'].to_csv(outputpath + gene_name + '_sa_result.trans.assoc.txt', sep='\t', index=False, header=True)
            # 删除inputpath + gene_name + '.txt_ck_result.assoc.txt'
            # os.remove(inputpath + gene_name + '.txt_ck_result.assoc.txt')
f1.close()

