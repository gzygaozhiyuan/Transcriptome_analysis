import argparse
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser()

parser.add_argument('--inputpath', '-ip', type=str, required=True, help='Input path[Required]')
parser.add_argument('--me', '-m', type=int, default=578746, help='Me[578746]')
parser.add_argument('--threshold', '-t', type=float, default=0.01, help='Threshold[Default: 0.01]')
parser.add_argument('--output', '-o', type=str, default='output.txt', help='Output file name[Default: output.txt]')
# parser.add_argument('--LD', '-l', type=int, default=10000, help='LD[Default: 10,000]')

args = parser.parse_args()

inputpath = args.inputpath
me = args.me
threshold = args.threshold
output = args.output
# LD = args.LD

if inputpath[-1] != '/':
    inputpath += '/'

# 获取inputpath下的所有文件名
files_ls = os.listdir(inputpath)

# 创建一个空的dataframe，用于存储trans-eQTL的结果，列名为：snp，chr，ps，p_wald
trans_eQTL = pd.DataFrame(columns=['rs', 'chr', 'ps', 'p_wald','gene','beta','se','af','n_obs'])

# 读取文件
for file in files_ls:
    if 'assoc.txt' in file:
        egwas = pd.read_csv(inputpath + file, sep='\t', header=0)
        # 将egwas中p_wald列值小于(threshold/me)的行提取出来
        egwas_need = egwas[egwas['p_wald'] < (threshold / me)]
        gene = file[:14]
        snp_num = egwas_need.shape[0]
        egwas_need['gene'] = gene
        if snp_num != 0:
            # 将egwas_need中rs列填入trans_eQTL的snp列，将chr列填入trans_eQTL的chr列，将ps列填入trans_eQTL的ps列，将p_wald列填入trans_eQTL的p_wald列
            trans_eQTL = pd.concat([trans_eQTL, egwas_need[['rs', 'chr', 'ps', 'p_wald','gene','beta','se','af','n_obs']]], ignore_index=True)
        else:
            print(gene + ' has no trans-eQTL')

# 保存结果为txt文件
trans_eQTL.to_csv(output+'_all.txt', sep='\t', index=False, header=True)