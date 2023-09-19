import argparse
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser()

parser.add_argument('--inputpath', '-ip', type=str, required=True, help='Input path[Required]')
parser.add_argument('--me', '-m', type=int, default=578746, help='Me[578746]')
parser.add_argument('--threshold', '-t', type=float, default=1, help='Threshold[Default: 1]')
parser.add_argument('--output', '-o', type=str, default='output.txt', help='Output file name[Default: output.txt]')

args = parser.parse_args()

inputpath = args.inputpath
me = args.me
threshold = args.threshold
output = args.output

if inputpath[-1] != '/':
    inputpath += '/'

# 获取inputpath下的所有文件名
files_ls = os.listdir(inputpath)

# 创建一个空的dataframe，用于存储cis-eQTL的结果，列名为：gene, lead_snp, snp_num, min_p, snp_ls
cis_eQTL = pd.DataFrame(columns=['gene', 'lead_snp', 'snp_num', 'min_p', 'snp_ls','lead_snp_beta','lead_snp_se','maf','n_obs'])

# 读取文件
for file in files_ls:
    if 'assoc.txt' in file:
        egwas = pd.read_csv(inputpath + file, sep='\t', header=0)
        gene = file[:14]
        # 将egwas中p_wald列值小于(threshold/me)的行提取出来
        egwas_need = egwas[egwas['p_wald'] < (threshold / me)]
        snp_num = egwas_need.shape[0]
        if snp_num != 0:
            egwas_min_p = egwas_need['p_wald'].min()
            snp_ls = '|'.join(egwas_need['rs'].tolist())
            lead_snp = egwas_need[egwas_need['p_wald'] == egwas_min_p]['rs'].tolist()[0]
            lead_snp_beta = egwas_need[egwas_need['p_wald'] == egwas_min_p]['beta'].tolist()[0]
            lead_snp_se = egwas_need[egwas_need['p_wald'] == egwas_min_p]['se'].tolist()[0]
            maf = egwas_need[egwas_need['p_wald'] == egwas_min_p]['af'].tolist()[0]
            n_obs = egwas_need[egwas_need['p_wald'] == egwas_min_p]['n_obs'].tolist()[0]
            cis_eQTL = cis_eQTL.append({'gene': gene, 'lead_snp': lead_snp, 'snp_num': snp_num, 'min_p': egwas_min_p, 'snp_ls': snp_ls,'lead_snp_beta':lead_snp_beta,'lead_snp_se':lead_snp_se,'maf':maf,'n_obs':n_obs}, ignore_index=True)
        else:
            print(gene + ' has no cis-eQTL')
# 保存结果为txt文件
cis_eQTL.to_csv(output, sep='\t', index=False, header=True)

