import numpy as np
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='manual to this script')

parser.add_argument('--inputpath','-ip',type=str,required=True,help='input path[required]')
parser.add_argument('--snpbin','-sb',type=str,required=True,help='snp to bin file[required]')
parser.add_argument('--me', '-m', type=int, default=578746, help='Me[578746]')
parser.add_argument('--threshold', '-t', type=float, default=1, help='Threshold[Default: 1]')
parser.add_argument('--output', '-o', type=str, default='output.txt', help='Output file name[Default: output.txt]')
# parser.add_argument('--LD', '-l', type=int, default=10000, help='LD[Default: 10,000]')

args = parser.parse_args()

inputpath = args.inputpath
snpbin = args.snpbin
me = args.me
threshold = args.threshold
output = args.output
# LD = args.LD

if inputpath[-1] != '/':
    inputpath += '/'

# 获取inputpath下的所有文件名
files_ls = os.listdir(inputpath)

# 读取snpbin文件
trans_eQTL_unique_sorted = pd.read_csv(snpbin, sep='\t', header=0)

data_dict = {}
for i in range(len(trans_eQTL_unique_sorted)):
    data_dict[trans_eQTL_unique_sorted['rs'][i]] = trans_eQTL_unique_sorted['bin'][i]

# 创建一个dataframe，用于存储result
result = pd.DataFrame(columns=['gene', 'bin', 'lead_snp','snp_num', 'p_value', 'se', 'beta', 'af', 'n_obs'])

# 读取egwas结果
for file in files_ls:
    if 'assoc.txt' in file:
        egwas1 = pd.read_csv(inputpath + file, sep='\t', header=0)
        egwas1_need = egwas1[egwas1['p_wald'] < (threshold / me)]
        if egwas1_need.shape[0] > 0:
            egwas1_need['bin'] = egwas1_need['rs'].map(data_dict)
            egwas1_need = egwas1_need.groupby('bin').filter(lambda x: len(x) >= 3)
            gene1 = file[:14]
            if egwas1_need.shape[0] > 0:
                for bin_type in egwas1_need['bin'].unique():
                    egwas_bin = egwas1_need[egwas1_need['bin'] == bin_type]
                    snp_num = len(egwas_bin)
                    p_value = egwas_bin['p_wald'].min()
                    se = egwas_bin[egwas_bin['p_wald'] == p_value]['se'].values[0]
                    beta = egwas_bin[egwas_bin['p_wald'] == p_value]['beta'].values[0]
                    af = egwas_bin[egwas_bin['p_wald'] == p_value]['af'].values[0]
                    n_obs = egwas_bin[egwas_bin['p_wald'] == p_value]['n_obs'].values[0]
                    lead_snp = egwas_bin[egwas_bin['p_wald'] == p_value]['rs'].values[0]
                    result = result.append({'gene': gene1, 'bin': bin_type, 'snp_num': snp_num, 'p_value': p_value, 'se': se, 'beta': beta, 'af': af, 'n_obs': n_obs,'lead_snp': lead_snp}, ignore_index=True)
            else:
                print('No bin in ' + gene1)
        else:
            print('No significant in ' + gene1)
result['pve'] = (result['beta']**2 * result['af'] * (1 - result['af']))/(2 * result['se']**2*result['n_obs']*result['af'] * (1 - result['af'])+result['beta']**2 * result['af'] * (1 - result['af']))

# 保存结果为txt文件
result.to_csv(output, sep='\t', index=False, header=True)