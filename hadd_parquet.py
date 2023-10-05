# -*- coding: utf-8 -*-
'''
merge parquet file with --hadd
command: $ python hadd_parquet.py --hadd a.parquet b.parquet c.parquet
'''
import argparse
import awkward as ak
import pandas as pd
# 创建 ArgumentParser 对象
parser = argparse.ArgumentParser()

# 定义列表参数
parser.add_argument("--hadd", nargs="+", type=str, help="一个或多个字符串")

# 解析命令行参数
args = parser.parse_args()

# 访问列表参数
print(args.hadd)

# 读取parquet
merge_parquet = []
for i in args.hadd:
    merge_parquet.append(ak.from_parquet(i))

# 将合并后的数据转换为 pandas DataFrame
result = ak.concatenate(merge_parquet)

# 保存合并后的数据到新的 .parquet 文件
ak.to_parquet(result, "merged_file.parquet")
