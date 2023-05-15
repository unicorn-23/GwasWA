#!/usr/bin/env python
import sys
import os
import subprocess
from myparsing import *
from myutils import *

# 分为几个部分 
# 1、动物全基因组测序（WGS）基础流程：包括fastq到vcf
# 2、GWAS流程和数据挖掘、可视化
# 3、GWAS后续工作：预测变异影响enformer、VEP等
# 首先判断使用哪个模块
# 可以单独使用某个功能，比如
# 1、单独做wgs部分，输入fastq得到vcf及其他中间文件
# 2、单独做gwas，输入vcf和表型，得到关联分析结果
# 3、单独输入变异，得到变异预测结果
# 也可以从头开始的pipeline，比如
# 输入fastq和表型文件，得到gwas结果和变异预测结果
# 输入vcf和表型文件，得到gwas结果和变异预测结果




argvals = None

def main(argvals=argvals):
    P = Parser(argvals)
    version = set_version_number()
    print(f"\ngwaswa v{version} \n")
    # 输入fastq
    args = sys.argv[1:]
    args.insert(0, 'python3')
    args.insert(1, os.path.join(os.path.dirname(__file__), 'mystarter.py'))



    process = subprocess.run(args)
    return process.returncode


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt as e:
        print("\n interrupted")
    finally:
        # 删除临时目录
        # shutil.rmtree("_vcf2gwas_temp", ignore_errors=True)
        pass