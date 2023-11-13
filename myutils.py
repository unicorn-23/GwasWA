#!/usr/bin/env python
import warnings
warnings.filterwarnings('ignore')
import shutil
from smtplib import LMTP
import time
import sys
import os
import stat
import subprocess
import threading
import re
from enformer.enformer_predict import *
from bs4 import BeautifulSoup
import html
import numpy as np
import json
import pandas as pd
import vcf
script_dir = os.path.dirname(os.path.abspath(__file__))

def cat_cmd(cmd_list):
    cmd_str = cmd_list[0]
    for i in range(1, len(cmd_list)):
        cmd_str = cmd_str + " && " + cmd_list[i]
    return cmd_str


def run_cmd(command, part):
    cmd_str = cat_cmd(command)
    # print("command########",cmd_str)
    p = subprocess.Popen(cmd_str, shell=True,
                         stdout=subprocess.DEVNULL, encoding="utf-8")
    p.wait()
    if p.returncode == 0:
        print(f'{part} succeeded')
        return True
    else:
        print(f'{part} failed')
        return False


def check_endswith(filename, str_list):
    flag = False
    for i in str_list:
        if filename.lower().endswith(i):
            flag = True
            break
    return flag


def get_fastq_file_path(fastq_list_path):
    tmp_single, tmp_paired = distinguish_single_paired(fastq_list_path)
    single_fastq_path_list = []
    paired_fastq_path_list = []

    for fastq in tmp_single:
        fastq_path = os.path.join(fastq_list_path, fastq)
        single_fastq_path_list.append(fastq_path)
    for paired_fastq in tmp_paired:
        tmp = []
        for fastq in paired_fastq:
            fastq_path = os.path.join(fastq_list_path, fastq)
            tmp.append(fastq_path)
        paired_fastq_path_list.append(tmp)

    return single_fastq_path_list, paired_fastq_path_list


def distinguish_single_paired(fastq_list_path):
    file_list = os.listdir(fastq_list_path)
    fastq_list_dir = []
    # 将文件夹下所有fq格式放入list
    for file_ in file_list:
        if check_endswith(file_, ['.fastq', ".fq", ".fastq.gz", ".fq.gz"]):
            fastq_list_dir.append(file_)
    # 判断是paired还是single
    single_fastq_list = []
    paired_fastq_list = []
    while len(fastq_list_dir) > 0:
        if (len(fastq_list_dir) == 1):
            single_fastq_list.append(fastq_list_dir[0])
            fastq_list_dir.pop(0)
            break

        f1 = fastq_list_dir[0]
        f1_id = re.match("[^_.]*", f1)
        f1_tail = re.search("_\d.", f1)
        if f1_id:
            tmp1 = None
            tmp2 = None
            tmp3 = True
            for i in range(1, len(fastq_list_dir)):
                f2 = fastq_list_dir[i]
                f2_id = re.search("[^_.]*", f2)
                f2_tail = re.search("_\d.", f2)
                if f2_id:
                    if f1_id.group() == f2_id.group():
                        if f1_tail:
                            if f2_tail:
                                paired_fastq_list.append([f1, f2])
                                tmp1 = fastq_list_dir[i]
                            else:
                                tmp2 = fastq_list_dir[i]
                        else:
                            tmp3 = False
                            break
            if tmp1 or tmp2:
                if tmp1:
                    fastq_list_dir.pop(fastq_list_dir.index(tmp1))
                if tmp2:
                    fastq_list_dir.pop(fastq_list_dir.index(tmp2))
            else:
                if tmp3:
                    single_fastq_list.append(fastq_list_dir[0])
        fastq_list_dir.pop(0)

    return single_fastq_list, paired_fastq_list


#################### Setting variables ####################
if sys.platform == 'darwin':
    sratoolkit_path = f"{script_dir}/sratoolkit/mac"
if sys.platform == 'linux':
    sratoolkit_path = f"{script_dir}/sratoolkit/linux"


def sratool_add_chmod():
    sratoolkit_list = os.listdir(sratoolkit_path)
    for tool in sratoolkit_list:
        tool_path = os.path.join(sratoolkit_path, tool)
        file_permission = os.stat(tool_path).st_mode
        os.chmod(tool_path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

# 设置权限需要
# prefetch工具


def get_prefetch_tool_path():
    return os.path.join(sratoolkit_path, "prefetch")

# vdb-validate工具


def get_vdb_validate_tool_path():
    return os.path.join(sratoolkit_path, "vdb-validate")


# fasterq-dump工具
def get_fasterq_dump_tool_path():
    return os.path.join(sratoolkit_path, "fasterq-dump")


# user-repository路径
def get_output_dir(P):
    out_dir = P.set_out_dir()
    return os.path.join(out_dir, "gwaswaOutput")

# 输出目录下的part目录


def get_wgs_output_path(P, part_name):
    output_dir = get_output_dir(P)
    wgs_output_path = os.path.join(output_dir, "wgs", part_name)
    os.makedirs(wgs_output_path, exist_ok=True)
    return wgs_output_path

# 输出目录下的part目录


def get_gwas_output_path(P, part_name):
    output_dir = get_output_dir(P)
    gwas_output_path = os.path.join(output_dir, "gwas", part_name)
    os.makedirs(gwas_output_path, exist_ok=True)
    return gwas_output_path

# 输出目录下的part目录


def get_assessment_output_path(P, part_name):
    output_dir = get_output_dir(P)
    assessment_output_path = os.path.join(output_dir, "assessment", part_name)
    os.makedirs(assessment_output_path, exist_ok=True)
    return assessment_output_path


#################### 1、  SRA数据下载 ####################
def do_download_sra(P, download_list):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []
    for srr in download_list:
        t = threading.Thread(target=download_content, args=(P, srr))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def download_content(P, srr):
    sra_output_path = get_wgs_output_path(P, "sra")
    prefetch_tool_path = get_prefetch_tool_path()
    print("downloading {} ...".format(srr))
    file_name = srr+".sra"
    sra_save_path = os.path.join(sra_output_path, file_name)
    prefetch_cmd = [
        f"{prefetch_tool_path} {srr} -o {sra_save_path}"
    ]
    wget_cmd = [
        f"wget https://sra-pub-run-odp.s3.amazonaws.com/sra/{srr}/{srr} --no-check-certificate -O {sra_save_path}"
    ]
    flag = run_cmd(prefetch_cmd, f"Download the {srr} using prefetch")
    if not flag:
        flag = run_cmd(wget_cmd, f"Download the {srr} using wget")


def check_sra(P):
    sra_output_path = get_wgs_output_path(P, "sra")
    file_list = os.listdir(sra_output_path)
    sra_list = []
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []
    for file_ in file_list:
        t = threading.Thread(target=check_content, args=(P, file_))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def check_content(P, file_):
    sra_output_path = get_wgs_output_path(P, "sra")
    validate_tool_path = get_vdb_validate_tool_path()
    if file_.endswith(".sra"):
        sra_file_path = os.path.join(sra_output_path, file_)
        print(f"check {file_}")
        check_cmd = [
            f"{validate_tool_path} {sra_file_path}"
        ]
        if not run_cmd(check_cmd, f"{file_} check"):
            err_log_txt = os.path.join(sra_output_path, "err_sra_log.txt")
            f = open(err_log_txt, 'a')
            f.write(file_+"\n")
            f.close()


def report_err_sra(P):
    sra_output_path = get_wgs_output_path(P, "sra")
    err_log_txt = os.path.join(sra_output_path, "err_sra_log.txt")
    if os.path.exists(err_log_txt):
        with open(err_log_txt, 'r') as f:
            for line in f:
                print(line.replace("\n", "  failed the inspection"))
        f.close()

#################### 2、  sra转fastq ####################


def get_sra_file_path(sra_list_path):
    sra_list = os.listdir(sra_list_path)
    sra_file_path_list = []
    for sra in sra_list:
        if sra.endswith(".sra"):
            sra_path = os.path.join(sra_list_path, sra)
            sra_file_path_list.append(sra_path)
    return sra_file_path_list


def do_sra_to_fastq(P, sra_file_path_list):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []
    for sra in sra_file_path_list:
        t = threading.Thread(target=fastq_content, args=(P, sra))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def fastq_content(P, sra):
    fasterq_dump_tool_path = get_fasterq_dump_tool_path()
    fastq_output_path = get_wgs_output_path(P, "raw")
    print(f"Convert {sra} to fastq start")
    nThrds = P.get_nThrds()
    cmd = f"{fasterq_dump_tool_path} --split-3"
    if nThrds:
        cmd += f" -e {nThrds}"
    fastq_cmd = [
        f"{cmd} {sra} -O {fastq_output_path}",
    ]
    if run_cmd(fastq_cmd, f"Convert {sra} to fastq"):
        if P.get_nosave():
            if os.path.exists(sra):
                os.remove(sra)


def zip_fastq(P):
    raw_path = get_wgs_output_path(P, "raw")
    tmp = []
    single_fastq_path_list, paired_fastq_path_list = get_fastq_file_path(
        raw_path)
    tmp.extend(single_fastq_path_list)
    count = 0
    thread_list = []
    nThrds = P.get_nThrds()
    for t in paired_fastq_path_list:
        tmp.extend(t)
    for fastq in tmp:
        if not fastq.endswith(".gz"):
            cmd = "bgzip"
            if nThrds:
                cmd += f" -@ {nThrds}"
            zip_cmd = [
                f"{cmd} -f {fastq}",
            ]
            print(f"zip .fastq file {fastq} start")

            t = threading.Thread(target=run_cmd, args=(
                zip_cmd, f"zip .fastq file {fastq}"))
            thread_list.append(t)
            t.start()
            if count < nThrds:
                count += 1
            if count == nThrds:
                for t in thread_list:
                    t.join()
                count = 0
                thread_list = []

#################### 3、  fastq质控 ####################


def do_qc(P, single_fastq_file_path_list, paired_fastq_file_path_list):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []
    for fastq_file in single_fastq_file_path_list:
        t = threading.Thread(target=qc_single_content, args=(P, fastq_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for fastq_file in paired_fastq_file_path_list:
        t = threading.Thread(target=qc_paired_content, args=(P, fastq_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def qc_single_content(P, fastq_file):
    fastq_output_path = get_wgs_output_path(P, "clean")
    phred = P.get_phred()
    quality = P.get_quality()
    length = P.get_length()
    stringency = P.get_stringency()
    error = P.get_error()
    nThrds = P.get_nThrds()
    cmd = "trim_galore --gzip"
    if nThrds:
        cmd += f" --cores {nThrds}"
    qc_cmd = [
        f"{cmd} --{phred} -q {quality} --length {length} --stringency {stringency} -e {error} {fastq_file} -o {fastq_output_path}"
    ]

    print(f"Quality control of {fastq_file} start")
    success = run_cmd(qc_cmd, f"Quality control of {fastq_file}")
    if success:
        if P.get_nosave():
            if os.path.exists(fastq_file):
                os.remove(fastq_file)


def qc_paired_content(P, fastq_file):
    fastq_output_path = get_wgs_output_path(P, "clean")
    phred = P.get_phred()
    quality = P.get_quality()
    length = P.get_length()
    stringency = P.get_stringency()
    error = P.get_error()
    fastq_file1, fastq_file2 = fastq_file
    nThrds = P.get_nThrds()
    cmd = "trim_galore --gzip"
    if nThrds:
        cmd += f" --cores {nThrds}"
    qc_cmd = [
        f"{cmd} --{phred} -q {quality} --length {length} --stringency {stringency} -e {error} --paired {fastq_file1} {fastq_file2} -o {fastq_output_path}"
    ]
    print(f"Quality control of {fastq_file}")
    if run_cmd(qc_cmd, f"Quality control of {fastq_file}"):
        if P.get_nosave():
            if os.path.exists(fastq_file1):
                os.remove(fastq_file1)
            if os.path.exists(fastq_file2):
                os.remove(fastq_file2)
#################### 4、  multiqc质量评估 ####################
def get_qc_file_path(clean_output_path):
    fastq_list = os.listdir(clean_output_path)
    fastq_file_path_list = []
    end_list = [".fq", ".fastq", ".fq.gz", "fastq.gz"]
    for fastq_file in fastq_list:
        if check_endswith(fastq_file, end_list):
            fastq_path = os.path.join(clean_output_path, fastq_file)
            fastq_file_path_list.append(fastq_path)
    return fastq_file_path_list


def do_fastqc(P, qc_file_path_list):
    fastqc_part = os.path.join("qualityEvaluation", "fastqc")
    fastqc_output_path = get_wgs_output_path(P, fastqc_part)
    nThrds = P.get_nThrds()
    for fastq_file in qc_file_path_list:
        cmd = "fastqc"
        if nThrds:
            cmd += f" -t {nThrds}"
        f_cmd = [
            f"{cmd} -o {fastqc_output_path} {fastq_file}"
        ]
        print(f"Use fastqc for quality evaluation of {fastq_file} start")
        run_cmd(f_cmd, f"Use fastqc for quality evaluation of {fastq_file}")


def do_multiqc(P):
    fastqc_part = os.path.join("qualityEvaluation", "fastqc")
    multiqc_part = os.path.join("qualityEvaluation", "multiqc")

    fastqc_output_path = get_wgs_output_path(P, fastqc_part)
    multiqc_output_path = get_wgs_output_path(P, multiqc_part)

    cmd = [
        f"multiqc {fastqc_output_path} -o {multiqc_output_path}"
    ]
    print("Use multiqc to evaluate the quality of fastqc's results start")
    run_cmd(cmd, f"Use multiqc to evaluate the quality of fastqc's results")


#################### 5、  建立基因组索引 ####################
def download_ref(P):
    accession = P.get_accession()
    taxon = P.get_taxon()
    refgenome = get_ref(P)
    if accession:
        ref_output_path = get_wgs_output_path(P, "ref")
        output_file = os.path.join(ref_output_path, f"{accession}.zip")
        cmd = [
            f"datasets download genome accession {accession} --filename {output_file}",
            f"unzip -j {output_file}  -d {ref_output_path}"
        ]
        print(f"Download reference genome {accession}  start")
        if run_cmd(cmd, f"Download reference genome {accession}"):
            file_list = os.listdir(ref_output_path)
            for file_ in file_list:
                if file_.endswith(".fna") or file_.endswith(".fa") or file_.endswith(".fasta"):
                    os.rename(os.path.join(ref_output_path, file_),
                              os.path.join(ref_output_path, "ref.fa"))
            if os.path.exists(output_file):
                os.remove(output_file)
    elif taxon:
        ref_output_path = get_wgs_output_path(P, "ref")
        output_file = os.path.join(ref_output_path, f"{taxon}.zip")
        cmd = [
            f"datasets download genome taxon {taxon} --filename {output_file}",
            f"unzip -j {output_file}  -d {ref_output_path}"
        ]
        print(f"Download reference genome {taxon} start")
        if run_cmd(cmd, f"Download reference genome {taxon}"):
            file_list = os.listdir(ref_output_path)
            for file_ in file_list:
                if file_.endswith(".fna") or file_.endswith(".fa") or file_.endswith(".fasta"):
                    os.rename(os.path.join(ref_output_path, file_),
                              os.path.join(ref_output_path, "ref.fa"))

            if os.path.exists(output_file):
                os.remove(output_file)

    if do_index(P, refgenome):
        if create_gatk_dict(P, refgenome):
            create_samtool_index(P, refgenome)


def do_index(P, refgenome):
    indexalgorithm = P.get_indexalgorithm()
    refname = os.path.basename(refgenome)
    output_file = os.path.join(os.path.dirname(refgenome), refname)
    cmd = [
        f"bwa index -a {indexalgorithm} -p {output_file} {refgenome}"
    ]
    print(f"Create an index of the reference genome for tool BWA start")
    return run_cmd(cmd, f"Create an index of the reference genome for tool BWA")


def create_gatk_dict(P, refgenome):
    match = re.search("\.fas?t?a?\.?g?z?", os.path.basename(refgenome))
    if match:
        postfix = match.group()
    output_file = refgenome.replace(postfix, ".dict")
    if not os.path.exists(output_file):
        cmd = [
            f"gatk CreateSequenceDictionary -R {refgenome} -O {output_file}"
        ]
        print("Create an index of the reference genome for tool GATK start")
        return run_cmd(cmd, f"Create an index of the reference genome for tool GATK")


def create_samtool_index(P, refgenome):
    cmd = [
        f"samtools faidx {refgenome}"
    ]
    print("Create an index of the reference genome for tool Samtools start")
    return run_cmd(cmd, f"Create an index of the reference genome for tool Samtools")


def get_ref(P):
    refgenome = P.get_refgenome()
    if not refgenome:
        ref_output_path = get_wgs_output_path(P, "ref")
        refgenome = os.path.join(ref_output_path, "ref.fa")
        if not os.path.exists(refgenome):
            print("Local reference genome needs to be set or downloaded")
    return refgenome
#################### 6、  bwa比对参考基因组 ####################


def do_align_bwa(P, single_fastq_file_path_list, paired_fastq_file_path_list):
    alignalgorithm = P.get_alignalgorithm()
    if alignalgorithm == "backtrack":
        aln(P, single_fastq_file_path_list, paired_fastq_file_path_list)
    elif alignalgorithm == "mem":
        mem(P, single_fastq_file_path_list, paired_fastq_file_path_list)
    elif alignalgorithm == "bwasw":
        bwasw(P, single_fastq_file_path_list, paired_fastq_file_path_list)


def mem(P, single_fastq_file_path_list, paired_fastq_file_path_list):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []
    for fastq_file in single_fastq_file_path_list:
        t = threading.Thread(target=mem_single_content, args=(P, fastq_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for fastq_file in paired_fastq_file_path_list:
        t = threading.Thread(target=mem_paired_content, args=(P, fastq_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def mem_single_content(P, fastq_file):
    align_output_path = get_wgs_output_path(P, "align")
    refgenome = get_ref(P)
    base_name = os.path.basename(fastq_file)
    file_name = re.match("[^_.]*", base_name).group()
    output_sam_file = os.path.join(align_output_path, file_name)+".sam"
    output_bam_file = os.path.join(align_output_path, file_name)+".bam"
    mem_cmd = make_mem_cmd(P)
    cmd = [
        f"{mem_cmd} -R '@RG\\tID:{file_name}\\tPL:unknown\\tLB:{file_name}\\tSM:{file_name}' {refgenome} {fastq_file} > {output_sam_file}",
    ]
    print(f"Alignment of {fastq_file} using algorithm mem start")
    if run_cmd(cmd, f"Alignment of {fastq_file} using algorithm mem"):
        if P.get_nosave():
            if os.path.exists(fastq_file):
                os.remove(fastq_file)
        sam_to_bam(P, output_sam_file, output_bam_file)


def mem_paired_content(P, fastq_file):
    align_output_path = get_wgs_output_path(P, "align")
    refgenome = get_ref(P)
    fastq_file1, fastq_file2 = fastq_file
    base_name = os.path.basename(fastq_file1)
    file_name = re.match("[^_.]*", base_name).group()
    output_sam_file = os.path.join(align_output_path, file_name)+".sam"
    output_bam_file = os.path.join(align_output_path, file_name)+".bam"
    mem_cmd = make_mem_cmd(P)
    cmd = [
        f"{mem_cmd} -R '@RG\\tID:{file_name}\\tPL:unknown\\tLB:{file_name}\\tSM:{file_name}' {refgenome} {fastq_file1} {fastq_file2} > {output_sam_file}"
    ]
    print(f"Alignment of {fastq_file} using algorithm mem start")
    if run_cmd(cmd, f"Alignment of {fastq_file} using algorithm mem"):
        if P.get_nosave():
            if os.path.exists(fastq_file1):
                os.remove(fastq_file1)
            if os.path.exists(fastq_file2):
                os.remove(fastq_file2)
        sam_to_bam(P, output_sam_file, output_bam_file)


def make_mem_cmd(P):
    cmd = "bwa mem"
    minSeedLen = P.get_minSeedLen()
    bandWidth = P.get_bandWidth()
    zDropoff = P.get_zDropoff()
    seedSplitRatio = P.get_seedSplitRatio()
    maxOcc = P.get_maxOcc()
    matchScore = P.get_matchScore()
    mmPenalty = P.get_mmPenalty()
    gapOpenPen = P.get_gapOpenPen()
    gapExtPen = P.get_gapExtPen()
    clipPen = P.get_clipPen()
    unpairPen = P.get_unpairPen()
    nThrds = P.get_nThrds()
    if minSeedLen:
        cmd += f" -k {minSeedLen}"
    if bandWidth:
        cmd += f" -w {bandWidth}"
    if zDropoff:
        cmd += f" -d {zDropoff}"
    if seedSplitRatio:
        cmd += f" -r {seedSplitRatio}"
    if maxOcc:
        cmd += f" -c {maxOcc}"
    if matchScore:
        cmd += f" -A {matchScore}"
    if mmPenalty:
        cmd += f" -B {mmPenalty}"
    if gapOpenPen:
        cmd += f" -O {gapOpenPen}"
    if gapExtPen:
        cmd += f" -E {gapExtPen}"
    if clipPen:
        cmd += f" -L {clipPen}"
    if unpairPen:
        cmd += f" -U {unpairPen}"
    if nThrds:
        cmd += f" -t {nThrds}"
    return cmd


def aln(P, single_fastq_file_path_list, paired_fastq_file_path_list):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []

    for fastq_file in single_fastq_file_path_list:
        t = threading.Thread(target=aln_single_content, args=(P, fastq_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for fastq_file in paired_fastq_file_path_list:
        t = threading.Thread(target=aln_paired_content, args=(P, fastq_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def aln_single_content(P, fastq_file):
    align_output_path = get_wgs_output_path(P, "align")
    refgenome = get_ref(P)
    base_name = os.path.basename(fastq_file)
    file_name = re.match("[^_.]*", base_name).group()
    output_sai_file = os.path.join(align_output_path, file_name)+".sai"
    output_sam_file = os.path.join(align_output_path, file_name)+".sam"
    output_bam_file = os.path.join(align_output_path, file_name)+".bam"
    single_cmd = make_aln_cmd(P)
    samse_cmd = make_samse_cmd(P)
    cmd = [
        f"{single_cmd} {refgenome} {fastq_file} > {output_sai_file}",
        f"{samse_cmd} -r '@RG\\tID:{file_name}\\tPL:unknown\\tLB:{file_name}\\tSM:{file_name}' {refgenome}  {output_sai_file} {fastq_file} > {output_sam_file}"
    ]
    print(f"Alignment of {fastq_file} using algorithm backtrack")
    if run_cmd(cmd, f"Alignment of {fastq_file} using algorithm backtrack"):
        if os.path.exists(output_sai_file):
            os.remove(output_sai_file)
        if P.get_nosave():
            if os.path.exists(fastq_file):
                os.remove(fastq_file)
        sam_to_bam(P, output_sam_file, output_bam_file)


def aln_paired_content(P, fastq_file):
    align_output_path = get_wgs_output_path(P, "align")
    refgenome = get_ref(P)
    fastq_file1, fastq_file2 = fastq_file
    base_name = os.path.basename(fastq_file1)
    file_name = re.match("[^_.]*", base_name).group()
    output_sai_file1 = os.path.join(align_output_path, file_name)+"_1.sai"
    output_sai_file2 = os.path.join(align_output_path, file_name)+"_2.sai"
    output_sam_file = os.path.join(align_output_path, file_name)+".sam"
    output_bam_file = os.path.join(align_output_path, file_name)+".bam"
    paired_cmd = make_aln_cmd(P)
    samse_cmd = make_sampe_cmd(P)
    cmd = [
        f"{paired_cmd} {refgenome} {fastq_file1} > {output_sai_file1}",
        f"{paired_cmd} {refgenome} {fastq_file2} > {output_sai_file2}",
        f"{samse_cmd} -r '@RG\\tID:{file_name}\\tPL:unknown\\tLB:{file_name}\\tSM:{file_name}' {refgenome}  {output_sai_file1} {output_sai_file2} {fastq_file2} {fastq_file1} > {output_sam_file}"
    ]
    print(f"Alignment of {fastq_file} using algorithm backtrack")
    if run_cmd(cmd, f"Alignment of {fastq_file} using algorithm backtrack"):
        if os.path.exists(output_sai_file1):
            os.remove(output_sai_file1)
        if os.path.exists(output_sai_file2):
            os.remove(output_sai_file2)
        if P.get_nosave():
            if os.path.exists(fastq_file1):
                os.remove(fastq_file1)
            if os.path.exists(fastq_file2):
                os.remove(fastq_file2)
        sam_to_bam(P, output_sam_file, output_bam_file)


def make_aln_cmd(P):
    cmd = "bwa aln"
    maxDiff = P.get_maxDiff()
    maxGapO = P.get_maxGapO()
    maxGapE = P.get_maxGapE()
    nDelTail = P.get_nDelTail()
    nIndelEnd = P.get_nIndelEnd()
    maxSeedDiff = P.get_maxSeedDiff()
    seedLen = P.get_seedLen()
    misMsc = P.get_misMsc()
    gapOsc = P.get_gapOsc()
    gapEsc = P.get_gapEsc()
    trimQual = P.get_trimQual()
    nThrds = P.get_nThrds()
    if maxDiff:
        cmd += f" -n {maxDiff}"
    if maxGapO:
        cmd += f" -o {maxGapO}"
    if maxGapE:
        cmd += f" -e {maxGapE}"
    if nDelTail:
        cmd += f" -d {nDelTail}"
    if nIndelEnd:
        cmd += f" -i {nIndelEnd}"
    if maxSeedDiff:
        cmd += f" -k {maxSeedDiff}"
    if seedLen:
        cmd += f" -l {seedLen}"
    if misMsc:
        cmd += f" -M {misMsc}"
    if gapOsc:
        cmd += f" -O {gapOsc}"
    if gapEsc:
        cmd += f" -E {gapEsc}"
    if trimQual:
        cmd += f" -q {trimQual}"
    if nThrds:
        cmd += f" -t {nThrds}"
    return cmd


def make_samse_cmd(P):
    cmd = "bwa samse"
    maxHitPaired = P.get_maxHitPaired()
    if maxHitPaired:
        cmd += f" -n {maxHitPaired}"
    return cmd


def make_sampe_cmd(P):
    cmd = "bwa sampe"
    maxHitPaired = P.get_maxHitPaired()
    maxInsSize = P.get_maxInsSize()
    maxOcc = P.get_maxOcc()
    maxHitDis = P.get_maxHitDis()
    if maxInsSize:
        cmd += f" -a {maxInsSize}"
    if maxOcc:
        cmd += f" -o {maxOcc}"
    if maxHitPaired:
        cmd += f" -n {maxHitPaired}"
    if maxHitDis:
        cmd += f" -N {maxHitDis}"

    return cmd


def bwasw(P, single_fastq_file_path_list, paired_fastq_file_path_list):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []
    for fastq_file in single_fastq_file_path_list:
        t = threading.Thread(target=bwasw_single_content, args=(P, fastq_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for fastq_file in paired_fastq_file_path_list:
        t = threading.Thread(target=bwasw_paired_content, args=(P, fastq_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def bwasw_single_content(P, fastq_file):
    align_output_path = get_wgs_output_path(P, "align")
    refgenome = get_ref(P)
    base_name = os.path.basename(fastq_file)
    file_name = re.match("[^_.]*", base_name).group()
    output_sam_file = os.path.join(align_output_path, file_name)+".sam"
    output_bam_file = os.path.join(align_output_path, file_name)+".bam"
    bwasw_cmd = make_bwasw_cmd(P)
    cmd = [
        f"{bwasw_cmd} {refgenome} {fastq_file} > {output_sam_file}",
    ]
    print(f"Alignment of {fastq_file} using algorithm bwasw")
    if run_cmd(cmd, f"Alignment of {fastq_file} using algorithm bwasw"):
        if P.get_nosave():
            if os.path.exists(fastq_file):
                os.remove(fastq_file)
        sam_to_bam(P, output_sam_file, output_bam_file)


def bwasw_paired_content(P, fastq_file):
    align_output_path = get_wgs_output_path(P, "align")
    refgenome = get_ref(P)
    fastq_file1, fastq_file2 = fastq_file
    base_name = os.path.basename(fastq_file1)
    file_name = re.match("[^_.]*", base_name).group()
    output_sam_file = os.path.join(align_output_path, file_name)+".sam"
    output_bam_file = os.path.join(align_output_path, file_name)+".bam"
    bwasw_cmd = make_bwasw_cmd(P)
    cmd = [
        f"{bwasw_cmd} {refgenome} {fastq_file1} {fastq_file2} > {output_sam_file}"
    ]
    print(f"Alignment of {fastq_file} using algorithm bwasw")
    if run_cmd(cmd, f"Alignment of {fastq_file} using algorithm bwasw"):
        if P.get_nosave():
            if os.path.exists(fastq_file1):
                os.remove(fastq_file1)
            if os.path.exists(fastq_file2):
                os.remove(fastq_file2)
        sam_to_bam(P, output_sam_file, output_bam_file)


def make_bwasw_cmd(P):
    cmd = "bwa bwasw"
    matchScore = P.getmatchScore()
    mmPen = P.getmmPen()
    bandWidth = P.getbandWidth()
    thres = P.getthres()
    hspIntv = P.gethspIntv()
    zBest = P.getzBest()
    nHspRev = P.getnHspRev()
    thresCoef = P.getthresCoef()
    gapOpenPen = P.get_gapOpenPen()
    gapExtPen = P.get_gapExtPen()
    nThrds = P.get_nThrds()

    if matchScore:
        cmd += f" -a {matchScore}"
    if mmPen:
        cmd += f" -b {mmPen}"
    if gapOpenPen:
        cmd += f" -q {gapOpenPen}"
    if gapExtPen:
        cmd += f" -r {gapExtPen}"
    if bandWidth:
        cmd += f" -w {bandWidth}"
    if thres:
        cmd += f" -T {thres}"
    if hspIntv:
        cmd += f" -s {hspIntv}"
    if zBest:
        cmd += f" -z {zBest}"
    if nHspRev:
        cmd += f" -N {nHspRev}"
    if thresCoef:
        cmd += f" -c {thresCoef}"
    if nThrds:
        cmd += f" -t {nThrds}"
    return cmd


def sam_to_bam(P, output_sam_file, output_bam_file):
    nThrds = P.get_nThrds()
    cmd = "samtools view -S"
    if nThrds:
        cmd += f" -@ {nThrds}"

    stbcmd = [
        f"{cmd} -b {output_sam_file} > {output_bam_file}"
    ]
    print(f"Converting {output_sam_file} to .bam start")
    if run_cmd(stbcmd, f"Converting {output_sam_file} to .bam"):
        if os.path.exists(output_sam_file):
            os.remove(output_sam_file)

#################### 7、  处理bam文件 ####################


def sort_bam(P, bam_file_path):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []
    bam_list = os.listdir(bam_file_path)
    for bam_file in bam_list:
        t = threading.Thread(target=sort_content,
                             args=(P, bam_file, bam_file_path))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def sort_content(P, bam_file, bam_file_path):
    processed_output_path = get_wgs_output_path(P, "processed")
    input_file = os.path.join(bam_file_path, bam_file)
    output_file = os.path.join(processed_output_path, bam_file)
    nThrds = P.get_nThrds()
    cmd = "samtools sort -m 4G -O bam"
    if nThrds:
        cmd += f" -@ {nThrds}"

    sort_cmd = [
        f"{cmd} -o {output_file} {input_file}"
    ]

    print(f"Sequencing for {bam_file} start")
    if run_cmd(sort_cmd, f"Sequencing for {bam_file}"):
        if P.get_nosave():
            if os.path.exists(input_file):
                os.remove(input_file)


def mark_pcr_repeat_picard(P):
    if P.get_delPCR():
        count = 0
        nThrds = P.get_nThrds()
        thread_list = []

        processed_output_path = get_wgs_output_path(P, "processed")
        file_list = os.listdir(processed_output_path)
        bam_list = []
        for file_ in file_list:
            if not file_.endswith("_marked.bam"):
                bam_list.append(file_.replace(".bam", ""))

        for bam_file in bam_list:
            t = threading.Thread(target=del_pcr_content, args=(P, bam_file))
            thread_list.append(t)
            t.start()
            if count < nThrds:
                count += 1
            if count == nThrds:
                for t in thread_list:
                    t.join()
                count = 0
                thread_list = []
        for t in thread_list:
            t.join()


def del_pcr_content(P, bam_file):
    processed_output_path = get_wgs_output_path(P, "processed")
    input_file = os.path.join(processed_output_path, bam_file)+".bam"
    output_file = os.path.join(processed_output_path, bam_file)+"_marked.bam"
    markdup_metrics_file = os.path.join(
        processed_output_path, bam_file)+"_marked_metrics.txt"
    cmd = [
        f"gatk MarkDuplicates \
            -I {input_file} \
            -O {output_file} \
            -M {markdup_metrics_file}"
    ]
    print(f"Marking PCR Duplicates in {bam_file} start")
    if run_cmd(cmd, f"Marking PCR Duplicates in {bam_file}"):
        if P.get_nosave():
            if os.path.exists(input_file):
                os.remove(input_file)


def index_for_bam(P):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []

    processed_output_path = get_wgs_output_path(P, "processed")
    file_list = os.listdir(processed_output_path)
    bam_list = []
    for file_ in file_list:
        if P.get_delPCR():
            if file_.endswith("_marked.bam"):
                bam_list.append(file_)
        else:
            bam_list.append(file_)

    for bam_file in bam_list:
        t = threading.Thread(target=index_bam_content, args=(P, bam_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def index_bam_content(P, bam_file):
    processed_output_path = get_wgs_output_path(P, "processed")
    input_file = os.path.join(processed_output_path, bam_file)
    cmd = "samtools index"
    nThrds = P.get_nThrds()
    if nThrds:
        cmd += f" -@ {nThrds}"

    index_bam_cmd = [
        f"{cmd} {input_file}"
    ]

    print(f"Building an index for {bam_file} start")
    run_cmd(index_bam_cmd, f"Building an index for {bam_file}")

# BaseRecalibrator 碱基质量重校正BQSR这里计算出了所有需要进行重校正的read和特征值，然后把这些信息输出为一份校准表文件（sample_name.recal.table）
# def do_base_recal(P):
#     align_output_path = get_wgs_output_path(P, "align")
#     file_list = os.listdir(align_output_path)
#     sorted_marked_bam_list = []
#     for file_ in file_list:
#         if file_.endswith("_sorted_marked_realign.bam"):
#             sorted_marked_bam_list.append(
#                 file_.replace("_sorted_marked_realign.bam", ""))

#     for bam_file in sorted_marked_bam_list:
#         input_file = os.path.join(
#             align_output_path, bam_file)+"_sorted_marked_realign.bam"
#         output_file = os.path.join(align_output_path, bam_file)+"_recal.table"
#         cmd = [
#             f"gatk BaseRecalibrator  \
#                 -R {hg38_ref}\
#                 -I {input_file} \
#                 -known {hg38_snp} \
#                 -known {hg38_indel} \
#                 -o {output_file}"
#         ]

#         print(bam_file, "base_recalibrator start")
#         run_cmd(cmd, f"{bam_file} base_recalibrator")


# PrintReads，这一步利用第一步得到的校准表文件（sample_name.recal.table）重新调整原来BAM文件中的碱基质量值，
# 并使用这个新的质量值重新输出一份新的BAM文件。
# def do_printreads(P):
#     align_output_path = get_wgs_output_path(P, "align")
#     file_list = os.listdir(align_output_path)
#     sorted_marked_bam_list = []
#     for file_ in file_list:
#         if file_.endswith("_sorted_marked_realign.bam"):
#             sorted_marked_bam_list.append(
#                 file_.replace("_sorted_marked_realign.bam", ""))

#     for bam_file in sorted_marked_bam_list:
#         input_file = os.path.join(
#             align_output_path, bam_file)+"_sorted_marked_realign.bam"
#         recal_table = os.path.join(align_output_path, bam_file)+"_recal.table"
#         output_file = os.path.join(
#             align_output_path, bam_file)+"_sorted_marked_realign_bqsr.bam"
#         cmd = [
#             f"gatk PrintReads  \
#                 -R {hg38_ref}\
#                 -I {input_file} \
#                 --BQSR {recal_table} \
#                 -o {output_file}"
#         ]
#         print(bam_file, "applyBQSR start")
#         run_cmd(cmd, f"{bam_file} applyBQSR")


#################### 8、  变异检测 ####################
# def add_read_group(P):

# java -jar picard.jar AddOrReplaceReadGroups \
#        I=input.bam \
#        O=output.bam \
#        RGID=4 \
#        RGLB=library1 \
#        RGPL=illumina \
#        RGPU=unit1 \
#        SORT_ORDER=coordinate \
#        RGSM=sample1

def get_info_from_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        chrs = set()
        for line in f:
            line = line.strip()
            if line[0] == '>':
                match = re.search("chro?m?o?s?o?m?e? ?.?\d?", line)
                if match:
                    chrs.add(match.group())
                    # chrs.append(line)
    return chrs


def do_gvcf(P, bam_file_path):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []

    file_list = os.listdir(bam_file_path)
    sorted_marked_bam_list = []

    for file_ in file_list:
        if file_.endswith("_marked.bam"):
            sorted_marked_bam_list.append(file_)
    if len(sorted_marked_bam_list) == 0:
        for file_ in file_list:
            if file_.endswith(".bam"):
                sorted_marked_bam_list.append(file_)

    for bam_file in sorted_marked_bam_list:
        t = threading.Thread(target=do_gvcf_content,
                             args=(P, bam_file, bam_file_path))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()


def do_gvcf_content(P, bam_file, bam_file_path):
    refgenome = get_ref(P)
    vcf_output_path = get_wgs_output_path(P, "gvcf")
    input_file = os.path.join(bam_file_path, bam_file)
    if bam_file.endswith("_marked.bam"):
        bam = bam_file.replace("_marked.bam", "")
    else:
        bam = bam_file.replace(".bam", "")
    output_file = os.path.join(vcf_output_path, bam)+"_g.vcf"
    cmd = [
        f"gatk HaplotypeCaller \
            -R {refgenome} \
            -I {input_file} \
            --emit-ref-confidence GVCF \
            -O {output_file}"
    ]
    print(f"Variant detection for file {bam_file} start")

    if run_cmd(cmd, f"Variant detection for file {bam_file}"):
        if P.get_nosave():
            if os.path.exists(bam_file):
                os.remove(bam_file)

#################### 9、  jointgenotype ####################

def do_jointgenotype(P,gvcf_file_path):
    file_dict_by_group = do_split_by_chr_gvcf(P,gvcf_file_path)
    chr_file_list = do_chr_gvcf_merge(P,file_dict_by_group)
    done_chr_file_list = do_joint_genotype(P,chr_file_list)
    
    do_merge_chr_vcf(P,done_chr_file_list)

def do_split_by_chr_gvcf(P,gvcf_file_path):
    gvcf_bychr_output_path = get_wgs_output_path(P, "gvcf_chr")

    gvcf_list = []
    file_list = os.listdir(gvcf_file_path)
    for file_ in file_list:
        if file_.endswith("_g.vcf"):
            gvcf_list.append(file_.replace("_g.vcf", ""))

    file_dict_by_group={}
    for sample_name in gvcf_list:
        input_file = os.path.join(gvcf_file_path, sample_name)+"_g.vcf"
        vcf_file = pd.read_csv(input_file, sep='\t', comment='#', header=None)
        grouped = vcf_file.groupby(0)

        with open(input_file) as f:
            comments = [line.strip() for line in f if line.startswith('#')]

        for name, group in grouped:
            chr_file_name=os.path.join(gvcf_bychr_output_path,f"{sample_name}_chr_{name}_g.vcf")
            if name not in file_dict_by_group:
                file_dict_by_group[name]=[]
            file_dict_by_group[name].append(chr_file_name)
            with open(chr_file_name, 'w') as f:
                f.write('\n'.join(comments) + '\n')
                group.to_csv(f, sep='\t', index=False, header=None)
    # file_dict_by_group:{'chr1':[path1,path2]}
    print("Divide gvcf files by chromosome done")
    return file_dict_by_group

def do_chr_gvcf_merge(P,file_dict_by_group):
    vcf_output_path = get_wgs_output_path(P, "vcf")

    count = 0
    nThrds = P.get_nThrds()
    thread_list = []
    chr_file_list=[]
    for key in file_dict_by_group.keys():
        output_file = os.path.join(vcf_output_path, f"chr_{key}_g.vcf")
        chr_file_list.append(output_file)

        t = threading.Thread(target=chr_gvcf_merge_content,args=(P, output_file,file_dict_by_group[key]))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()

    return chr_file_list

def chr_gvcf_merge_content(P,output_file,file_list):
    refgenome = get_ref(P)
    v_cmd = ""
    for gvcf_file in file_list:
        v_cmd = v_cmd + f" -V {gvcf_file}"

    nMem=P.get_nMem()
    if nMem:
        mem_cmd=f"--java-options -Xmx{nMem}"
    else:
        mem_cmd=""
    cmd = [
        f"gatk {mem_cmd} CombineGVCFs   \
            -R {refgenome} \
            {v_cmd}\
            -O {output_file}"
    ]
    print("Merge gvcf files start")
    if run_cmd(cmd, f"Merge gvcf files"):
        if P.get_nosave():
            if os.path.exists(input_file):
                os.remove(input_file)


# def do_gvcf_merge(P, gvcf_file_path):
#     vcf_output_path = get_wgs_output_path(P, "vcf")
#     file_list = os.listdir(gvcf_file_path)
#     refgenome = get_ref(P)
#     gvcf_list = []
#     for file_ in file_list:
#         if file_.endswith("_g.vcf"):
#             gvcf_list.append(file_.replace("_g.vcf", ""))
#     v_cmd = ""
#     for gvcf_file in gvcf_list:
#         input_file = os.path.join(gvcf_file_path, gvcf_file)+"_g.vcf"
#         v_cmd = v_cmd + f"-V {input_file} "

#     output_file = os.path.join(vcf_output_path, "merged_g.vcf")
#     nMem=P.get_nMem()
#     if nMem:
#         mem_cmd=f"--java-options -Xmx{nMem}"
#     else:
#         mem_cmd=""
#     cmd = [
#         f"gatk {mem_cmd} CombineGVCFs   \
#             -R {refgenome} \
#             {v_cmd}\
#             -O {output_file}"
#     ]
#     print("do_gvcf_merge start")
#     if run_cmd(cmd, f"do_gvcf_merge"):
#         if P.get_nosave():
#             os.remove(input_file)

# def do_split_by_chr(P):
#     vcf_output_path = get_wgs_output_path(P, "vcf")
#     input_file = os.path.join(vcf_output_path, "merged_g.vcf")

#     vcf_file = pd.read_csv(input_file, sep='\t', comment='#', header=None)
#     # # 将第一列进行因子化处理
#     # vcf_file[0] = pd.factorize(vcf_file[0])[0] + 1
#     # group by chromosome number
#     grouped = vcf_file.groupby(0)

#     with open(input_file) as f:
#         comments = [line.strip() for line in f if line.startswith('#')]


#     # 按染色体号存放
#     chr_file_list=[]
#     for name, group in grouped:
#         chr_file_name=os.path.join(vcf_output_path,f"chr_{name}_g.vcf")
#         chr_file_list.append(chr_file_name)
#         with open(chr_file_name, 'w') as f:
#             f.write('\n'.join(comments) + '\n')
#             group.to_csv(f, sep='\t', index=False, header=None)

#     do_joint_genotype(P,chr_file_list)


def do_joint_genotype(P,chr_file_list):
    count = 0
    nThrds = P.get_nThrds()
    thread_list = []
    done_chr_file_list=[]
    for chr_file in chr_file_list:    
        output_file = chr_file.replace("_g.vcf",".vcf")
        done_chr_file_list.append(output_file)
        t = threading.Thread(target=joint_content,args=(P, chr_file,output_file))
        thread_list.append(t)
        t.start()
        if count < nThrds:
            count += 1
        if count == nThrds:
            for t in thread_list:
                t.join()
            count = 0
            thread_list = []
    for t in thread_list:
        t.join()
    return done_chr_file_list

    # zip_vcf(output_file)
    # do_index_vcf(output_file)


def joint_content(P,chr_file,output_file):
    refgenome = get_ref(P)
    nMem=P.get_nMem()
    if nMem:
        mem_cmd=f"--java-options -Xmx{nMem}"
    else:
        mem_cmd=""
    cmd = [
        f"gatk {mem_cmd} GenotypeGVCFs \
            -R {refgenome} \
            -V {chr_file} \
            -O {output_file}"
    ]
    print("Joint genotype start")
    run_cmd(cmd, f"Joint genotype")





def do_merge_chr_vcf(P, done_chr_file_list):
    vcf_output_path = get_wgs_output_path(P, "vcf")
    merge_cmd=""
    for done_chr_file in done_chr_file_list:
        merge_cmd+=f"-I {done_chr_file} "

    output_name = os.path.join(vcf_output_path, "genotype.vcf")
    cmd = [
        f'gatk MergeVcfs {merge_cmd} -O {output_name}',
    ]
    print("Merge VCF files for all chromosomes start")
    run_cmd(cmd, f"Merge VCF files for all chromosomes")

def zip_vcf(filename):
    cmd = [
        f"bgzip -f {filename}",
    ]
    print(f"zip {filename} start")
    run_cmd(cmd, f"zip {filename}")


def do_index_vcf(output_file):
    file_name = output_file+".gz"
    cmd = [
        f"tabix -p vcf {file_name}",
    ]
    print(f"Index {file_name} start")
    run_cmd(cmd, f"Index {file_name}")


#################### 10、  vcf质控 ####################
def do_vcf_qc(P, vcf_file):
    select_snp(P, vcf_file)
    select_indel(P, vcf_file)
    do_merge_snp_indel(P, vcf_file)


def select_snp(P, vcf_file):
    vcf_snp_file = vcf_file.replace(".vcf", "_snp.vcf")
    cmd = [
        f"gatk SelectVariants \
            -select-type SNP \
            -V {vcf_file} \
            -O {vcf_snp_file}"
    ]
    print(vcf_file, "select_snp start")
    run_cmd(cmd, f"{vcf_file} select_snp")
    do_snp_filter(P, vcf_snp_file)


def do_snp_filter(P, vcf_snp_file):
    output_name = vcf_snp_file.replace("_snp.vcf", "_snp_filter.vcf")
    snpQUAL = P.get_snpQUAL()
    snpQD = P.get_snpQD()
    snpMQ = P.get_snpMQ()
    snpFS = P.get_snpFS()
    snpSOR = P.get_snpSOR()
    snpMQRankSum = P.get_snpMQRankSum()
    snpReadPosRankSum = P.get_snpReadPosRankSum()

    cmd = [
        f'gatk VariantFiltration \
            -V {vcf_snp_file} \
            --filter-expression "QUAL < {snpQUAL} || QD < {snpQD} || MQ < {snpMQ} || FS > {snpFS} || SOR > {snpSOR} || MQRankSum < {snpMQRankSum} || ReadPosRankSum < {snpReadPosRankSum}" \
            --filter-name "Filter" \
            -O {output_name}',
    ]
    print(vcf_snp_file, "do_snp_filter start")
    run_cmd(cmd, f"{vcf_snp_file} do_snp_filter")


def select_indel(P, vcf_file):
    vcf_indel_file = vcf_file.replace(".vcf", "_indel.vcf")
    cmd = [
        f"gatk SelectVariants \
            -select-type INDEL \
            -V {vcf_file} \
            -O {vcf_indel_file}"
    ]
    print(vcf_file, "select_indel start")
    run_cmd(cmd, f"{vcf_file} select_indel")
    do_indel_filter(P, vcf_indel_file)


def do_indel_filter(P, vcf_indel_file):
    output_name = vcf_indel_file.replace("_indel.vcf", "_indel_filter.vcf")
    indelQUAL = P.get_indelQUAL()
    indelQD = P.get_indelQD()
    indelFS = P.get_indelFS()
    indelSOR = P.get_indelSOR()
    indelMQRankSum = P.get_indelMQRankSum()
    indelReadPosRankSum = P.get_indelReadPosRankSum()
    cmd = [
        f'gatk VariantFiltration \
            -V {vcf_indel_file} \
            --filter-expression "QUAL < {indelQUAL} || QD < {indelQD} || FS > {indelFS} || SOR > {indelSOR} || MQRankSum < {indelMQRankSum} || ReadPosRankSum < {indelReadPosRankSum}" \
            --filter-name "Filter" \
            -O {output_name}',
    ]
    print(vcf_indel_file, "do_indel_filter start")
    run_cmd(cmd, f"{vcf_indel_file} do_indel_filter")


def do_merge_snp_indel(P, vcf_file):
    vcf_snp_filter_file = vcf_file.replace(".vcf", "_snp_filter.vcf")
    vcf_indel_filter_file = vcf_file.replace(".vcf", "_indel_filter.vcf")

    output_name = vcf_file.replace(".vcf", "_filter.vcf")
    cmd = [
        f'gatk MergeVcfs \
            -I {vcf_snp_filter_file} \
            -I {vcf_indel_filter_file} \
            -O {output_name}',
    ]
    print(vcf_file, "do_merge_snp_indel start")
    success = run_cmd(cmd, f"{vcf_file} do_merge_snp_indel")

    if success:
        if P.get_nosave():
            remove_tmp_file(vcf_file)
            vcf_idx = vcf_file.replace(".vcf.gz", ".vcf.idx")
            vcf_tbi = vcf_file.replace(".vcf.gz", ".vcf.gz.tbi")
            os.remove(vcf_file)
            os.remove(vcf_idx)
            os.remove(vcf_tbi)


def remove_tmp_file(vcf_file):
    vcf_snp_file = vcf_file.replace(".vcf", "_snp.vcf")
    vcf_indel_file = vcf_file.replace(".vcf", "_indel.vcf")

    vcf_snp_filter_file = vcf_file.replace(".vcf", "_snp_filter.vcf")
    vcf_indel_filter_file = vcf_file.replace(".vcf", "_indel_filter.vcf")
    cmd = [
        f'rm -f {vcf_snp_file}* {vcf_indel_file}* {vcf_snp_filter_file}* {vcf_indel_filter_file}*',
    ]
    print(vcf_file, "remove_tmp_file start")
    run_cmd(cmd, f"{vcf_file} remove_tmp_file")


########################################################################################################################
########################################################################################################################
########################################################################################################################
#################### 第2部分 GWAS ####################
# #################### 2.0、  基因型填充 ####################
def vcf_add_id(P,vcf_file):
    vcf_input_file=vcf_file
    if vcf_input_file.endswith(".gz"):
        cmd=[
            f"bgzip -d {vcf_input_file}"
        ]
        vcf_input_file=vcf_file.replace(".gz","")
        run_cmd(cmd, f"Add id for mutation")

    vcf_reader = vcf.Reader(filename=vcf_input_file)
    add_id_vcf_file=vcf_input_file.replace(".vcf","_addid.vcf")
    vcf_writer = vcf.Writer(open(add_id_vcf_file, 'w'), vcf_reader)

    line=0
    for record in vcf_reader:
        record.ID = str(record.CHROM) + "_"+ str(record.POS)
        line+=1
        vcf_writer.write_record(record)
    vcf_writer.close()
    return add_id_vcf_file

def select_atcg(P,vcf_file):
    transvcf_output_path = get_gwas_output_path(P, "transvcf")
    output_name = os.path.join(transvcf_output_path, "only_atcg.vcf")
    only_atcg_path=os.path.join(transvcf_output_path, "only_atcg")
    cmd = [
        f"plink --vcf {vcf_file} --list-duplicate-vars ids-only suppress-first --allow-extra-chr --out {only_atcg_path}",
        f"plink --vcf {vcf_file} --alleleACGT --snps-only just-acgt --exclude {only_atcg_path}.dupvar --recode vcf --out {only_atcg_path} --allow-extra-chr"
    ]
    print("Keep the ATCG start")
    run_cmd(cmd, f"Keep the ATCG")
    return only_atcg_path+".vcf"

   
def impute_vcf(P, vcf_file):

    add_id_vcf_file = vcf_add_id(P,vcf_file)
    select_atcg_file=select_atcg(P,add_id_vcf_file)

    beagle_tool = f"{script_dir}/beagle/beagle.18May20.d20.jar"
    transvcf_output_path = get_gwas_output_path(P, "transvcf")
    output_name = os.path.join(transvcf_output_path, "genotype")
    nThrds = P.get_nThrds()
    nMem = P.get_nMem()
    impute_cmd=f"java"
    if nMem:
        impute_cmd+=f" -Xmx{nMem}"
    impute_cmd+=f" -jar {beagle_tool} gt={select_atcg_file}"
    if nThrds:
        impute_cmd+=f" nthreads={nThrds}"
    cmd = [
        f"{impute_cmd} out={output_name}",
    ]

    print("Genotype imputation start")
    run_cmd(cmd, f" Genotype imputation")

# #################### 2.1、  vcf转ped、map\bim、bed、fam ####################

def trans_vcf(P, vcf_file):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_name = os.path.join(gwasqc_output_path, "data")
    if vcf_file.endswith(".vcf") or vcf_file.endswith(".vcf.gz"):
        cmd = [
            f"plink --vcf {vcf_file} --make-bed --allow-extra-chr --out {output_name}"
        ]
    print(f"Convert {vcf_file} to bfile start")
    run_cmd(cmd, f"Convert {vcf_file} to bfile")
 

def merge_pheno_fam(P):
    phenotypefile = P.get_phenotypefile()
    if phenotypefile:
        gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
        fam = os.path.join(gwasqc_output_path, "data.fam")

        fam_data = []
        print("Merge phenotype into .fam file start")
        with open(fam, 'r') as f:
            for line in f:
                row = line.replace("\n", "").split(" ")
                fam_data.append(row)
        f.close()

        phenotype_data = []
        with open(phenotypefile, 'r') as f:
            for line in f:
                row = line.replace("\n", "").split(" ")
                if len(row)==1:
                    row = line.replace("\n", "").split("\t")
                phenotype_data.append(row)
        f.close()

        res = []
        for i in fam_data:
            for j in phenotype_data:
                if i[0] == j[0] and i[1] == j[1]:
                    temp = i
                    temp[5] = j[2]
                    res.append(temp)

        with open(fam, 'w')as f:
            for line in res:
                f.write(" ".join(line)+"\n")
        f.close()
        print("Merge phenotype into .fam file success")


# #################### 2.2、  gwas质控 ####################
def change_chr_to_num(P):
    pass
def do_gwasqc(P,bfile_dir):
    file_list = os.listdir(bfile_dir)
    bfile_name=""
    for file_ in file_list:
        if file_.endswith(".bed"):
            bfile_name = os.path.join(bfile_dir, file_.replace(".bed", ""))
            break
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    qc_bfile_file = os.path.join(gwasqc_output_path, "data")

    # snp缺失
    snpmiss = P.get_snpmiss()
    if snpmiss:
        gwasqc_snpmissingness(P,bfile_name)
        bfile_name=qc_bfile_file
    # ind缺失
    indmiss = P.get_indmiss()
    if indmiss:
        gwasqc_indmissingness(P,bfile_name)
        bfile_name=qc_bfile_file
    # 性别检查
    checksex = P.get_checksex()
    if checksex:
        gwasqc_check_sex(P,bfile_name)
        bfile_name=qc_bfile_file
    imputesex=P.get_imputesex()
    if imputesex:
        impute_sex(P,bfile_name)
        bfile_name=qc_bfile_file
    rmproblemsex=P.get_rmproblemsex()
    if rmproblemsex:
        remove_problem_sex(P,bfile_name)
        bfile_name=qc_bfile_file
    # maf
    maf=P.get_maf()
    if maf:
        gwasqc_maf(P,bfile_name)
        bfile_name=qc_bfile_file
    # hwe
    hwe=P.get_hwe()
    if hwe:
        gwasqc_hwe(P,bfile_name)
        bfile_name=qc_bfile_file   
    hweall=P.get_hweall()
    if hweall:
        gwasqc_hweall(P,bfile_name)
        bfile_name=qc_bfile_file
    
    indep=P.get_indep()
    indepPairwise=P.get_indepPairwise()
    indepPairphase=P.get_indepPairphase()
    if indep or indepPairwise or indepPairphase:
        gwasgc_LDprune(P,bfile_name)
        bfile_name=qc_bfile_file
    # heterozygosity
    heterozygosity=P.get_heterozygosity()
    if heterozygosity:
        gwasqc_heterozygosity(P,bfile_name)
        bfile_name=qc_bfile_file

    # gwasqc_removerelationship(P)
    remove_gwas_qc_tmp_file(P)

def draw_snp_missingness(P,bfile_name,flag=""):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_name = os.path.join(gwasqc_output_path, "data")
    ind_miss_file = output_name+".imiss"
    snp_miss_file = output_name+".lmiss"
    hh_file = output_name+".hh"
    log_file = output_name+".log"
    snpmiss = P.get_snpmiss()
    snp_miss_out_file = os.path.join(gwasqc_output_path, f"snpMissingness{flag}.png")
    cmd = [
        f"plink --bfile {bfile_name} --missing --allow-extra-chr --allow-no-sex --out {output_name}",
        f"Rscript --no-save {script_dir}/rscript/draw_miss_snp.R {snp_miss_file} {snpmiss} {snp_miss_out_file}",
    ]
    print("draw_snp_missingness start")
    if run_cmd(cmd, f"draw_snp_missingness"):
        if os.path.exists(hh_file):
            os.remove(hh_file)
        if os.path.exists(log_file):
            os.remove(log_file)
        if os.path.exists(ind_miss_file):
            os.remove(ind_miss_file)
        if os.path.exists(snp_miss_file):
            os.remove(snp_miss_file)

def draw_ind_missingness(P,bfile_name,flag=""):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_name = os.path.join(gwasqc_output_path, "data")
    ind_miss_file = output_name+".imiss"
    snp_miss_file = output_name+".lmiss"
    hh_file = output_name+".hh"
    log_file = output_name+".log"
    ind_miss_out_file = os.path.join(gwasqc_output_path, f"individualMissingness{flag}.png")
    indmiss = P.get_indmiss()
    cmd = [
        f"plink --bfile {bfile_name} --missing --allow-extra-chr --allow-no-sex --out {output_name}",
        f"Rscript --no-save {script_dir}/rscript/draw_miss_ind.R {ind_miss_file} {indmiss} {ind_miss_out_file}",
    ]
    print("draw_ind_missingness start")
    if run_cmd(cmd, f"draw_ind_missingness"):
        if os.path.exists(hh_file):
            os.remove(hh_file)
        if os.path.exists(log_file):
            os.remove(log_file)
        if os.path.exists(ind_miss_file):
            os.remove(ind_miss_file)
        if os.path.exists(snp_miss_file):
            os.remove(snp_miss_file)


def gwasqc_snpmissingness(P,bfile_name):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_name = os.path.join(gwasqc_output_path, "data")
    snpmiss = P.get_snpmiss()
    cmd = [
        f"plink --bfile {bfile_name} --geno {snpmiss} --make-bed --allow-extra-chr --allow-no-sex --out {output_name}",
    ]
    draw_snp_missingness(P,bfile_name,f"_before_{snpmiss}")
    print("gwasqc_snpmissingness start")
    if run_cmd(cmd, f"gwasqc_snpmissingness"):
        draw_snp_missingness(P,output_name,f"_after_{snpmiss}")

def gwasqc_indmissingness(P,bfile_name):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_name = os.path.join(gwasqc_output_path, "data")
    indmiss = P.get_indmiss()
    cmd = [
        f"plink --bfile {bfile_name} --mind {indmiss} --make-bed --allow-extra-chr --allow-no-sex --out {output_name}",
    ]
    draw_ind_missingness(P,bfile_name,f"_before_{indmiss}")
    print("gwasqc_indmissingness start")
    if run_cmd(cmd, f"gwasqc_indmissingness"):
        draw_ind_missingness(P,output_name,f"_after_{indmiss}")


# bed文件中要有sex那列或者fam文件中要有
def gwasqc_check_sex(P,bfile_name,flag="",remove=True):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    sexcheck_file = output_file+".sexcheck"
    gender_check_file = os.path.join(gwasqc_output_path, f"gender_check{flag}.png")

    cmd = [
        f"plink --bfile {bfile_name} --check-sex --out {output_file}",
        f"Rscript --no-save {script_dir}/rscript/draw_gender_check.R {sexcheck_file} {gender_check_file}",
    ]

    print("gwasqc_check_sex start")
    if run_cmd(cmd, f"gwasqc_check_sex") and remove:
        if os.path.exists(sexcheck_file):
            os.remove(sexcheck_file)

def impute_sex(P,bfile_name):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")

    cmd = [
        f'plink --bfile {bfile_name} --impute-sex --make-bed --out {output_file}',
    ]

    gwasqc_check_sex(P,bfile_name,flag="_before_impute")
    print("impute_sex start")
    if run_cmd(cmd, f"impute_sex"):
        gwasqc_check_sex(P,bfile_name,flag="_after_impute")

def remove_problem_sex(P,bfile_name):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    problem_sex = os.path.join(gwasqc_output_path, "problem_sex.txt")
    sexcheck_file = output_file+".sexcheck"
    print_cmd="'{print$1,$2}'"
    cmd = [
        f"grep 'PROBLEM'{sexcheck_file} | awk {print_cmd}> {problem_sex}",
        f"plink --bfile {bfile_name} --remove {{problem_sex}} --make-bed --out {output_file}"
    ]

    gwasqc_check_sex(P,bfile_name,flag="_before_remove",remove=False)
    print("remove_problem_sex start")
    if run_cmd(cmd, f"remove_problem_sex"):
        gwasqc_check_sex(P,bfile_name,flag="after_remove")

# 过滤稀有突变
def gwasqc_maf(P,bfile_name):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    maf = P.get_maf()

    cmd = [
        f"plink --bfile {bfile_name} --maf {maf} --make-bed --allow-extra-chr --allow-no-sex --out {output_file}",
    ]
    draw_maf(P,bfile_name,flag="_before_filter")
    print("gwasqc_maf start")
    if run_cmd(cmd, f"gwasqc_maf"):
        draw_maf(P,output_file,flag="_after_filter")

def draw_maf(P,bfile_name,flag=""):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    maf_out = os.path.join(gwasqc_output_path, "maf_check")
    maf = P.get_maf()
    maf_frq_file = maf_out+".frq"
    maf_hh_file = maf_out+".hh"
    maf_log_file = maf_out+".log"
    maf_nosex_file = maf_out+".nosex"
    maf_distribution_file = os.path.join(gwasqc_output_path, f"maf_distribution{flag}.png")

    cmd = [
        f"plink --bfile {bfile_name} --freq --allow-extra-chr --allow-no-sex --out {maf_out}",
        f"Rscript --no-save {script_dir}/rscript/draw_maf_check.R {maf_frq_file} {maf} {maf_distribution_file} ",
    ]
    print("gwasqc_maf start")
    if run_cmd(cmd, f"gwasqc_maf"):
        if os.path.exists(maf_hh_file):
            os.remove(maf_hh_file)
        if os.path.exists(maf_frq_file):
            os.remove(maf_frq_file)
        if os.path.exists(maf_log_file):
            os.remove(maf_log_file)
        if os.path.exists(maf_nosex_file):
            os.remove(maf_nosex_file)

def gwasqc_hwe(P,bfile_name):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    hwe = P.get_hwe()

    cmd = [
        f"plink --bfile {bfile_name} --hwe {hwe} --make-bed --allow-extra-chr --allow-no-sex --out {output_file}",
    ]
    draw_hwe(P,bfile_name,flag="_before_filter")
    print("gwasqc_hwe start")
    if run_cmd(cmd, f"gwasqc_hwe"):
        draw_hwe(P,bfile_name,flag="_after_filter")

def gwasqc_hweall(P,bfile_name):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    hweall = P.get_hweall()

    cmd = [
        f"plink --bfile {bfile_name} --hwe {hweall} --hwe-all --make-bed --allow-extra-chr --allow-no-sex --out {output_file}",
    ]
    draw_hwe(P,bfile_name,flag="_before_filter")
    print("gwasqc_hwe start")
    if run_cmd(cmd, f"gwasqc_hwe"):
        draw_hwe(P,bfile_name,flag="_after_filter")

def draw_hwe(P,bfile_name,flag=""):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    hwe = P.get_hwe()
    hwe_file = bfile_name+".hwe"
    hwe_hh_file = bfile_name+".hh"
    hwe_log_file = bfile_name+".log"
    hwe_nosex_file = bfile_name+".nosex"
    hwe_out_file = os.path.join(gwasqc_output_path, f"hwe{flag}.png")


    cmd = [
        f"plink --bfile {bfile_name} --hardy --allow-extra-chr --allow-no-sex --out {output_file}",
        f"Rscript --no-save {script_dir}/rscript/draw_hwe.R {hwe_file} {hwe} {hwe_out_file}",
    ]

    print("draw_hwe start")
    if run_cmd(cmd, f"draw_hwe"):
        if os.path.exists(hwe_file):
            os.remove(hwe_file)
        if os.path.exists(hwe_hh_file):
            os.remove(hwe_hh_file)
        if os.path.exists(hwe_log_file):
            os.remove(hwe_log_file)
        if os.path.exists(hwe_nosex_file):
            os.remove(hwe_nosex_file)


def gwasgc_LDprune(P,bfile_name):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    indepSNP = os.path.join(gwasqc_output_path, "data.prune.in")

    indep_cmd=""
    if P.get_indep():
        indep=' '.join(P.get_indep())
        indep_cmd=f"--indep {indep}"
    if P.get_indepPairwise():
        indepPairwise=' '.join(P.get_indepPairwise())
        indep_cmd=f"--indep-pairwise {indepPairwise}"
    if P.get_indepPairphase():
        indepPairphase=' '.join(P.get_indepPairphase())
        indep_cmd=f"--indep-pairphase {indepPairphase}"
    temp_log_file = output_file+".log"
    temp_out_file = output_file+".prune.out"
    temp_hh_file = output_file+".hh"
    temp_nosex_file = output_file+".nosex"

    cmd = [
        f"plink --bfile {bfile_name} {indep_cmd} --allow-extra-chr --allow-no-sex --out {output_file}",
        f"plink --bfile {output_file} --make-bed --extract {indepSNP} --allow-extra-chr --allow-no-sex --out {output_file}",
    ]
    if run_cmd(cmd, f"gwasgc_LDprune"):
        # if os.path.exists(indepSNP):
        #     os.remove(indepSNP)
        if os.path.exists(temp_log_file):
            os.remove(temp_log_file)
        if os.path.exists(temp_out_file):
            os.remove(temp_out_file)
        if os.path.exists(temp_hh_file):
            os.remove(temp_hh_file)
        if os.path.exists(temp_nosex_file):
            os.remove(temp_nosex_file)



def gwasqc_heterozygosity(P,bfile_name):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    het_file = output_file+".het"
    all_out_het = os.path.join(gwasqc_output_path, "all_out_het.txt")
    heterozygosity = P.get_heterozygosity()

    out_het_snp = os.path.join(gwasqc_output_path, "out_het_snp.txt")

    awk_cmd = "{print $1,$2}"

    cmd = [
        f"Rscript --no-save {script_dir}/rscript/heterozygosity_outliers_list.R {het_file} {all_out_het} {heterozygosity}",
        f"sed 's/\"//g' {all_out_het} | awk '{awk_cmd}' > {out_het_snp}",
        f"plink --bfile {bfile_name} --remove {out_het_snp} --make-bed --allow-extra-chr --allow-no-sex --out {output_file}",
    ]

    draw_heterozygosity(P,bfile_name,flag="_before_filter",remove=False)
    print("gwasqc_heterozygosity start")
    if run_cmd(cmd, f"gwasqc_heterozygosity"):
        draw_heterozygosity(P,bfile_name,flag="_after_filter")
        if os.path.exists(out_het_snp):
            os.remove(out_het_snp)
        if os.path.exists(het_file):
            os.remove(het_file)






def draw_heterozygosity(P,bfile_name,flag="",remove=True):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    output_file = os.path.join(gwasqc_output_path, "data")
    heterozygosity = P.get_heterozygosity()
    het_file = output_file+".het"
    temp_log_file = output_file+".log"
    temp_nosex_file = output_file+".nosex"


    het_out_file = os.path.join(gwasqc_output_path, f"heterozygosity{flag}.png")
 

    cmd = [
        f"plink --bfile {output_file} --het --allow-extra-chr --allow-no-sex --out {output_file}",
        f"Rscript --no-save {script_dir}/rscript/draw_heterozygosity.R {het_file} {heterozygosity} {het_out_file}",
    ]

    print("draw_heterozygosity start")
    if run_cmd(cmd, f"draw_heterozygosity") and remove:
        if os.path.exists(temp_log_file):
            os.remove(temp_log_file)
        if os.path.exists(temp_nosex_file):
            os.remove(temp_nosex_file)


# def gwasqc_removerelationship(P,bfile_name):
#     if P.get_removerelationship():
#         gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
#         input_file = os.path.join(gwasqc_output_path, "data")
#         pihat = P.get_pihat()
#         pihat_out = os.path.join(
#             gwasqc_output_path, "pihat_min_"+pihat.replace(".", "_"))
#         awk_cmd = "{if ($8 > 0.9) print $0}"
#         pihat_genome = pihat_out+".genome"
#         zoom_pihat = os.path.join(gwasqc_output_path, "zoom_pihat.genome")

#         relatedness = os.path.join(gwasqc_output_path, "relatedness.png")
#         zoom_relatedness = os.path.join(
#             gwasqc_output_path, "zoom_relatedness.png")
#         hist_relatedness = os.path.join(
#             gwasqc_output_path, "hist_relatedness.png")

#         cmd = [
#             f"plink --bfile {input_file} --genome --min {pihat} --allow-extra-chr --allow-no-sex --out {pihat_out}",
#             f"awk '{awk_cmd}' {pihat_genome} > {zoom_pihat}",
#             f"Rscript --no-save gwaswa/rscript/Relatedness.R {relatedness} {pihat_genome} {zoom_relatedness} {zoom_pihat} {hist_relatedness} {pihat_genome}",
#             f"plink --bfile {input_file} --filter-founders --make-bed --allow-extra-chr --allow-no-sex --out {input_file}"
#         ]

#         print("gwasqc_removerelationship start")
#         run_cmd(cmd, f"gwasqc_removerelationship")


def remove_gwas_qc_tmp_file(P):
    gwasqc_output_path = get_gwas_output_path(P, "gwasqc")
    file_list = os.listdir(gwasqc_output_path)
    for file_ in file_list:
        if file_.endswith("~"):
            os.remove(os.path.join(gwasqc_output_path, file_))


# #################### 2.3、  pca分析 ####################
def do_population(P,bfile_dir):
    file_list = os.listdir(bfile_dir)
    for file_ in file_list:
        if file_.endswith(".bed"):
            bfile_name = os.path.join(bfile_dir, file_.replace(".bed", ""))
            break
    do_pca(P,bfile_dir)    
    do_admixture(P,bfile_dir)
    get_bestk_clean_cve(P)
    # set_admixture_info(P) 不需要了
    draw_cv_error(P)
    draw_admixture(P,bfile_name)
    draw_pca(P,bfile_name)

def do_pca(P, bfile_dir):
    file_list = os.listdir(bfile_dir)
    pcanum = P.get_pcanum()
    for file_ in file_list:
        if file_.endswith(".bed"):
            input_file = os.path.join(bfile_dir, file_.replace(".bed", ""))
            break

    pca_output_path = get_gwas_output_path(P, "pca")
    output_name = os.path.join(pca_output_path, "pca")

    eigenvec = output_name+".eigenvec"
    eigenvec_filename = os.path.join(pca_output_path, "eigenvec.xls")
    eigenval = output_name+".eigenval"
    eigenval_filename = os.path.join(pca_output_path, "eigenval.xls")

    cmd = [
        f"plink --bfile {input_file} --pca {pcanum} --allow-extra-chr --allow-no-sex --out {output_name}",
        f"Rscript {script_dir}/rscript/r_pca.R {eigenvec} {eigenvec_filename} {eigenval} {eigenval_filename}"
    ]
    print("do_pca start")
    run_cmd(cmd, f"do_pca")


def do_admixture(P, bfile_dir):
    file_list = os.listdir(bfile_dir)
    for file_ in file_list:
        if file_.endswith(".bed"):
            input_file = os.path.join(bfile_dir, file_)
            break
    groupnum = P.get_groupnum()
    nThrds = P.get_nThrds()
    pca_output_path = os.path.abspath(get_gwas_output_path(P, "pca"))

    input_file = os.path.abspath(input_file)
    if groupnum:
        out_path = os.path.join(pca_output_path, f"log{groupnum}.out")
        admixture_cmd=f"admixture --cv {input_file} {groupnum}"
        if nThrds:
            admixture_cmd+=f" -j{nThrds}"
        cmd = [
            f"cd {pca_output_path}",
            f"{admixture_cmd} | tee {out_path}",
        ]
        print("do_admixture start")
        run_cmd(cmd, f"do_admixture")
    else:
        for i in range(2, 20):
            out_path = os.path.join(pca_output_path, f"log{i}.out")
            admixture_cmd=f"admixture --cv {input_file} {i}"
            if nThrds:
                admixture_cmd+=f" -j{nThrds}"
            cmd = [
                f"cd {pca_output_path}",
                f"{admixture_cmd} | tee {out_path}",
            ]
            print(f"do_admixture_groupnum{i} start")
            run_cmd(cmd, f"do_admixture_groupnum{i}")


# 获取最佳k值,和clean_cve文件
def get_bestk_clean_cve(P):
    pca_output_path = get_gwas_output_path(P, "pca")
    cross_validation_error_file = os.path.join(
        pca_output_path, "cross_validation_error.txt")
    groupnum = P.get_groupnum()
    if not groupnum:
        out_path = os.path.join(pca_output_path, "log*.out")
    else:
        out_path = os.path.join(pca_output_path, f"log{groupnum}.out")
    cmd = [
        f"grep -h CV {out_path} |sort -nk4  -t ' ' > {cross_validation_error_file}"
    ]
    print("grep_CV start")
    run_cmd(cmd, f"grep_CV")

    # 读cverror，提取 k 和 value
    cvef = open(cross_validation_error_file, 'r')
    cvef_lines = cvef.readlines()
    k_v_list = ["K\tError\n"]
    bestk = 0
    bestv = 9999
    for line in cvef_lines:
        k, v = line.split("K=")[1].split("): ")
        k_v_list.append("\t".join([k, v]))
        if float(v) < bestv:
            bestv = float(v)
            bestk = k
    cvef.close()
    # 将 k、v写到clean文件中
    clean_cve_file = os.path.join(pca_output_path, f"cve_clean_{bestk}.txt")
    cvef_clean = open(clean_cve_file, "w")
    for line in k_v_list:
        cvef_clean.write(line+'\n')
    cvef_clean.close()
    print("cvef_clean success")


# # 设置每个样本的分组、颜色、形状，用来画pca的图
    # def set_admixture_info(P):
    #     pca_output_path = get_gwas_output_path(P, "pca")
    #     file_list = os.listdir(pca_output_path)
    #     for file_ in file_list:
    #         if file_.startswith("cve_clean_"):
    #             clean_cve_file = os.path.join(pca_output_path, file_)
    #             break
    #     bestk = re.search("\d+$", clean_cve_file.replace(".txt", "")).group()
    #     for file_ in file_list:
    #         if file_.endswith(f"{bestk}.Q"):
    #             q_file = os.path.join(pca_output_path, file_)
    #             break
    #     out_file = os.path.join(pca_output_path, "best_adinfo.csv")
    #     cmd = [
    #         f"Rscript gwaswa/rscript/r_adinfo.R {q_file} {out_file}",
    #     ]
    #     print("set_admixture_info start")
    #     run_cmd(cmd, f"set_admixture_info")


# 画出各个k的err值 折线图
def draw_cv_error(P):
    pca_output_path = get_gwas_output_path(P, "pca")
    file_list = os.listdir(pca_output_path)
    # 编写命令行
    for file_ in file_list:
        if file_.startswith("cve_clean_"):
            input_file = os.path.join(pca_output_path, file_)
            break
    cve_image = os.path.join(pca_output_path, "cve.png")
    cmd = [
        f"Rscript {script_dir}/rscript/draw_cv_error.R {input_file} {cve_image}"
    ]
    print("draw_cv_error start")
    run_cmd(cmd, f"draw_cv_error")

# 画出群体结构 柱状图
def draw_admixture(P,bfile_name):
    pca_output_path = get_gwas_output_path(P, "pca")
    file_list = os.listdir(pca_output_path)
    fam1=os.path.join(pca_output_path,"fam1.txt")
    
    for file_ in file_list:
        if file_.endswith(".Q"):
            q_file = os.path.join(pca_output_path, file_)
            k = re.search("\d+$", q_file.replace(".Q", "")).group()
            out_file = q_file.replace(f".{k}.Q", f"_{k}_admixture.png")
            cmd = [
                f"cut -f 1 -d ' ' {bfile_name}.fam > {fam1}",
                f"Rscript {script_dir}/rscript/draw_admixture.R {fam1} {q_file} {k} {out_file} "
            ]
            print("draw_admixture start")
            run_cmd(cmd, f"draw_admixture")


def draw_pca(P,bfile_name):
    pca_output_path = get_gwas_output_path(P, "pca")
    eigenvec_xls = os.path.join(pca_output_path, "eigenvec.xls")
    eigenval_xls = os.path.join(pca_output_path, "eigenval.xls")
    pca_image = os.path.join(pca_output_path, "pca.png")

    file_list = os.listdir(pca_output_path)
    for file_ in file_list:
        if file_.startswith("cve_clean_"):
            clean_cve_file = os.path.join(pca_output_path, file_)
            break
    bestk = re.search("\d+$", clean_cve_file.replace(".txt", "")).group()

    best_Q_file=os.path.join(pca_output_path,os.path.basename(bfile_name))+f".{bestk}.Q"


    # 用R绘图
    cmd = [
        f"Rscript {script_dir}/rscript/draw_pca.R {eigenvec_xls} {eigenval_xls} {best_Q_file} {bestk} {pca_image}"
    ]
    print("draw_pca start")
    run_cmd(cmd, f"draw_pca")


##################################################    kinship ########


# 计算亲缘关系
def do_bed_to_grm(P, clean_file_path):
    file_list = os.listdir(clean_file_path)
    for file_ in file_list:
        if file_.endswith(".bed"):
            input_file = os.path.join(
                clean_file_path, file_.replace(".bed", ""))
            break

    kinship_output_path = get_gwas_output_path(P, "kinship")
    output_file = os.path.join(kinship_output_path, "kinship")
    cmd = [
        f"gcta64 --make-grm-gz --bfile {input_file} --out {output_file} --autosome"
    ]
    print("do_bed_to_grm start")
    run_cmd(cmd, f"do_bed_to_grm")

# 整理kinship


def arrange_kinship(P):
    kinship_output_path = get_gwas_output_path(P, "kinship")

    grm_gz = os.path.join(kinship_output_path, "kinship.grm.gz")
    grm_id = os.path.join(kinship_output_path, "kinship.grm.id")

    out_file = os.path.join(kinship_output_path, "kinship.txt")
    cmd = [
        f"Rscript {script_dir}/rscript/r_kinship.R {grm_gz} {grm_id} {out_file}"
    ]
    print("arrange_kinship start")
    run_cmd(cmd, f"arrange_kinship")


# 画出kinship热图
def draw_kinship(P):
    kinship_output_path = get_gwas_output_path(P, "kinship")

    kinship_file = os.path.join(kinship_output_path, "kinship.txt")
    out_image = os.path.join(kinship_output_path, "kinship.png")
    cmd = [
        f"Rscript {script_dir}/rscript/draw_kinship.R {kinship_file} {out_image}",
    ]
    print("draw_kinship start")
    run_cmd(cmd, f"draw_kinship")


# #################### 2.5、  关联分析 ####################
def do_associate_lm(P, clean_file_path):
    file_list = os.listdir(clean_file_path)
    for file_ in file_list:
        if file_.endswith(".bed"):
            input_file = os.path.join(
                clean_file_path, file_.replace(".bed", ""))
            break
    association_output_path = get_gwas_output_path(
        P, os.path.join("association", "lm"))

    assoc_txt = os.path.join(association_output_path, "result.assoc.txt")
    man_image = os.path.join(association_output_path, "man.png")
    qq_image = os.path.join(association_output_path, "qq.png")
    cmd = [
        f"gemma -bfile {input_file} -lm 1 -outdir {association_output_path}",
        f"Rscript {script_dir}/rscript/draw_man.R {assoc_txt} {man_image} {qq_image}",
    ]
    print("do_associate_lm start")
    run_cmd(cmd, f"do_associate_lm")


def do_associate_lmm(P, clean_file_path):
    file_list = os.listdir(clean_file_path)
    for file_ in file_list:
        if file_.endswith(".bed"):
            input_file = os.path.join(
                clean_file_path, file_.replace(".bed", ""))
            break    
    association_output_path = get_gwas_output_path(P, os.path.join("association", "lmm"))
    assoc_txt = os.path.join(association_output_path, "result.assoc.txt")
    man_image = os.path.join(association_output_path, "man.png")
    qq_image = os.path.join(association_output_path, "qq.png")

    pcafile = P.get_pcafile()
    kinshipfile = get_kinship_file(P)
    lmm_cmd=f"gemma -bfile {input_file}"
    if pcafile:
        merge_pca_fam_to_covar(P, pcafile)
        covar = os.path.join(association_output_path, "covar.txt")
        lmm_cmd+=f" -c {covar}"
    if kinshipfile:
        lmm_cmd+=f" -k {kinshipfile}"
    cmd = [
        f"gemma -bfile {input_file} -gk 1 -outdir {association_output_path}",
        f"{lmm_cmd} -lmm 1 -outdir {association_output_path}",
        f"Rscript {script_dir}/rscript/draw_man.R {assoc_txt} {man_image} {qq_image}",
    ]
    print("do_associate_lmm start")
    run_cmd(cmd, f"do_associate_lmm")



# def get_pca_file(P):
#     pcafile = P.get_pcafile()
#     if pcafile:
#         pca_file = pcafile
#     else:
#         pca_output_path = get_gwas_output_path(P, "pca")
#         pca_file = os.path.join(pca_output_path, "pca.eigenvec")

#     if os.path.exists(pca_file):
#         return pca_file
#     return None

def get_kinship_file(P):
    kinshipfile = P.get_kinshipfile()
    if kinshipfile:
        kinship_file = kinshipfile
    else:
        kinship_output_path = get_gwas_output_path(P, os.path.join("association", "lmm"))
        kinship_file = os.path.join(kinship_output_path, "result.cXX.txt")
    return kinship_file

def merge_pca_fam_to_covar(P, pca_file):

    lmm_output_path = get_gwas_output_path(P, os.path.join("association", "lmm"))
    covar = os.path.join(lmm_output_path, "covar.txt")

    qc_output_path = get_gwas_output_path(P, "gwasqc")
    fam = os.path.join(qc_output_path, "data.fam")

    new = {}
    # with open(fam, 'r') as f:
    #     for line in f:
    #         tmp = ["1"]
    #         row = line.replace("\n", "").split(" ")
    #         tmp.append(row[4])
    #         new[row[1]] = tmp
    # f.close()

    with open(pca_file, 'r') as f:
        for line in f:
            row = line.replace("\n", "").split(" ")
            # new[row[1]].extend(row[2:])
            new[row[1]] = row[2:]
    f.close()

    # print(new.values())

    with open(covar, 'w')as f:
        for line in new.values():
            f.write(" ".join(line)+"\n")
    f.close()
    print("merge_pca_fam_to_covar success")



def do_selectsnp(P):
    assocfile = P.get_assocfile()
    pvaluelimit = P.get_pvaluelimit()
    snp_output_path = get_gwas_output_path(P,"selectsnp")
    
    remarkable_snp = os.path.join(snp_output_path, f"remarkable_snp_{pvaluelimit}.txt")
    variant_file = os.path.join(snp_output_path, f"remarkable_variant_{pvaluelimit}.vcf")
    
    cmd = [
        f"Rscript {script_dir}/rscript/select_snp.R {assocfile} {remarkable_snp} {pvaluelimit}",
        f"Rscript {script_dir}/rscript/make_variant.R {remarkable_snp} {variant_file}",
    ]
    print("do_selectsnp start")
    run_cmd(cmd, f"do_selectsnp")


########################################################################################################################
########################################################################################################################
########################################################################################################################
#################### 第3部分 预测 ####################
# #################### 3.1、  enformer ####################
def check_hg38(fasta_file):
    if not os.path.exists(fasta_file):
        cmd = [
            f"wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > {fasta_file}",
        ]
        print("download hg38 start")
        run_cmd(cmd, f"download hg38")

def enformer_predict(P, file_path):
    fasta_file = f"{script_dir}/enformer/data/hg38.fa"
    check_hg38(fasta_file)
    model_path = f"{script_dir}/enformer/data/enformer_1"
    transform_path = f'{script_dir}/enformer/data/trained_model.pkl'
    targets_txt = f'{script_dir}/enformer/data/targets_human.txt'
    enformer_output_path = get_assessment_output_path(P, "enformer")
    model = Enformer(model_path)
    fasta_extractor = FastaStringExtractor(fasta_file)

    variant_list = []
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                temp=line.replace("\n", "").split(" ")
                if len(temp)>1:
                    variant_list.append(temp)
                else:
                    # print(temp)
                    variant_list.append(line.replace("\n", "").split("\t"))
    f.close()

    snp_tracks={}
    for chrm, ps, rs, ref, alt in variant_list:
        if chrm.startswith("ch"):
            chrom = chrm
        else:
            chrom = "chr"+chrm

        variant = kipoiseq.Variant(chrom, int(ps), ref, alt, id=rs)
        # variant_path = os.path.join(enformer_output_path, f"{chrom}_{ps}")
        # os.makedirs(variant_path, exist_ok=True)
        print("assessing", variant)
        interval = kipoiseq.Interval(variant.chrom, variant.start, variant.start).resize(SEQUENCE_LENGTH)
        seq_extractor = kipoiseq.extractors.VariantSeqExtractor(reference_sequence=fasta_extractor)
        center = interval.center() - interval.start

        reference = seq_extractor.extract(interval, [], anchor=center)
        alternate = seq_extractor.extract(interval, [variant], anchor=center)

        reference_prediction = model.predict_on_batch(one_hot_encode(reference)[np.newaxis])['human'][0]
        alternate_prediction = model.predict_on_batch(one_hot_encode(alternate)[np.newaxis])['human'][0]

        variant_track = np.zeros_like(reference_prediction[:, 0], dtype=bool)
        variant_track[variant_track.shape[0] // 2] = True


        target = pd.read_table(targets_txt)
        for i in range(len(target)):
            tracks = [
                reference_prediction[:, i].tolist(),
                alternate_prediction[:, i].tolist(),
                (alternate_prediction[:, i]-reference_prediction[:, i]).tolist()
            ]
            snp_dir=os.path.join(enformer_output_path,f"{variant.chrom}_{variant.pos}")
            if not os.path.exists(snp_dir):
                os.makedirs(snp_dir)
            with open(os.path.join(snp_dir,f'{i}.js'), 'w') as f:
                f.write('var tracks=')
                json.dump(tracks, f)
            f.close()
        
    print("assessment summary done")


    # out_csv_file=os.path.join(enformer_output_path, "allFeaturesZscore.csv")
    # variant_score(model_path, transform_path,file_path,out_csv_file,targets_txt)

# def variant_score(model_path, transform_path,file_path,out_csv_file,targets_txt):
#     SEQUENCE_LENGTH = 393216
#     enformer_score_variants_all = EnformerScoreVariantsNormalized(model_path, transform_path)
#     it = variant_centered_sequences(file_path, sequence_length=SEQUENCE_LENGTH,
#                                     gzipped=False, chr_prefix='chr')
#     example_list = []
#     df_targets = pd.read_csv(targets_txt, sep='\t')
#     for i, example in enumerate(it):
#         if i >= 20:
#             break
#         variant_scores = enformer_score_variants_all.predict_on_batch(
#             {k: v[tf.newaxis] for k, v in example['inputs'].items()})[0]
#         variant_scores = {f'{i}_{name[:20]}': score for i, (name, score) in enumerate(
#             zip(df_targets.description, variant_scores))}
#         example_list.append({**example['metadata'],
#                             **variant_scores})
#         if i % 2 == 0:
#             print(f'Done {i}')
#     df = pd.DataFrame(example_list)
#     df.to_csv(out_csv_file)
#     print("Done")


def predict(model, fasta_extractor, variant,targets_txt, variant_path):
    # 以变异位点为中心，前后取SEQUENCE_LENGTH长度的序列
    interval = kipoiseq.Interval(variant.chrom, variant.start, variant.start).resize(SEQUENCE_LENGTH)
    seq_extractor = kipoiseq.extractors.VariantSeqExtractor(reference_sequence=fasta_extractor)
    center = interval.center() - interval.start

    reference = seq_extractor.extract(interval, [], anchor=center)
    alternate = seq_extractor.extract(interval, [variant], anchor=center)

    reference_prediction = model.predict_on_batch(one_hot_encode(reference)[np.newaxis])['human'][0]
    alternate_prediction = model.predict_on_batch(one_hot_encode(alternate)[np.newaxis])['human'][0]

    variant_track = np.zeros_like(reference_prediction[:, 0], dtype=bool)
    variant_track[variant_track.shape[0] // 2] = True


    target = pd.read_table(targets_txt)
    for i in range(len(target)):
        tracks = [
            reference_prediction[:, i].tolist(),
            alternate_prediction[:, i].tolist(),
            (alternate_prediction[:, i]-reference_prediction[:, i]).tolist()
        ]
        snp_dir=f"{script_dir}/enformer/{variant.chrom}_{variant.pos}"
        if not os.path.exists(snp_dir):
            os.makedirs(snp_dir)
        with open(f'{script_dir}/enformer/{variant.chrom}_{variant.pos}/{i}.js', 'w') as f:
            f.write('var tracks=')
            json.dump(tracks, f)
        


# def plot_tracks(tracks, interval, out_file, height=1,):
#     fig, axes = plt.subplots(len(tracks), 1, figsize=(20, height * len(tracks)), sharex=True)
#     for ax, (title, y) in zip(axes, tracks.items()):
#         # 将预测值放在区间内
#         ax.fill_between(np.linspace(interval.start, interval.end, num=len(y)), y)
#         # 字典的key
#         ax.set_title(title)
#         sns.despine(top=True, right=True, bottom=True)
#     ax.set_xlabel(str(interval))
#     plt.tight_layout()
#     plt.savefig(out_file)
#     plt.close()

# #################### 3.2、  vep ####################


def vep_predict(P, file_path):
    vep_output_path = get_assessment_output_path(P, "vep")
    output_file = os.path.join(vep_output_path, os.path.basename(
        file_path).replace(".vcf", "_vep.out.vcf"))

    nThrds = P.get_nThrds()
    vep_cmd = f"vep -i {file_path} -o {output_file} "
    if nThrds:
        vep_cmd += f" --fork {nThrds}"
    species=P.get_species()
    assembly=P.get_assembly()
    if species=="homo_sapiens":
        if assembly:
            cmd = [
                f"{vep_cmd} --force_overwrite --symbol --sift b --polyphen b  --database  --species {species} --assembly {assembly}",
            ]
        else:
            cmd = [
                f"{vep_cmd} --force_overwrite --symbol --sift b --polyphen b  --database  --species {species}",
            ]
    else:
        if assembly:
            cmd = [
                f"{vep_cmd} --force_overwrite --symbol --sift b --is_multispecies 1 --database  --species {species} --genomes --assembly {assembly}",
            ]
        else:
            cmd = [
                f"{vep_cmd} --force_overwrite --symbol --sift b --is_multispecies 1 --database  --species {species} --genomes ",
            ]

    print("vep_predict start")
    run_cmd(cmd, f"vep_predict")


########################################################################################################################
###############################################         VEP        #####################################################
# 这里数据做：统计
def getConsequenceDictAndsiftAndPolyPhenDict(vcffile,snpfile):
    snpList = getSNPList(snpfile)
    vepout = open(vcffile, 'r')

    consequenceDict = {
        "transcript_ablation": [],
        "splice_acceptor": [],
        "splice_donor": [],
        "stop_gained": [],
        "frameshift": [],
        "stop_lost": [],
        "start_lost": [],
        "transcript_amplification": [],
        "inframe_insertion": [],
        "inframe_deletion": [],
        "missense": [],
        "protein_altering": [],
        "splice_region": [],
        "splice_donor_5th_base": [],
        "splice_donor_region": [],
        "splice_polypyrimidine_tract": [],
        "incomplete_terminal_codon": [],
        "start_retained": [],
        "stop_retained": [],
        "synonymous": [],
        "coding_sequence": [],
        "mature_miRNA": [],
        "5_prime_UTR": [],
        "3_prime_UTR": [],
        "non_coding_transcript_exon": [],
        "intron": [],
        "NMD_transcript": [],
        "non_coding_transcript": [],
        "upstream_gene": [],
        "downstream_gene": [],
        "TFBS_ablation": [],
        "TFBS_amplification": [],
        "TF_binding_site": [],
        "regulatory_region_ablation": [],
        "regulatory_region_amplification": [],
        "feature_elongation": [],
        "regulatory_region": [],
        "feature_truncation": [],
        "intergenic": [],
    }
    siftAndPolyPhenDict = {
        "deleterious": [],
        "deleterious_low_confidence": [],
        "tolerated_low_confidence": [],
        "tolerated": [],
        "probably_damaging": [],
        "possibly_damaging": [],
        "benign": [],
        "unknown": [],
    }
    for line in vepout:
        if not line.startswith('##') and not line.startswith('#'):
            # Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
            columns = line.strip().split('\t')
            snpid = columns[0]
            pos = columns[1]
            position = "chr"+pos
            allele = columns[2]
            affected_Gene = columns[3]
            feature = columns[4]
            feature_type = columns[5]
            consequence = columns[6].replace("_variant", "")
            cDNA_position = columns[7]
            cDS_position = columns[8]
            protein_position = columns[9]
            amino_acids = columns[10]
            codons = columns[11]
            existing_variation = columns[12]
            extra = columns[13]


            impact = ""
            pattern = r'IMPACT=[^;]*'
            match = re.search(pattern, extra)
            if match:
                impact = match.group(0).replace("IMPACT=", "")

            symbol = ""
            pattern = r'SYMBOL=[^;]*'
            match = re.search(pattern, extra)
            # print(match)
            if match:
                symbol = match.group(0).replace("SYMBOL=", "")

            sift = ""
            pattern = r'SIFT=[^;]*'
            match = re.search(pattern, extra)
            if match:
                sift = match.group(0).replace("SIFT=", "").replace("("," (")

            polyphen = ""
            pattern = r'PolyPhen=[^;]*'
            match = re.search(pattern, extra)
            if match:
                polyphen = match.group(0).replace("PolyPhen=", "").replace("("," (")
            
            for i in snpList:
                if i[1]==snpid:
                    ref = i[2]
                    alt = i[3]

            title=[
                "SNP_ID",
                "Position",
                "Ref",
                "Alt",
                "Affected_Gene_ID",
                "Gene_Symbol",
                "Feature_ID",
                "Feature_Type",
                "Impact",
                "SIFT",
                "PolyPhen",
                "Consequence",
                ]
            
            obj = {
                "SNP_ID": snpid,
                "Position": position,
                "Ref":ref,
                "Alt":alt,
                "Affected_Gene_ID": affected_Gene,
                "Gene_Symbol": symbol,
                "Feature_ID": feature,
                "Feature_Type": feature_type,
                "Impact": impact,
                "SIFT": sift,
                "PolyPhen": polyphen,
                "Consequence": consequence,
            }
            # 这里consequenceDict，根据结果类型存放了每条数据
            for i in consequence.split(","):
                consequenceDict[i].append(obj)
            # 这里siftAndPolyPhenDict，根据sift和poly评分划分数据
            pattern = r'^\w+'
            match = re.search(pattern, sift)
            if match:
                siftCons = match.group(0)
                siftAndPolyPhenDict[siftCons].append(obj)

            pattern = r'^\w+'
            match = re.search(pattern, polyphen)
            if match:
                polyphenCons = match.group(0)
                siftAndPolyPhenDict[polyphenCons].append(obj)
            # print("consequenceDict",consequenceDict)
            # print("siftAndPolyPhenDict",siftAndPolyPhenDict)
    vepout.close()
    return consequenceDict, siftAndPolyPhenDict ,title,snpList


# 这里数据做：输入了哪些变异
def getSNPList(snpfile_path):
    snpfile = open(snpfile_path, 'r')
    snpList = []
    for line in snpfile:
        if not line.startswith('##') and not line.startswith('#'):
            # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            columns = line.strip().split(' ')
            if len(columns)==1:
                columns = line.strip().split("\t")
            columns = line.strip().split(' ')
            chrom = "chr"+columns[0]
            pos = columns[1]
            position = chrom+":"+pos
            snpid = columns[2]
            ref = columns[3]
            alt = columns[4]
            snpList.append([position, snpid, ref, alt])

    snpfile.close()
    return snpList


def assess_summary(P,snpfile):
    species=P.get_species()
    # 导入到pie里的数据：
    assessment_output_path = get_assessment_output_path(P, "")
    output_html=os.path.join(assessment_output_path,'assessment_summary.html')

    vep_output_path = get_assessment_output_path(P, "vep")
    vcffile = os.path.join(vep_output_path, os.path.basename(snpfile).replace(".vcf", "_vep.out.vcf"))
    
    consequenceDict, siftAndPolyPhenDict,objtitle,snpList = getConsequenceDictAndsiftAndPolyPhenDict(vcffile,snpfile)
    consequencePie = []
    for key in consequenceDict.keys():
        name = key
        value = len(consequenceDict[key])
        if value > 0:
            consequencePie.append({ "value":value,"name":name})


    siftcount=0
    for item in siftAndPolyPhenDict:
        siftcount+=len(siftAndPolyPhenDict[item])
        

    siftAndPolyPhenPie = []
    for key in siftAndPolyPhenDict.keys():
        name = key
        if name!="unknown":
            value = len(siftAndPolyPhenDict[key])
            if name=="tolerated_low_confidence":
                color="#82B0D2"
            elif name=="tolerated":
                color="#8ECFC9"
            elif name=="deleterious_low_confidence":
                color="#FFBE7A"
            elif name=="deleterious":
                color="#FA7F6F"
            elif name=="benign":
                color="#8ECFC9"
            elif name=="possibly_damaging":
                color="#FFBE7A"
            elif name=="probably_damaging":
                color="#FA7F6F"
            if value > 0:
                siftAndPolyPhenPie.append({ "value":value,"name":name,"itemStyle":{"color":color}})

    ########################################################################################################################
    #############################################         enformer        ###################################################
    if species=="homo_sapiens":
        targets_txt = f'{script_dir}/enformer/data/targets_human.txt'
        # json_file='my_list.json'
        snp_tracks=""
        # with open('my_list.json', 'r') as f:
        #     snp_tracks = json.load(f)



        # npydata='snp_tracks.npz'

        # data = np.load(npydata,allow_pickle=True)
        # json_data = {}
        # for key in data:
        #     json_data[key] = data[key].tolist()

        # json_data=json.dumps(json_data)
        # print(json_data)
        # # 将JSON格式的数据写入文件
        # with open(jsondata, 'w') as f:
        #     f.write('testDemo(')
        #     json.dump(json_data, f)
        #     f.write(')')
        # f.close()

    #############################################         enformer 下拉菜单       ###################################################
    
        selectHTML=""
        for snp in snpList:
            selectHTML+=f'<option value="{snp[0]}">{snp[0]}</option>'

        target = pd.read_table(targets_txt)
        targetList=[]
        optionHtml=""
        for i in range(len(target)):
            targetList.append([i,target['description'][i]])

        for tag in targetList:
            optionHtml+= f"<option value='{tag[0]}'>{tag[0]}.{tag[1]}</option>"
        init_snp=""+snpList[0][0].replace(":","_")+"/0"

    ########################################################################################################################

    # 创建HTML页面
    html_content = """
    <!DOCTYPE html>
    <html lang="en">

    <head>
        <meta charset="UTF-8">
        <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <script src="https://cdnjs.cloudflare.com/ajax/libs/echarts/5.4.1/echarts.min.js"></script>
        <!-- CSS -->
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@4.6.2/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-xOolHFLEh07PJGoPkLv1IbcEPTNtaed2xpHsD9ESMhqIYd0nLMwNLD69Npy4HI+N" crossorigin="anonymous">
        <script src="https://cdn.jsdelivr.net/npm/jquery@3.5.1/dist/jquery.slim.min.js" integrity="sha384-DfXdz2htPH0lsSSs5nCTpuj/zy4C+OGpamoFVy38MVBnE+IbbVYUew+OrCXaRkfj" crossorigin="anonymous"></script>
        <script src="https://cdn.jsdelivr.net/npm/bootstrap@4.6.2/dist/js/bootstrap.bundle.min.js" integrity="sha384-Fy6S3B9q64WdZWQUiU+q4/2Lc9npb8tCaSX9FK7E8HnRr0Jz8D6OP9dO5Vg3Q9ct" crossorigin="anonymous"></script>
        <title>Summary of Assessment Variant Effects</title>
        <style>
            thead{
                position:sticky;
                top:-1px;
                background-color:#f2f2f2;
                backface-visibility:hidden;
            }
            .table,.table-striped thead{
                border-bottom:1px solid #dee2e6;
            }
            .table,.table-striped thead th{
                border:0;
            }
        </style>
    </head>

    <body>
    <div style="margin:20px">

    <!-- #####################################
    ########## Consequence的card ########
    ##################################### -->
        <div class="card">
            <h5 class="card-header">Consequences Summary</h5>
            <div class="card-body">
                <div class="row">
                    <div class="col-4">
                        <div id="consequencePie" style="height:400px;"></div>
                    </div>
                    <div class="col-8">
                    <h5 id="consequenceTableTitle"></h5>
                    <div class="table table-striped" style="height:400px;overflow: auto">
                        <table id="consequenceTable" style="width: 100%;font-size:12px"></table>
                    </div>
                </div>
            </div>
        </div>
        </div>
    """
    if siftcount>0:
        html_content+="""
    <!-- ############################################
    ########## siftAndPolyPhen的card ########
    ############################################ -->
        <div class="card">
            <h5 class="card-header">SIFT & PolyPhen of Coding Variants Summary</h5>
            <div class="card-body">
                <div class="row">
                    <div class="col-4">
                        <div id="siftAndPolyPhenPie" style="height:400px;"></div>
                    </div>
                    <div class="col-8">
                    <h5 id="siftAndPolyPhenTableTitle"></h5>
                    <div class="table table-striped" style="height:400px;overflow: auto">
                        <table id="siftAndPolyPhenTable" style="width: 100%;font-size:12px"></table>
                    </div>
                </div>
            </div>
        </div>
        </div>

    """
    if species=="homo_sapiens":
        html_content+="""
        <!-- ############################################
        ########## enformerchart的card ########
        ############################################ -->
            <div class="card">
                <h5 class="card-header">Assessment of Variant Impact on Regulatory Activity in Tissues and Cells</h5>
                <div class="row" style="margin-top:15px">
                    <div class="col-1">
                        <p style="line-height:40px;margin-left:20px;font-size:18px;">SNP:</p>
                    </div>
                    <div class="col-2">
                        <select class="custom-select" id="select1" onchange='changeSelect(this.value)'>
                            """ +f"{str(selectHTML)}"+ """
                        </select>
                    </div>
                    <div class="col-2">
                        <p style="line-height:40px;margin-left:20px;font-size:18px;">Tissues and cell types:</p>
                    </div>
                    <div class="col-7">
                        <select class="custom-select" id="select2" onchange='changeSelect(this.value)'>
                            """ +f"{str(optionHtml)}"+ """
                        </select>
                    </div>
                </div>
                <hr style="margin:0">
                <div class="card-body">
                    <div class="row">
                        <div class="col">
                            <div id="enformerchart" style="height:800px;"></div>
                        </div>
                    </div>
                </div>
            </div>
            </div>
        """
    
    html_content+="""
    </div>
    <script type="text/javascript">
    // ####################  VEP data ####################
        let consequenceDict=""" +f"{str(consequenceDict)}"+ """
        let consequencePie=""" +f"{str(consequencePie)}"+ """
        let objtitle=""" +f"{str(objtitle)}"+ """
    // ###################################################
    """
    if siftcount>0:
        html_content+="""
        let siftAndPolyPhenDict=""" +f"{str(siftAndPolyPhenDict)}"+ """
        let siftAndPolyPhenPie=""" +f"{str(siftAndPolyPhenPie)}"+ """
    """
    if species=="homo_sapiens":
        html_content+="""
    // ####################  enformer data ####################
        // var snp_tracks=""" +f"{str(snp_tracks)}"+ """
        var snpList=""" +f"{str(snpList)}"+ """
        var SEQUENCE_LENGTH=393216

    """
    html_content+="""
        window.onload = function() {
            initConsequence()
    """
    if siftcount>0:
        html_content+="""
            initSiftAndPolyPhen()
        """
    if species=="homo_sapiens":
        html_content+="""
            initEnformer()
        """
    html_content+="""
        };
    """


    html_content+="""
        createConsequencePie()
    """
    if siftcount>0:
        html_content+="""
        createSiftAndPolyPhenPie()
        """
    html_content+="""
        function initConsequence(){
        // ##################################################
        // ############ 初始化Consequence的默认表格 ############
        // ##################################################
            // 创建标题
            let myconsequenceTableTitle = document.getElementById("consequenceTableTitle");
            myconsequenceTableTitle.textContent=consequencePie[0]["name"]

            // 创建表头（thead）元素
            let myconsequenceTable = document.getElementById("consequenceTable");
            
            let thead = myconsequenceTable.createTHead();
            let tr = thead.insertRow();
            for (let i = 0; i < objtitle.length; i++) {
                let temp =i
                let th = document.createElement("th");
                th.textContent = objtitle[temp];
                tr.appendChild(th);
            }
            // 创建新的表格主体
            let tbody = myconsequenceTable.createTBody();
            for (let i = 0; i<consequenceDict[consequencePie[0]["name"]].length;i++){
                let row = tbody.insertRow();
                for (let j = 0; j < objtitle.length; j++) {
                    let cell = row.insertCell();
                    let cellContent=consequenceDict[consequencePie[0]["name"]][i][objtitle[j]]
                    // 表格内容
                    if(cellContent==""){
                        cell.textContent="-"
                    }else{
                        cell.textContent=cellContent
                    }
                    // impact颜色之类的
                    if(objtitle[j]=="Impact"){
                        if(cellContent=="HIGH"){
                            cell.style.color="#FA7F6F"
                        }else if(cellContent=="MODERATE"){
                            cell.style.color="#FFBE7A"
                        }else if(cellContent=="LOW"){
                            cell.style.color="#8ECFC9"
                        }else if(cellContent=="MODIFIER"){
                            cell.style.color="#BEB8DC"
                        }
                    }else if(objtitle[j]=="SIFT"){
                        if(cellContent.startsWith("deleterious_low_confidence")){
                            cell.style.color="#FFBE7A"
                            cell.style.fontWeight="bold"
                        }else if(cellContent.startsWith("deleterious")){
                            cell.style.color="#FA7F6F"
                            cell.style.fontWeight="bold"
                        }
                    }else if(objtitle[j]=="PolyPhen"){
                        if(cellContent.startsWith("possibly_damaging")){
                            cell.style.color="#FFBE7A"
                            cell.style.fontWeight="bold"
                        }else if(cellContent.startsWith("probably_damaging")){
                            cell.style.color="#FA7F6F"
                            cell.style.fontWeight="bold"
                        }
                    }
                }
            }

        }

        function initSiftAndPolyPhen(){
        // ########################################################
        // ############ 初始化siftAndPolyPhen默认表格 ############
        // #######################################################
            // 创建标题
            let mysiftAndPolyPhenTableTitle = document.getElementById("siftAndPolyPhenTableTitle");
            mysiftAndPolyPhenTableTitle.textContent=siftAndPolyPhenPie[0]["name"]

            // 创建表头（thead）元素
            let mysiftAndPolyPhenTable = document.getElementById("siftAndPolyPhenTable");
        
            let thead = mysiftAndPolyPhenTable.createTHead();
            let tr = thead.insertRow();
            for (let i = 0; i < objtitle.length; i++) {
                let temp =i
                let th = document.createElement("th");
                th.textContent = objtitle[temp];
                tr.appendChild(th);
            }
            // 创建新的表格主体
            let tbody = mysiftAndPolyPhenTable.createTBody();
            for (let i = 0; i<siftAndPolyPhenDict[siftAndPolyPhenPie[0]["name"]].length;i++){
                let row = tbody.insertRow();
                for (let j = 0; j < objtitle.length; j++) {
                    let cell = row.insertCell();
                    let cellContent=siftAndPolyPhenDict[siftAndPolyPhenPie[0]["name"]][i][objtitle[j]]
                    // 表格内容
                    if(cellContent==""){
                        cell.textContent="-"
                    }else{
                        cell.textContent=cellContent
                    }
                    // impact颜色之类的
                    if(objtitle[j]=="Impact"){
                        if(cellContent=="HIGH"){
                            cell.style.color="#FA7F6F"
                        }else if(cellContent=="MODERATE"){
                            cell.style.color="#FFBE7A"
                        }else if(cellContent=="LOW"){
                            cell.style.color="#8ECFC9"
                        }else if(cellContent=="MODIFIER"){
                            cell.style.color="#BEB8DC"
                        }
                    }else if(objtitle[j]=="SIFT"){
                        if(cellContent.startsWith("deleterious_low_confidence")){
                            cell.style.color="#FFBE7A"
                            cell.style.fontWeight="bold"
                        }else if(cellContent.startsWith("deleterious")){
                            cell.style.color="#FA7F6F"
                            cell.style.fontWeight="bold"
                        }
                    }else if(objtitle[j]=="PolyPhen"){
                        if(cellContent.startsWith("possibly_damaging")){
                            cell.style.color="#FFBE7A"
                            cell.style.fontWeight="bold"
                        }else if(cellContent.startsWith("probably_damaging")){
                            cell.style.color="#FA7F6F"
                            cell.style.fontWeight="bold"
                        }
                    }
                }
            }

        }
        function initEnformer(){
        // ##############################################
        // ############ 初始化Enformer        ############
        // ##############################################
            let snp = document.getElementById('select1').value;
            let tag = document.getElementById('select2').value;
            let pos = Number(snp.split(':')[1])
            let start = pos - SEQUENCE_LENGTH/2 - 1
            let end = pos + SEQUENCE_LENGTH/2 + 1

            

            let intervalLinspace=linspace(start,end,tracks[0].length,pos)
            createEnformerChart(tracks,intervalLinspace)
        }

        function linspace(start, end, num ,center) {
            const tickInterval = Math.round((end - start) / (num - 1));
            // 计算总刻度数
            var totalTicks = num+1;

            // 计算中心点到起始值和终止值的刻度数
            var startTicks = num/2;
            var endTicks = num/2;

            // 初始化刻度数组
            var ticks = [];

            // 向左添加刻度值
            for (var i = 0; i < startTicks; i++) {
                ticks.push(center - (i + 1) * tickInterval);
            }
            
            ticks=ticks.reverse()
            // 向右添加刻度值
            for (var i = 0; i < endTicks; i++) {
                ticks.push(center + (i + 1) * tickInterval);
            }

            return ticks
        }
    // ########################################################

        function createConsequencePie(){    
        // #########################################
        // ############ Consequence饼图 ############
        // #########################################
            let myconsequencePie = echarts.init(document.getElementById('consequencePie'));
            consequencePie[0]["selected"]=true
            let consequencePieOption = {
                // 设置鼠标移上去之后显示的tip
                tooltip: {
                    position:"right",
                    formatter:function(obj){
                        let impact=consequenceDict[obj.data.name][0]["Impact"]
                        let color=""
                        let content=""
                        if(impact=="HIGH"){
                            color="FA7F6F"
                            content="The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function or triggering nonsense mediated decay."
                        }else if(impact=="MODERATE"){
                            color="#FFBE7A"
                            content="A non-disruptive variant that might change protein effectiveness."
                        }else if(impact=="LOW"){
                            color="#8ECFC9"
                            content="Assumed to be mostly harmless or unlikely to change protein behaviour."
                        }else if(impact=="MODIFIER"){
                            color="#BEB8DC"
                            content="Usually non-coding variants or variants affecting non-coding genes, where predictions are difficult or there is no evidence of impact."
                        }
                        
                        return(
                            '<div style="padding-bottom: 7px;margin-bottom: 7px;">' +
                            '<div style="font-size:14px">Consequence : <span style="font-size:15px;font-weight:bold">' + obj.data.name + '</sapn></div>' +
                            '<div style="font-size:14px">Count : <span style="font-size:15px;font-weight:bold">' + obj.data.value + '</sapn></div>' +
                            '<div style="font-size:14px">Impact : <span style="font-size:15px;font-weight:bold;color:'+color+'">' + impact + '</sapn></div>' +
                            '<div style="font-size:14px;width:300px;">Interpretation : <span style="white-space: normal;font-size:15px;">' + content + '</sapn></div>' +
                            '</div>'
                            )
                    }
                },
                // 设置饼图样式
                series: [
                    {
                        type: 'pie',
                        radius: ['20%','55%'],
                        data: consequencePie,
                        itemStyle: {
                            borderRadius: 10,
                            borderColor: '#fff',
                            borderWidth: 2
                        },
                        selectedOffset: 20, // 选中的扇区偏移量
                        selectedMode: "single",
                        avoidLabelOverlap: false,
                    }
                ]
            };
            myconsequencePie.setOption(consequencePieOption);
            // 设置点击一个扇形之后改变右侧表格的内容
            myconsequencePie.on('click', function (params) {
                // 修改标题
                let myconsequenceTableTitle = document.getElementById("consequenceTableTitle");
                myconsequenceTableTitle.textContent=params["name"]

                // 清空table
                let myconsequenceTable = document.getElementById("consequenceTable");
                while (myconsequenceTable.rows.length > 0) {
                    myconsequenceTable.deleteRow(0);
                }

                // 创建表头（thead）元素
                let thead = myconsequenceTable.createTHead();
                let tr = thead.insertRow();
                for (let i = 0; i < objtitle.length; i++) {
                    let temp =i
                    let th = document.createElement("th");
                    th.textContent = objtitle[temp];
                    tr.appendChild(th);
                }
                // 创建新的表格主体
                let tbody = myconsequenceTable.createTBody();
                for (let i = 0; i<consequenceDict[params["name"]].length;i++){
                    let row = tbody.insertRow();
                    for (let j = 0; j < objtitle.length; j++) {
                        let cell = row.insertCell();
                        let cellContent=consequenceDict[params["name"]][i][objtitle[j]]
                        // 表格内容
                        if(cellContent==""){
                            cell.textContent="-"
                        }else{
                            cell.textContent=cellContent
                        }
                        // impact颜色之类的
                        if(objtitle[j]=="Impact"){
                            if(cellContent=="HIGH"){
                                cell.style.color="#FA7F6F"
                            }else if(cellContent=="MODERATE"){
                                cell.style.color="#FFBE7A"
                            }else if(cellContent=="LOW"){
                                cell.style.color="#8ECFC9"
                            }else if(cellContent=="MODIFIER"){
                                cell.style.color="#BEB8DC"
                            }
                        }else if(objtitle[j]=="SIFT"){
                            if(cellContent.startsWith("deleterious_low_confidence")){
                                cell.style.color="#FFBE7A"
                                cell.style.fontWeight="bold"
                            }else if(cellContent.startsWith("deleterious")){
                                cell.style.color="#FA7F6F"
                                cell.style.fontWeight="bold"
                            }
                        }else if(objtitle[j]=="PolyPhen"){
                            if(cellContent.startsWith("possibly_damaging")){
                                cell.style.color="#FFBE7A"
                                cell.style.fontWeight="bold"
                            }else if(cellContent.startsWith("probably_damaging")){
                                cell.style.color="#FA7F6F"
                                cell.style.fontWeight="bold"
                            }
                        }
                    }
                }
            });
        }




        function createSiftAndPolyPhenPie(){
        // ################################################
        // ############ siftAndPolyPhenDict饼图 ############
        // ################################################
            let mysiftAndPolyPhenPie = echarts.init(document.getElementById('siftAndPolyPhenPie'));
            siftAndPolyPhenPie[0]["selected"]=true
            let siftAndPolyPhenPieOption = {
                // 设置鼠标移上去之后显示的tip
                tooltip: {
                    position:"right",
                    formatter:function(obj){
                        let color=""
                        let content=""

                        let siftList=["tolerated_low_confidence","tolerated","deleterious_low_confidence","deleterious"]
                        let type=""
                        let socre=""
                        if(siftList.includes(obj.data.name)){
                            type="SIFT"
                            socre=siftAndPolyPhenDict[obj.data.name][0]["SIFT"]
                        }else{
                            type="PolyPhen"
                            socre=siftAndPolyPhenDict[obj.data.name][0]["PolyPhen"]
                        }

                        if(socre.startsWith("tolerated_low_confidence")){
                            color="#82B0D2"
                            content="The variant is predicted to have no significant impact on protein function, but the confidence in this prediction is low due to limitations in sequence alignment or protein structure information."
                        }else if(socre.startsWith("tolerated")){
                            color="#8ECFC9"
                            content="The variant is predicted to have no significant impact on protein function. These variants are usually found in highly variable regions and are less likely to be disease-causing."
                        }else if(socre.startsWith("deleterious_low_confidence")){
                            color="#FFBE7A"
                            content="The variant is predicted to have a significant impact on protein function, but the confidence in this prediction is low due to limitations in sequence alignment or protein structure information."
                        }else if(socre.startsWith("deleterious")){
                            color="#FA7F6F"
                            content="The variant is predicted to have a significant impact on protein function, potentially causing changes to protein structure or function and leading to disease or other biological effects."
                        }else if(socre.startsWith("benign")){
                            color="#8ECFC9"
                            content="The substitution is predicted to have a low probability of causing an impact on protein function, and is usually found in highly variable regions where substitutions are less likely to cause disease."
                        }else if(socre.startsWith("possibly_damaging")){
                            color="#FFBE7A"
                            content="The substitution is predicted to have a moderate probability of causing an impact on protein function, but the effect may not be as significant or severe."
                        }else if(socre.startsWith("probably_damaging")){
                            color="#FA7F6F"
                            content="The substitution is predicted to have a high probability of causing a significant impact on protein function, potentially leading to disease or other biological effects."
                        }

                        return(
                            '<div style="padding-bottom: 7px;margin-bottom: 7px;">' +
                            '<div style="font-size:14px">Socre Type : <span style="font-size:15px;font-weight:bold">' + type + '</sapn></div>' +
                            '<div style="font-size:14px">Impact : <span style="font-size:15px;font-weight:bold;color:'+color+'">' + obj.data.name + '</sapn></div>' +
                            '<div style="font-size:14px">Count : <span style="font-size:15px;font-weight:bold">' + obj.data.value + '</sapn></div>' +
                            '<div style="font-size:14px;width:300px;">Interpretation : <span style="white-space: normal;font-size:15px;">' + content + '</sapn></div>' +
                            '</div>'
                            )
                    }
                },
                // 设置饼图样式
                series: [
                    {
                        type: 'pie',
                        radius: ['20%','55%'],
                        data: siftAndPolyPhenPie,
                        itemStyle: {
                            borderRadius: 10,
                            borderColor: '#fff',
                            borderWidth: 2
                        },
                        selectedOffset: 20, // 选中的扇区偏移量
                        selectedMode: "single",
                        avoidLabelOverlap: false,
                    }
                ]
            };
            mysiftAndPolyPhenPie.setOption(siftAndPolyPhenPieOption);
            // 设置点击一个扇形之后改变右侧表格的内容
            mysiftAndPolyPhenPie.on('click', function (params) {
                // 修改标题
                let mysiftAndPolyPhenTableTitle = document.getElementById("siftAndPolyPhenTableTitle");
                mysiftAndPolyPhenTableTitle.textContent=params["name"]

                // 清空table
                let mysiftAndPolyPhenTable = document.getElementById("siftAndPolyPhenTable");
                while (mysiftAndPolyPhenTable.rows.length > 0) {
                    mysiftAndPolyPhenTable.deleteRow(0);
                }

                // 创建表头（thead）元素
                let thead = mysiftAndPolyPhenTable.createTHead();
                let tr = thead.insertRow();
                for (let i = 0; i < objtitle.length; i++) {
                    let temp =i
                    let th = document.createElement("th");
                    th.textContent = objtitle[temp];
                    tr.appendChild(th);
                }
                // 创建新的表格主体
                let tbody = mysiftAndPolyPhenTable.createTBody();
                for (let i = 0; i<siftAndPolyPhenDict[params["name"]].length;i++){
                    let row = tbody.insertRow();
                    for (let j = 0; j < objtitle.length; j++) {
                        let cell = row.insertCell();
                        let cellContent=siftAndPolyPhenDict[params["name"]][i][objtitle[j]]
                        // 表格内容
                        if(cellContent==""){
                            cell.textContent="-"
                        }else{
                            cell.textContent=cellContent
                        }
                        // impact颜色之类的
                        if(objtitle[j]=="Impact"){
                            if(cellContent=="HIGH"){
                                cell.style.color="#FA7F6F"
                            }else if(cellContent=="MODERATE"){
                                cell.style.color="#FFBE7A"
                            }else if(cellContent=="LOW"){
                                cell.style.color="#8ECFC9"
                            }else if(cellContent=="MODIFIER"){
                                cell.style.color="#BEB8DC"
                            }
                        }else if(objtitle[j]=="SIFT"){
                            if(cellContent.startsWith("deleterious_low_confidence")){
                                cell.style.color="#FFBE7A"
                                cell.style.fontWeight="bold"
                            }else if(cellContent.startsWith("deleterious")){
                                cell.style.color="#FA7F6F"
                                cell.style.fontWeight="bold"
                            }
                        }else if(objtitle[j]=="PolyPhen"){
                            if(cellContent.startsWith("possibly_damaging")){
                                cell.style.color="#FFBE7A"
                                cell.style.fontWeight="bold"
                            }else if(cellContent.startsWith("probably_damaging")){
                                cell.style.color="#FA7F6F"
                                cell.style.fontWeight="bold"
                            }
                        }
                    }
                }
            });
        }

        function createEnformerChart(tracks,intervalLinspace) {
        // ################################################
        // ############    myEnformerChart     ############
        // ################################################
            let myEnformerChart = echarts.init(document.getElementById('enformerchart'));
            let ref=tracks[0]
            let alt=tracks[1]
            let alt_ref=tracks[2]
            let enformerchartOption = {
                tooltip: {
                    trigger: 'axis',
                    position: function (pt) {
                    return [pt[0], '10%'];
                    }
                },
                title: [{
                    left: 'center',
                    top: '0%',
                    text: 'Variant effect assessment (ref & alt)'
                },{
                    left: 'center',
                    top: '49%',
                    text: 'Effect of variant on regulatory activity (alt - ref)'
                }
                ],
                legend: [
                    {
                        top: '4%',
                        data: ['ref','alt'],
                        show: true // 显示上图的legend
                    },
                    {
                        top: '53%',
                        data: [
                            {
                                name:'alt - ref'
                            }
                        ],
                        show: true // 显示下图的legend
                    }
                ],
                
                grid: [
                    {
                        left: '3%',
                        right: '3%',
                        top: '7%',
                        height: '35%',
                        containLabel: true,
                    },
                    {
                        left: '1%',
                        right: '3%',
                        top: '57%',
                        height: '35%',
                        containLabel: true,
                    }
                ],
                xAxis: [
                    {   
                        name: 'pos',
                        type: 'category',
                        boundaryGap: false,
                        data:intervalLinspace,
                        axisLabel: {
                            interval: intervalLinspace.length/7,
                            color: '#999',
                            
                            formatter: function(value, index) {
                                if (index == 452) {
                                    return '{value|' + value + '}';
                                } else {
                                    return value;
                                }
                            },
                            rich: {
                                value: {
                                    color: 'red'
                                }
                            }
                        },
                    },
                    {
                        name: 'pos',
                        type: 'category',
                        boundaryGap: false,
                        data:intervalLinspace,
                        axisLabel: {
                            interval: intervalLinspace.length/7
                        },
                        gridIndex: 1
                    },
                ],
                yAxis: [
                    {
                        name: 'score',
                        type: 'value',
                    },
                    {
                        name:'score',
                        type: 'value',
                        gridIndex: 1
                    }
                ],
                series: [
                    {
                        name:'ref',
                        type: 'line',
                        tack: 'Total',
                        data: ref,
                        markPoint: {
                            data: [
                                {   
                                    xAxis: intervalLinspace.length/2,
                                    yAxis: 0,
                                    symbolRotate:180,
                                    itemStyle: { // 标记点样式
                                        color:"#FA7F6F",
                                        
                                    },
                                    
                                    label: { // 标记点标签样式
                                        show: true, // 显示标签
                                        position: 'insideBottom', // 标签位置
                                        formatter: function(params) { // 自定义标签内容
                                            return 'Variant';
                                        },
                                        fontWeight: 'bold', // 标签文字加粗
                                        color: 'white', // 标签文字颜色
                                        backgroundColor: '#FA7F6F', // 标签背景颜色
                                        padding: [5, 10], // 标签内边距
                                        borderRadius: 5 // 标签圆角
                                    }
                                    
                                }
                            ]
                        },
                    },
                    {
                        name:'alt',
                        type: 'line',
                        tack: 'Total',
                        data: alt,
                    },
                    {
                        name:'alt - ref',
                        type: 'line',
                        data: alt_ref,
                        xAxisIndex: 1, // 指定该系列所属的x轴为第一个x轴
                        yAxisIndex: 1, // 指定该系列所属的y轴为第一个y轴
                    }
                ],
                toolbox: {
                    feature: {
                    saveAsImage: {}
                    }
                },
                dataZoom: [
                    {
                        show: true,
                        xAxisIndex: [0, 1]
                    },
                ],

            };


            myEnformerChart.setOption(enformerchartOption);
        }

        function changeSelect(){
        // ################################################
        // ###############   changeSelect  ###############
        // ################################################
            let snp = document.getElementById('select1').value;
            let tag = document.getElementById('select2').value;

            var script = document.createElement("script");
            
            script.type="text/javascript"
            script.src = "./enformer/" + snp.replace(":","_") + "/"+tag+".js"; // JavaScript 文件的路径和名称
            script.onload = function() {
                let pos = Number(snp.split(':')[1])
                let start = pos - SEQUENCE_LENGTH/2 - 1
                let end = pos + SEQUENCE_LENGTH/2 + 1

                
                let intervalLinspace=linspace(start,end,tracks[0].length,pos)
                createEnformerChart(tracks,intervalLinspace)
            };
            document.body.appendChild(script); // 添加到 body 元素中
        }
    </script>
    """
    if species=="homo_sapiens":
        html_content+="""
    <script type="text/javascript" src="./enformer/"""+init_snp+""".js"></script>
    """
    html_content+="""    
    </body>
    </html>
    """

    # 将HTML页面写入文件
    with open(output_html, 'w') as f:
        f.write(html_content)
    f.close()