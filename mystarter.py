#!/usr/bin/env python
from myparsing import *
from myutils import *

argvals = sys.argv[1:]
# print(argvals)
#################### Setting variables ####################
P = Parser(argvals)
output_dir = get_output_dir(P)
#给sratool添加权限
sratool_add_chmod()
#################### 流程判断 ####################
step = P.get_step()
if step:
    os.makedirs(output_dir, exist_ok=True)
# print(step)
#################### 第一部分 WGS ####################
flag = 0
sra = P.get_sra()
sralist = P.get_sralist()
if sra or sralist:
    flag = 1
sradir = P.get_sradir()
if sradir:
    flag = 2
rawfastqdir = P.get_rawfastqdir()
if rawfastqdir:
    flag = 3
fastqdir = P.get_fastqdir()
if fastqdir:
    flag = 4
cleanfastqdir = P.get_cleanfastqdir()
if cleanfastqdir:
    flag = 5
bamdir= P.get_bamdir()
if bamdir:
    flag = 6
processedbamdir= P.get_processedbamdir()
if processedbamdir:
    flag = 7
gvcfdir= P.get_gvcfdir()
if gvcfdir:
    flag = 8
vcffile = P.get_vcffile()
if vcffile:
    flag = 9


########################################################################################################################
########################################################################################################################
########################################################################################################################
#################### 1、  SRA数据下载 ####################
if step in ["downloadsra", "wgsall"]:
    download_list = []
    if sra:
        download_list.extend(sra)
    if sralist:
        f = open(sralist, 'r')
        sra_list = f.read().split("\n")
        download_list.extend(sra_list)

    do_download_sra(P, download_list)
#################### 对文件进行校验
    check_sra(P)
    report_err_sra(P)
#################### 2、  sra转fastq ####################
if step in ["sratofastq", "wgsall"]:
    sra_file_path_list = []
    if flag < 2:
        sra_path = get_wgs_output_path(P, "sra")
        tmp_sra_file_path_list = get_sra_file_path(sra_path)
        sra_file_path_list.extend(tmp_sra_file_path_list)
    if flag == 2:
        tmp_sra_file_path_list = get_sra_file_path(sradir)
        sra_file_path_list.extend(tmp_sra_file_path_list)

    do_sra_to_fastq(P, sra_file_path_list)
    zip_fastq(P)
#################### 3、  fastq质控 ####################
if step in ["readsqc", "wgsall"]:
    single_fastq_file_path_list = []
    paired_fastq_file_path_list = []
    
    if flag < 3:
        raw_path = get_wgs_output_path(P, "raw")
        tmp_single, tmp_paired = get_fastq_file_path(raw_path)
        single_fastq_file_path_list.extend(tmp_single)
        paired_fastq_file_path_list.extend(tmp_paired)
    if flag == 3:
        tmp_single, tmp_paired = get_fastq_file_path(rawfastqdir)
        single_fastq_file_path_list.extend(tmp_single)
        paired_fastq_file_path_list.extend(tmp_paired)
    do_qc(P, single_fastq_file_path_list, paired_fastq_file_path_list)

#################### 4、  fastqc+multiqc质量评估 ####################

if step in ["qualityevaluation", "wgsall"]:
    qc_file_path_list = []
    if flag < 4:
        clean_output_path = get_wgs_output_path(P, "clean")
        tmp_qc_file_path_list = get_qc_file_path(clean_output_path)
        qc_file_path_list.extend(tmp_qc_file_path_list)
    if flag == 4:
        tmp_qc_file_path_list = get_qc_file_path(fastqdir)
        qc_file_path_list.extend(tmp_qc_file_path_list)
    do_fastqc(P, qc_file_path_list)
    do_multiqc(P)


#################### 5、  建立基因组索引 ####################
if step in ["downloadref", "wgsall"]:
    download_ref(P)
#################### 6、  bwa比对参考基因组 ####################
if step in ["align", "wgsall"]:
    single_fastq_file_path_list = []
    paired_fastq_file_path_list = []
    if flag < 5:
        clean_path = get_wgs_output_path(P, "clean")
        tmp_single, tmp_paired = get_fastq_file_path(clean_path)
        single_fastq_file_path_list.extend(tmp_single)
        paired_fastq_file_path_list.extend(tmp_paired)
    if flag == 5:
        tmp_single, tmp_paired = get_fastq_file_path(cleanfastqdir)
        single_fastq_file_path_list.extend(tmp_single)
        paired_fastq_file_path_list.extend(tmp_paired)
    do_align_bwa(P, single_fastq_file_path_list, paired_fastq_file_path_list)

#################### 7、  处理bam文件 ####################
if step in ["dealbam", "wgsall"]:
    if flag < 6:
        bam_file_path = get_wgs_output_path(P, "align")
    if flag == 6:
        bam_file_path = bamdir
    delPCR=P.get_delPCR()

    sort_bam(P,bam_file_path)
    if delPCR:
        mark_pcr_repeat_picard(P)
    index_for_bam(P)


#################### 8、  变异检测 ####################
if step in ["detect", "wgsall"]:
    if flag < 7:
        bam_file_path = get_wgs_output_path(P, "processed")
    if flag == 7:
        bam_file_path = processedbamdir
    
    do_gvcf(P,bam_file_path)

#################### 9、  jointgenotype ####################
if step in ["jointgenotype", "wgsall"]:
    if flag < 8:
        gvcf_file_path = get_wgs_output_path(P, "gvcf")
    if flag == 8:
        gvcf_file_path = gvcfdir

    do_jointgenotype(P,gvcf_file_path)

#################### 10、  vcf质控 ####################
if step in ["vcfqc", "wgsall"]:
    vcf_file_path_list = []
    if flag < 9:
        vcf_output_path = get_wgs_output_path(P, "vcf")
        vcf_file = os.path.join(vcf_output_path, "genotype.vcf.gz")
    if flag == 9:
        vcf_file = vcffile
    do_vcf_qc(P, vcf_file)
########################################################################################################################
########################################################################################################################
########################################################################################################################
#################### 第2部分 GWAS ####################
genotypefile = P.get_genotypefile()
if genotypefile:
    flag = 11
bfiledir = P.get_bfiledir()
if bfiledir:
    flag = 12
cleanbfiledir = P.get_cleanbfiledir()
if cleanbfiledir:
    flag = 13
# #################### 2.0、  基因型填充 ####################
if step in ["impute", "gwasall"]:
    if flag == 11:
        vcf_file = genotypefile
    impute_vcf(P, vcf_file)
    
# #################### 2.1、  vcf转ped、map\bim、bed、fam ####################
if step in ["transvcf", "gwasall"]:
    if flag < 11:
        vcf_output_path = get_gwas_output_path(P, "transvcf")
        vcf_file = os.path.join(vcf_output_path, "genotype.vcf.gz")
    if flag == 11:
        vcf_file = genotypefile
    trans_vcf(P, vcf_file)
    merge_pheno_fam(P)




# #################### 2.2、  gwas质控 ####################
if step in ["gwasqc", "gwasall"]:
    species = P.get_species()
    if flag < 12:
        gwasqc_path = get_gwas_output_path(P, "gwasqc")
    if flag == 12:
        gwasqc_path = bfiledir

    do_gwasqc(P,gwasqc_path)

# #################### 2.3、  pca分析 ####################
if step in ["pca", "gwasall"]:
    if flag < 13:
        clean_file_path = get_gwas_output_path(P, "gwasqc")
    if flag == 13:
        clean_file_path = cleanbfiledir

    do_population(P,clean_file_path)


# #################### 2.4、  kinship ####################
if step in ["kinship", "gwasall"]:
    if flag < 13:
        clean_file_path = get_gwas_output_path(P, "gwasqc")
    if flag == 13:
        clean_file_path = cleanbfiledir

    do_bed_to_grm(P, clean_file_path)
    arrange_kinship(P)
    draw_kinship(P)

# #################### 2.5、  关联分析 ####################
if step in ["association", "gwasall"]:
    if flag < 13:
        clean_file_path = get_gwas_output_path(P, "gwasqc")
    if flag == 13:
        clean_file_path = cleanbfiledir
    if P.get_lm():
        do_associate_lm(P, clean_file_path)
    if P.get_lmm():
        do_associate_lmm(P,clean_file_path)

        
if step in ["selectsnp", "gwasall"]:
    do_selectsnp(P)


########################################################################################################################
########################################################################################################################
########################################################################################################################
#################### 第3部分 预测 ####################
snpfile = P.get_snpfile()
species=P.get_species()
    
if snpfile:
    flag = 20

# #################### 3.1、  enformer ####################
# if step in ["enformer"]:
#     if species=="homo_sapiens":
#         if flag < 20:
#             lm_path = get_gwas_output_path(P, os.path.join("association", "lm"))
#             lmm_path = get_gwas_output_path(P, os.path.join("association", "lmm"))
#             lm_snp = os.path.join(lm_path, "variant_lm.vcf")
#             lmm_snp = os.path.join(lmm_path, "variant_lmm.vcf")
#             if os.path.exists(lm_snp):
#                 enformer_predict(P, lm_snp)
#             if os.path.exists(lmm_snp):
#                 enformer_predict(P, lmm_snp)
#         if flag == 20:
#             enformer_predict(P, snpfile)
#     else:
#         print("Only homo_sapiens can use the enformer")

# # #################### 3.2、  vep ####################
# if step in ["vep"]:
#     if flag < 20:
#         lm_path = get_gwas_output_path(P, os.path.join("association", "lm"))
#         lmm_path = get_gwas_output_path(P, os.path.join("association", "lmm"))
#         lm_snp = os.path.join(lm_path, "variant_lm.vcf")
#         lmm_snp = os.path.join(lmm_path, "variant_lmm.vcf")
#         if os.path.exists(lm_snp):
#             vep_predict(P, lm_snp)
#         if os.path.exists(lmm_snp):
#             vep_predict(P, lmm_snp)

#     if flag == 20:
#         vep_predict(P, snpfile)
# #################### 3.3、  整合结果 ####################
if step in ["assess"]:
    if flag == 20:
        vep_predict(P, snpfile)
        if species=="homo_sapiens":
            enformer_predict(P, snpfile)

        assess_summary(P,snpfile)
