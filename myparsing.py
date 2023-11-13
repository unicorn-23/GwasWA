#!/usr/bin/env python
import argparse
import os

def set_version_number():
    return "1.0"

def getArgs(argv=None):
    default_output_dir = os.getcwd()
    parser = argparse.ArgumentParser(
        prog="gwaswa", #指定程序的名称
        formatter_class=argparse.RawTextHelpFormatter,  # 用于指定帮助文档格式
        description='''
        1、python gwaswa --step wgsall --sralist example/wgs/srr_list.txt --accession GCF_000157115.2
        2、python gwaswa --step gwasall --genotypefile example/gwas/example.vcf.gz --phenotypefile example/gwas/pheno.txt --checksex --removerelationship

        ''',
        # epilog="For a detailed description of all options, please refer to the manual."# 帮助文档尾部的信息
        )
    version = set_version_number()
    parser.add_argument(
        '--version', action='version', version='%(prog)s Version {}'.format(version)
    )
    parser.add_argument(
        "-o", "--output", metavar="<path>", type=str, default=default_output_dir, 
        help="change the output directory \ndefault: %(default)s\ndirectory will be created if non-existent"
    )
    parser.add_argument(
        "--nosave", action="store_true",
        help="是否存储中间文件"
    )
    parser.add_argument(
        "--nThrds", metavar="<str>", type=str,
        help='''多线程数量'''
    )
    parser.add_argument(
        "--nMem", metavar="<str>", type=str,
        help='''最大占用内存.'''
    )
    #################### 选择流程 ####################
    parser.add_argument(
        "--step", metavar="<str>", type=str,
        choices=[

            "wgsall",
            "downloadsra",
            "sratofastq",
            "readsqc",
            "qualityevaluation",
            "downloadref",
            "align",
            "dealbam",
            "detect",
            "jointgenotype",
            "vcfqc",

            # "gwasall",
            "transvcf",
            "impute",
            "gwasqc",
            "pca",
            "kinship",
            "association",
            "selectsnp",
            
            "assess",
            # "enformer",
            # "vep",
            ],
        help='''
        download，只运行下载sra
        sratofastqonly，只运行sra转换为fastq格式
        readsqc,只运行fastq格式的质控
        multiqc,只运行fastq格式的质量评估
        downloadref,只运行参考基因组索引
        align,只运行参考基因组索引
        dealbam,只对bwa比对参考基因组生成的bam文件进行处理
        realignonly，只对gvcf文件进行合并重比对
        vcfqc，只对vcf进行质控


        transvcf，只运行vcf转bfile
        gwasqc，只运行SNP、sample level质量控制
        pca，只运行主成分分析
        kinship，只运行亲缘矩阵
        association，只运行关联分析

        assess，编译预测全流程
        enformer，只运行enformer
        vep，只运行vep
        '''
    )

########################################################################################################################
########################################################################################################################
########################################################################################################################
    #################### 第一部分 WGS ####################
#################### 1、  SRA数据下载 ####################
    parser.add_argument(
        "--sra", metavar="<str>", type=str, nargs="+", help="例如:--sra SRR3418005 SRR3418019"
    )
    parser.add_argument(
        "--sralist", metavar="<filename>", type=str,
        help='''例如:--sralist path/srr_list.txt
        需要一个txt文件srr_list.txt,每行一个SRRid
        '''
    )
    
#################### 2、  sra转fastq ####################
    parser.add_argument(
        "--sradir", metavar="<path>", type=str,
        help='''例如:--sradir Output/wgs/sra
        需要sra文件目录
        '''
    )
#################### 3、  fastq质控 ####################
    parser.add_argument(
        "--rawfastqdir", metavar="<path>", type=str,
        help='''例如:--rawfastqdir path/fastq/
        需要fastq文件目录
        '''
    )
    parser.add_argument(
        "--quality",metavar="<int>",type=int,default=20,
        help='''设定Phred quality score阈值,从reads中修剪低质量末端。
        默认 Phred 分数:%(default)s
        ''')

    parser.add_argument(
        "--phred", metavar="<str>", type=str,default="phred33",
        choices=["phred33","phred64"],
        help='''phred33、phred64
        '''
    )
    parser.add_argument(
        "--length",metavar="<int>",type=int,default=20,
        help='''丢弃由于质量或适配器修整而变得短于length<INT>的读取。值0有效地禁用此行为。
        默认:%(default)sbp
        ''')
    parser.add_argument(
        "--stringency",metavar="<int>",type=int,default=1,
        help='''设定与接头重叠的序列,设定可以忍受的前后adapter重叠的碱基数,可以适度放宽,因为后一个adapter几乎不可能被测序仪读到。
        默认:%(default)s（非常苛刻）。
        ''')
    parser.add_argument(
        "--error",metavar="<float>",type=float,default=0.1,
        help='''最大允许错误率（错误数除以匹配区域的长度）
        默认:%(default)s
        ''')
#################### 4、  qc质量评估 ####################
    parser.add_argument(
        "--fastqdir", metavar="<path>", type=str,
        help='''例如:--fastqdir path/clean_fastq/
        需要fastq文件目录
        '''
    )

#################### 5、  建立基因组索引 ####################
    parser.add_argument(
        "--accession", metavar="<str>", type=str,
        help='''NCBI Reference sequence,e.g.,
        --accession GCF_000001405.40
        '''
    )
    parser.add_argument(
        "--taxon", metavar="<str>", type=str,
        help='''NCBI Taxonomy ID 或分类名称
        '''
    )
    parser.add_argument(
        "--refgenome", metavar="<filename>", type=str,
        help='''Local reference sequence path,e.g.,
        --refgenome .../ref.fa
        '''
    )
    parser.add_argument(
        "--indexalgorithm", metavar="<str>", type=str, default="bwtsw",
        choices=["is", "bwtsw"],
        help='''Algorithm for constructing index
        Default: %(default)s
        Available metrics: %(choices)s
        '''
    )
    
#################### 6、  比对参考基因组-bwa ####################
    parser.add_argument(
        "-afa","--alignalgorithm", metavar="<str>", type=str, default="mem",
        choices=["mem","bwasw","backtrack"],
        help='''Algorithm for alignment
        mem(首推):在reads长度在70bp-1Mbp范围时,推荐本算法。支持long-reads,split alignment。
        bwasw:在reads具有频繁的gap时,比对更敏感,推荐本算法。reads长度一般为70bp-1Mbp,支持long-reads,split alignment。
        backtrack:reads长度<70bp时,推荐本算法,建议输入reads长度 < 100bp。
        Default: %(default)s
        Available metrics: %(choices)s
        '''
    )

    parser.add_argument(
        "--cleanfastqdir", metavar="<path>", type=str,
        help='''存放需要比对的文件的目录'''
    )
    # parser.add_argument(
    #     "-rgpl", metavar="<str>", type=str,default="UNKNOWN",
    #     choices=["ILLUMINA","SLX","SOLEXA","SOLID","454","LS454","COMPLETE","PACBIO","IONTORRENT","CAPILLARY","HELICOS","UNKNOWN"],
    #     help='''
    #         测序平台
    #     '''
    # )
#################### mem参数
    parser.add_argument(
        "--minSeedLen", metavar="<str>", type=str,
        help='''Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20.'''
    )
    parser.add_argument(
        "--bandWidth", metavar="<str>", type=str,
        help='''Band width. Essentially, gaps longer than INT will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option.'''
    )
    parser.add_argument(
        "--zDropoff", metavar="<str>", type=str,
        help='''Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment.'''
    )
    parser.add_argument(
        "--seedSplitRatio", metavar="<str>", type=str,
        help='''Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy.'''
    )
    parser.add_argument(
        "--maxOcc", metavar="<str>", type=str,
        help='''Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter.'''
    )
    parser.add_argument(
        "--matchScore", metavar="<str>", type=str,
        help='''Matching score.'''
    )
    parser.add_argument(
        "--mmPenalty", metavar="<str>", type=str,
        help='''Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}'''
    )
    parser.add_argument(
        "--gapOpenPen", metavar="<str>", type=str,
        help='''Gap open penalty.'''
    )
    parser.add_argument(
        "--gapExtPen", metavar="<str>", type=str,
        help='''Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap)'''
    )
    parser.add_argument(
        "--clipPen", metavar="<str>", type=str,
        help='''Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted.'''
    )
    parser.add_argument(
        "--unpairPen", metavar="<str>", type=str,
        help='''Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing.'''
    )

#################### aln参数
# https://bio-bwa.sourceforge.net/bwa.shtml
    parser.add_argument(
        "--maxDiff", metavar="<str>", type=str,
        help='''
            Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths.
        '''
    )
    parser.add_argument(
        "--maxGapO", metavar="<str>", type=str,
        help='''
            Maximum number of gap opens.
        '''
    )
    parser.add_argument(
        "--maxGapE", metavar="<str>", type=str,
        help='''
            Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps)
        '''
    )
    parser.add_argument(
        "--nDelTail", metavar="<str>", type=str,
        help='''
            Disallow a long deletion within INT bp towards the 3’-end
        '''
    )

    parser.add_argument(
        "--nIndelEnd", metavar="<str>", type=str,
        help='''
            Disallow an indel within INT bp towards the ends
        '''
    )
    parser.add_argument(
        "--maxSeedDiff", metavar="<str>", type=str,
        help='''
            Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for ‘-k 2’.
        '''
    )

    parser.add_argument(
        "--seedLen", metavar="<str>", type=str,
        help='''
            Maximum edit distance in the seed.
        '''
    )

    parser.add_argument(
        "--misMsc", metavar="<str>", type=str,
        help='''
            Mismatch penalty. BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc).
        '''
    )
    parser.add_argument(
        "--gapOsc", metavar="<str>", type=str,
        help='''
            Gap open penalty.
        '''
    )
    parser.add_argument(
        "--gapEsc", metavar="<str>", type=str,
        help='''
            Gap extension penalty.
        '''
    )
    parser.add_argument(
        "--trimQual", metavar="<str>", type=str,
        help='''
            	Parameter for read trimming.
        '''
    )

#################### samse
    parser.add_argument(
        "--maxHitPaired", metavar="<str>", type=str,
        help='''
            	Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written.
        '''
    )
#################### sampe
    parser.add_argument(
        "--maxInsSize", metavar="<str>", type=str,
        help='''
            	Maximum insert size for a read pair to be considered being mapped properly. Since 0.4.5, this option is only used when there are not enough good alignment to infer the distribution of insert sizes.
        '''
    )

    parser.add_argument(
        "--maxHitDis", metavar="<str>", type=str,
        help='''
            	Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written.
        '''
    )
#################### bwasw参数
    parser.add_argument(
        "--mmPen", metavar="<str>", type=str,
        help='''
            Mismatch penalty
        '''
    )
    parser.add_argument(
        "--thres", metavar="<str>", type=str,
        help='''
            Minimum score threshold divided by a
        '''
    )
    parser.add_argument(
        "--hspIntv", metavar="<str>", type=str,
        help='''
            Coefficient for threshold adjustment according to query length. Given an l-long query, the threshold for a hit to be retained is a*max{T,c*log(l)}.[5.5]
        '''
    )
    parser.add_argument(
        "--zBest", metavar="<str>", type=str,
        help='''
            Z-best heuristics. Higher -z increases accuracy at the cost of speed.[1]
        '''
    )
    parser.add_argument(
        "--nHspRev", metavar="<str>", type=str,
        help='''
            Maximum SA interval size for initiating a seed. Higher -s increases accuracy at the cost of speed. [3]
        '''
    )
    parser.add_argument(
        "--thresCoef", metavar="<str>", type=str,
        help='''
            Coefficient for threshold adjustment according to query length. Given an l-long query, the threshold for a hit to be retained is a*max{T,c*log(l)}. [5.5]
        '''
    )
#################### 7、  处理bam文件 ####################
    parser.add_argument(
        "--bamdir", metavar="<path>", type=str,
        help='''存放bam文件的目录'''
    )

    parser.add_argument(
        "--delPCR", action="store_true"
    )

#################### 8、  变异检测 ####################
    parser.add_argument(
        "--processedbamdir", metavar="<path>", type=str,
        help='''存放处理过后bam文件的目录'''
    )
#################### 9、  jointgenotype ####################
    parser.add_argument(
        "--gvcfdir", metavar="<path>", type=str,
        help='''存放gvcf文件的目录'''
    )

#################### 10、  vcf质控 ####################
    parser.add_argument(
        "--vcffile", metavar="<filename>", type=str,
        help='''vcf.gz文件目录'''
    )
    parser.add_argument(
        "--snpQUAL", metavar="<str>", type=str,default="30.0",
        help='''
           snpQUAL
        '''
    )
    parser.add_argument(
        "--snpQD", metavar="<str>", type=str,default="2.0",
        help='''
            QualByDepth（QD）QD是变异质量值（Quality）除以覆盖深度（Depth）得到的比值。这里的变异质量值就是VCF中QUAL的值——用来衡量变异的可靠程度
        '''
    )
    parser.add_argument(
        "--snpMQ", metavar="<str>", type=str,default="40.0",
        help='''
            FisherStrandRMSMappingQuality
        '''
    )
    parser.add_argument(
        "--snpFS", metavar="<str>", type=str,default="60.0",
        help='''
            StrandOddsRatio
        '''
    )
    parser.add_argument(
        "--snpSOR", metavar="<str>", type=str,default="3.0",
        help='''
            StrandOddsRatio
        '''
    )
    parser.add_argument(
        "--snpMQRankSum", metavar="<str>", type=str,default="-12.5",
        help='''
            MQRankSum
        '''
    )
    parser.add_argument(
        "--snpReadPosRankSum", metavar="<str>", type=str,default="-8.0",
        help='''
            ReadPosRankSum
        '''
    )

    parser.add_argument(
        "--indelQUAL", metavar="<str>", type=str,default="30.0",
        help='''
           indelQUAL
        '''
    )
    parser.add_argument(
        "--indelQD", metavar="<str>", type=str,default="2.0",
        help='''
            indelQD
        '''
    )
    parser.add_argument(
        "--indelFS", metavar="<str>", type=str,default="200.0",
        help='''
            indelFS
        '''
    )
    parser.add_argument(
        "--indelSOR", metavar="<str>", type=str,default="10.0",
        help='''
            indelSOR
        '''
    )
    parser.add_argument(
        "--indelMQRankSum", metavar="<str>", type=str,default="-12.5",
        help='''
            indelMQRankSum
        '''
    )
    parser.add_argument(
        "--indelReadPosRankSum", metavar="<str>", type=str,default="-8.0",
        help='''
            indelReadPosRankSum
        '''
    )
########################################################################################################################
########################################################################################################################
########################################################################################################################
    #################### 第2部分 GWAS ####################
# #################### 2.1、  vcf转ped、map\bim、bed、fam ####################
    parser.add_argument(
        "--genotypefile", metavar="<filename>", type=str,
        help='''基因型文件'''
    )
    parser.add_argument(
        "--phenotypefile", metavar="<filename>", type=str,
        help='''表型文件'''
    )
# #################### 2.2、  gwas质控 ####################
    parser.add_argument(
        "--bfiledir", metavar="<path>", type=str,
        help='''
        bed文件目录
        '''
    )
    parser.add_argument(
        "--snpmiss", metavar="<path>", type=str,
        help='''
        snpmiss
        '''
    )
    parser.add_argument(
        "--atgc", action="store_true",
        help='''
        只保留atgc
        '''
    )
    parser.add_argument(
        "--indmiss", metavar="<str>", type=str,
        help='''
            indmiss
        '''
    )
    parser.add_argument(
        "--checksex", action="store_true",
        help="是否检查性别"
    )
    parser.add_argument(
        "--rmproblemsex", action="store_true",
        help="删除性别有问题的个体"
    )
    parser.add_argument(
        "--imputesex", action="store_true",
        help="根据基因型信息将性别进行填补"
    )
    parser.add_argument(
        "--maf", metavar="<str>", type=str,
        help='''
            maf
        '''
    )
    parser.add_argument(
        "--hwe", metavar="<str>", type=str,
        help='''
            对照组过滤hwe, --hwe 1e-6 
        '''
    )
    parser.add_argument(
        "--hweall", metavar="<str>", type=str,
        help='''
            过滤全部样本, --hweall 1e-6 
        '''
    )
    parser.add_argument(
        "--indep", metavar="<str>", type=str,nargs="+",
        help='''
            --indep <window size>['kb'] <step size (variant ct)> <VIF threshold>
            --indep  50 5 2
        '''
    )
    parser.add_argument(
        "--indepPairwise", metavar="<str>", type=str,nargs="+",
        help='''
            --indepPairwise  <window size>['kb'] <step size (variant ct)> <r^2 threshold>
            --indepPairwise  500 50 0.2
        '''
    )
    parser.add_argument(
        "--indepPairphase", metavar="<str>", type=str,nargs="+",
        help='''
            --indepPairphase  <window size>['kb'] <step size (variant ct)> <r^2 thresh>
        '''
    )

    parser.add_argument(
        "--heterozygosity", metavar="<str>", type=str,
        help='''
            我们建议删除样品杂合率平均值中偏离±3 SD的个体。
            --heterozygosity 3
        '''
    )
    # parser.add_argument(
    #     "--pihat", metavar="<str>", type=str,default="0.2",
    #     help='''删除PI_HAT值低于给定截止值的行'''
    # )
    # parser.add_argument(
    #     "--removerelationship", action="store_true",
    #     help='''删除亲缘关系个体
    #     比如我们分析的群体里面有亲子关系的个体，想要进行分析，不需要做这一步的筛选。'''
    # )


# #################### 2.3、  pca分析 ####################
    parser.add_argument(
        "--cleanbfiledir", metavar="<path>", type=str,
        help='''
        qc过后的bed、bim、fam文件目录
        '''
    )

    parser.add_argument(
        "--pcanum", metavar="<str>", type=str,default="6",
        help='''
            pcanum
        '''
    )
    parser.add_argument(
        "--groupnum", metavar="<str>", type=str,
        help='''
            groupnum，如果不设置分群数量，将从2~6之间，迭代cv error值最小的num作为分群数量
        '''
    )


# #################### 2.4、  kinship ####################

# #################### 2.5、  关联分析 ####################

    parser.add_argument(
        "--lm", action="store_true",
        help="广义线性模型"
    )
    parser.add_argument(
        "--lmm", action="store_true",
        help="混合线性模型"
    )
    parser.add_argument(
        "--pcafile", metavar="<filename>", type=str,
        help='''pcafile'''
    )
    parser.add_argument(
        "--kinshipfile", metavar="<filename>", type=str,
        help='''kinshipfile'''
    )
    parser.add_argument(
        "--pvaluelimit", metavar="<filename>", type=str,default="1e-7",
        help='''pvaluelimit'''
    )
# #################### 2.6、  selectsnp ####################
    parser.add_argument(
        "--assocfile",  metavar="<filename>", type=str,
        help="result.assoc.txt"
    )
########################################################################################################################
########################################################################################################################
########################################################################################################################
#################### 第3部分 预测 ####################
    parser.add_argument(
        "--snpfile", metavar="<filename>", type=str,
        help='''
        保存变异的染色体号，位置，ref，alt
        '''
    )

    parser.add_argument(
        "--species", metavar="<str>", type=str,
        help='''是动物或植物，qc中的check-sex'''
    )


    parser.add_argument(
        "--assembly", metavar="<str>", type=str,
        help='''数据库版本'''
    )



    return parser.parse_args(argv)
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
class Parser:
    def __init__(self, args):
        self.args = getArgs(args)
    def set_out_dir(self):
        return self.args.output
    def get_nosave(self):
        return self.args.nosave
    def get_nThrds(self):
        return self.args.nThrds
    def get_nMem(self):
        return self.args.nMem
        
    #################### 选择流程 ####################
    def get_step(self):
        return self.args.step


    #################### 第一部分 WGS ####################
#################### 1、  SRA数据下载 ####################

    def get_sra(self):
        return self.args.sra
    def get_sralist(self):
        return self.args.sralist
    def get_sradir(self):
        return self.args.sradir
#################### 2、  sra转fastq ####################
    def get_rawfastqdir(self):
        return self.args.rawfastqdir
#################### 3、  fastq质控 ####################
    def get_quality(self):
        return self.args.quality
    def get_phred(self):
        return self.args.phred
    def get_length(self):
        return self.args.length
    def get_stringency(self):
        return self.args.stringency
    def get_error(self):
        return self.args.error
#################### 4、  multiqc质量评估 ####################
    def get_fastqdir(self):
        return self.args.fastqdir

#################### 5、  建立基因组索引 ####################
    def get_accession(self):
        return self.args.accession
    def get_taxon(self):
        return self.args.taxon
        
    def get_refgenome(self):
        return self.args.refgenome
    def get_indexalgorithm(self):
        return self.args.indexalgorithm
    def get_refname(self):
        return self.args.refname
#################### 6、  比对参考基因组-bwa ####################
    def get_alignalgorithm(self):
        return self.args.alignalgorithm

    def get_cleanfastqdir(self):
        return self.args.cleanfastqdir
    # def get_rgpl(self):
    #     return self.args.rgpl



        
################### mem
    def get_minSeedLen(self):
        return self.args.minSeedLen
    def get_bandWidth(self):
        return self.args.bandWidth
    def get_zDropoff(self):
        return self.args.zDropoff
    def get_seedSplitRatio(self):
        return self.args.seedSplitRatio
    def get_maxOcc(self):
        return self.args.maxOcc
    def get_matchScore(self):
        return self.args.matchScore
    def get_mmPenalty(self):
        return self.args.mmPenalty
    def get_gapOpenPen(self):
        return self.args.gapOpenPen
    def get_gapExtPen(self):
        return self.args.gapExtPen
    def get_clipPen(self):
        return self.args.clipPen
    def get_unpairPen(self):
        return self.args.unpairPen
################### aln
    def get_maxDiff(self):
        return self.args.maxDiff
    def get_maxGapO(self):
        return self.args.maxGapO
    def get_maxGapE(self):
        return self.args.maxGapE
    def get_nDelTail(self):
        return self.args.nDelTail
    def get_nIndelEnd(self):
        return self.args.nIndelEnd
    def get_maxSeedDiff(self):
        return self.args.maxSeedDiff
    def get_seedLen(self):
        return self.args.seedLen
    def get_misMsc(self):
        return self.args.misMsc
    def get_gapOsc(self):
        return self.args.gapOsc
    def get_gapEsc(self):
        return self.args.gapEsc
    def get_trimQual(self):
        return self.args.trimQual
################### samse
    def get_maxHitPaired(self):
        return self.args.maxHitPaired
################### sampe
    def get_maxInsSize(self):
        return self.args.maxInsSize
    def get_maxHitDis(self):
        return self.args.maxHitDis
################### bwasw
    def get_mmPen(self):
        return self.args.mmPen
    def get_thres(self):
        return self.args.thres
    def get_hspIntv(self):
        return self.args.hspIntv
    def get_zBest(self):
        return self.args.zBest
    def get_nHspRev(self):
        return self.args.nHspRev
    def get_thresCoef(self):
        return self.args.thresCoef

#################### 7、  处理bam文件 ####################
    def get_bamdir(self):
        return self.args.bamdir
    def get_delPCR(self):
        return self.args.delPCR        

#################### 8、  变异检测 ####################
    def get_processedbamdir(self):
        return self.args.processedbamdir
        
#################### 9、  jointgenotype ####################
    def get_gvcfdir(self):
        return self.args.gvcfdir
#################### 10、  vcf质控 ####################

    def get_vcffile(self):
        return self.args.vcffile
    def get_snpQUAL(self):
        return self.args.snpQUAL
    def get_snpQD(self):
        return self.args.snpQD
    def get_snpMQ(self):
        return self.args.snpMQ
    def get_snpFS(self):
        return self.args.snpFS
    def get_snpSOR(self):
        return self.args.snpSOR
    def get_snpMQRankSum(self):
        return self.args.snpMQRankSum
    def get_snpReadPosRankSum(self):
        return self.args.snpReadPosRankSum
    def get_indelQUAL(self):
        return self.args.indelQUAL
    def get_indelQD(self):
        return self.args.indelQD
    def get_indelFS(self):
        return self.args.indelFS
    def get_indelSOR(self):
        return self.args.indelSOR
    def get_indelMQRankSum(self):
        return self.args.indelMQRankSum
    def get_indelReadPosRankSum(self):
        return self.args.indelReadPosRankSum
########################################################################################################################
########################################################################################################################
########################################################################################################################
#################### 第2部分 GWAS ####################
# #################### 2.1、  vcf转ped、map\bim、bed、fam ####################
    def get_genotypefile(self):
        return self.args.genotypefile
    def get_phenotypefile(self):
        return self.args.phenotypefile
# #################### 2.2、  基因型填充 ####################

# #################### 2.2、  gwas质控 ####################
    def get_bfiledir(self):
        return self.args.bfiledir
    def get_atgc(self):
        return self.args.atgc
    def get_snpmiss(self):
        return self.args.snpmiss
    def get_indmiss(self):
        return self.args.indmiss
    def get_checksex(self):
        return self.args.checksex
    def get_rmproblemsex(self):
        return self.args.rmproblemsex
    def get_imputesex(self):
        return self.args.imputesex
    def get_maf(self):
        return self.args.maf
    def get_hwe(self):
        return self.args.hwe
    def get_hweall(self):
        return self.args.hweall

    def get_indep(self):
        return self.args.indep
    def get_indepPairwise(self):
        return self.args.indepPairwise
    def get_indepPairphase(self):
        return self.args.indepPairphase
        
    def get_heterozygosity(self):
        return self.args.heterozygosity
    # def get_removerelationship(self):
    #     return self.args.removerelationship
    # def get_pihat(self):
    #     return self.args.pihat

# #################### 2.3、  pca分析 ####################
    def get_cleanbfiledir(self):
        return self.args.cleanbfiledir
    def get_pcanum(self):
        return self.args.pcanum
    def get_groupnum(self):
        return self.args.groupnum



# #################### 2.4、  kinship ####################

# #################### 2.5、  关联分析 ####################

    def get_lm(self):
        return self.args.lm
    def get_lmm(self):
        return self.args.lmm
    def get_pcafile(self):
        return self.args.pcafile
    def get_kinshipfile(self):
        return self.args.kinshipfile
    def get_pvaluelimit(self):
        return self.args.pvaluelimit
    def get_assocfile(self):
        return self.args.assocfile

########################################################################################################################
########################################################################################################################
########################################################################################################################
#################### 第3部分 预测 ####################
    def get_species(self):
        return self.args.species
    def get_assembly(self):
        return self.args.assembly
    def get_snpfile(self):
        return self.args.snpfile
    