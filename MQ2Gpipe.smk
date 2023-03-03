#!/usr/bin/env python
from concurrent.futures import thread
from runpy import run_module
import yaml
import re

##========= Globals ==========
configfile: 'config.yaml'

## Set samples
SAMPLEFILES = yaml.load(open(config['sample_list'], 'r'), Loader=yaml.FullLoader)
SAMPLES     = sorted(SAMPLEFILES.keys())    # get sample names
OUTPUTDIR   = config["output_dir"]          # get output dir

## Set reference file
DNA = config["dna"]
GTF = config["gtf"]
VCF = config["vcf"]
DICT = config["genome_dict"]
BWA_INDEX_PREFIX = config["genome_bwa_index"]

VEP_DB = config["vep_db"]
VEP_VER = config["vep_ver"]

##======== Rules ============
rule all:
    input:
    ## ========================================================================================= ##
    ## -------------------------- Module One: Data Clean Up      ------------------------------- ##
    ## ========================================================================================= ##
     # Step 00: Prepare data
       expand( OUTPUTDIR + "./DataCleanUp/00.DataPrepare/{sample}.R1.fq.gz", sample=SAMPLES ),
       expand( OUTPUTDIR + "./DataCleanUp/00.DataPrepare/{sample}.R2.fq.gz", sample=SAMPLES ),
     # Step 01: Quality Control
       expand( OUTPUTDIR + "./DataCleanUp/01.FastqFilter/{sample}/{sample}.R1.fq.gz", sample=SAMPLES ),
       expand( OUTPUTDIR + "./DataCleanUp/01.FastqFilter/{sample}/{sample}.R2.fq.gz", sample=SAMPLES ),
     # Step 02: Map to Genome
       expand( OUTPUTDIR + "./DataCleanUp/02.Map2Genome/{sample}.bam", sample=SAMPLES ),
     # Step 03: AddOrReplaceReadGroups
       expand( OUTPUTDIR + "./DataCleanUp/03.AddOrReplaceReadGroups/{sample}.add_rg.bam", sample=SAMPLES ),
     # Step 04: MarkDuplicates
       expand( OUTPUTDIR + "./DataCleanUp/04.MarkDuplicates/{sample}.dedupped.bam", sample=SAMPLES ),
       expand( OUTPUTDIR + "./DataCleanUp/04.MarkDuplicates/{sample}.dedupped.metrics", sample=SAMPLES ),
     # Step 05: BaseRecalibrator
       expand( OUTPUTDIR + "./DataCleanUp/05.BaseRecalibrator/{sample}.recal.tab", sample=SAMPLES ),
     # Step 06: ApplyBQSR
       expand( OUTPUTDIR + "./DataCleanUp/06.ApplyBQSR/{sample}.BQSR.bam", sample=SAMPLES ),
       expand( OUTPUTDIR + "./DataCleanUp/06.ApplyBQSR/{sample}.BQSR.bai", sample=SAMPLES ),
    # ========================================================================================= ##
    # -------------------------- Module Two: Variants Discover -------------------------------- ##
    # ========================================================================================= ##
      # Step 01: HaplotypeCaller
        expand( OUTPUTDIR + "./VarDiscover/01.HaplotypeCaller/{sample}.g.vcf.gz", sample=SAMPLES ),
      # Step 02: Consolidate GVCFs -- CombineGVCFs
        OUTPUTDIR + "./VarDiscover/02.CombineGVCFs/gvcf.list",
        OUTPUTDIR + "./VarDiscover/02.CombineGVCFs/cohort.script",
        OUTPUTDIR + "./VarDiscover/02.CombineGVCFs/cohort.g.vcf.gz",
      # Step 03: join-call -- GenotypeGVCFs
        OUTPUTDIR + "./VarDiscover/03.GenotypeGVCFs/RawGeno.script",
        OUTPUTDIR + "./VarDiscover/03.GenotypeGVCFs/RawGeno.vcf.gz",
        OUTPUTDIR + "./VarDiscover/03.GenotypeGVCFs/RawGeno.ok",
      # Step 04: SelectVariants
        OUTPUTDIR + "./VarDiscover/04.SelectVariants/SNP/RawGeno.SNP.script",
        OUTPUTDIR + "./VarDiscover/04.SelectVariants/SNP/RawGeno.SNP.vcf.gz",
        OUTPUTDIR + "./VarDiscover/04.SelectVariants/SNP/RawGeno.SNP.ok",
        OUTPUTDIR + "./VarDiscover/04.SelectVariants/INDEL/RawGeno.Indel.script",
        OUTPUTDIR + "./VarDiscover/04.SelectVariants/INDEL/RawGeno.Indel.vcf.gz",
        OUTPUTDIR + "./VarDiscover/04.SelectVariants/INDEL/RawGeno.Indel.ok",
      # Step 05: Variant Filtration
        OUTPUTDIR + "./VarDiscover/05.VariantFiltration/SNP/snps_filtered.script",
        OUTPUTDIR + "./VarDiscover/05.VariantFiltration/SNP/snps_filtered.vcf.gz",
        OUTPUTDIR + "./VarDiscover/05.VariantFiltration/SNP/snps_filtered.ok",
        OUTPUTDIR + "./VarDiscover/05.VariantFiltration/INDEL/indels_filtered.script",
        OUTPUTDIR + "./VarDiscover/05.VariantFiltration/INDEL/indels_filtered.vcf.gz",
        OUTPUTDIR + "./VarDiscover/05.VariantFiltration/INDEL/indels_filtered.ok",
      # Step 06: Merge SNP and Indel
        OUTPUTDIR + "./VarDiscover/06.MergeVcfs/SNP_Indel.Filtered.src",
        OUTPUTDIR + "./VarDiscover/06.MergeVcfs/SNP_Indel.Filtered.vcf.gz",
        OUTPUTDIR + "./VarDiscover/06.MergeVcfs/SNP_Indel.Filtered.ok",
      # Step 07: Fetch Passed vcf 
        OUTPUTDIR + "./VarDiscover/07.PurePassedVCF/SNP_Indel.Passed.vcf.gz",
        OUTPUTDIR + "./VarDiscover/07.PurePassedVCF/SNP_Indel.Passed.vcf.gz.tbi",
      # Step 08: Rename Passed vcf
        OUTPUTDIR + "./VarDiscover/08.VarAddReName/SNP_Indel_AddName.vcf.gz",
      # Step 09: Calculate het rate
        OUTPUTDIR + "./VarDiscover/09.VariantStat/calc_het_rate.script",
        OUTPUTDIR + "./VarDiscover/09.VariantStat/variants.hwe.gz",
        OUTPUTDIR + "./VarDiscover/09.VariantStat/high_het_rate_var.list",
        OUTPUTDIR + "./VarDiscover/09.VariantStat/high_het_rate_var.ok",
      # Step 10: Filter And Fetch Final Genotype
        OUTPUTDIR + "./VarDiscover/10.FinalGenotype/Genotype.vcf.gz",
    ## ========================================================================================= ##
    ## -------------------------- Module Three: Variants Annotation ---------------------------- ##
    ## ========================================================================================= ##
      # Step 01: Annotation
        OUTPUTDIR + "./Annotaion/01.VepAnnotation/Annotation.script",
        OUTPUTDIR + "./Annotaion/01.VepAnnotation/VarWithAnno.vcf.gz",
        OUTPUTDIR + "./Annotaion/01.VepAnnotation/Annotation_Summary.html",
        OUTPUTDIR + "./Annotaion/01.VepAnnotation/Annotation_Waring.txt",
        OUTPUTDIR + "./Annotaion/01.VepAnnotation/Variant_Annotation.ok",
    ## ========================================================================================= ##
    ## -------------------------- Module Four: Variants Imputation ----------------------------- ##
    ## ========================================================================================= ##
      # Step 01: Beagle Imputation
        OUTPUTDIR + "./Imputation/01.BeagleImpute/BeagleImputed.vcf.gz",
# ............................................................................................... #
# =========================================================================== #
# ``````````````````````````````````````````````````````````````````````````` #
#                      Part 01 Data Preparation                               #
# ........................................................................... #
# =========================================================================== #
# Step 00: Prepare data
rule DataCleanUp_00_PrepareFastqFile:
    input:
        lambda wildcards:SAMPLEFILES[wildcards.sample]
    output:
       R1 =  OUTPUTDIR + "./DataCleanUp/00.DataPrepare/{sample}.R1.fq.gz",
       R2 =  OUTPUTDIR + "./DataCleanUp/00.DataPrepare/{sample}.R2.fq.gz"
    threads:
        1
    shell:
        """
        ln -sf {input[0]} {output.R1} && ln -sf {input[1]} {output.R2}
        """

# Step 01: Quality Control
rule DataCleanUp_01_FastqFilter:
    input:
        R1 =  OUTPUTDIR + "./DataCleanUp/00.DataPrepare/{sample}.R1.fq.gz",
        R2 =  OUTPUTDIR + "./DataCleanUp/00.DataPrepare/{sample}.R2.fq.gz"
    output:
        R1   = OUTPUTDIR + "./DataCleanUp/01.FastqFilter/{sample}/{sample}.R1.fq.gz",
        R2   = OUTPUTDIR + "./DataCleanUp/01.FastqFilter/{sample}/{sample}.R2.fq.gz",
        json = OUTPUTDIR + "./DataCleanUp/01.FastqFilter/{sample}/{sample}.json",
        html = OUTPUTDIR + "./DataCleanUp/01.FastqFilter/{sample}/{sample}.html"
    message:
        "Begin to filter fastq!"
    log:
        OUTPUTDIR + "./AllLogs/DataCleanUp/01.FastqFilter/{sample}.FastqFilter.log"
    benchmark:
        OUTPUTDIR + "./Benchmark/DataCleanUp/01.FastqFilter/{sample}.benchmark"
    params:
        "--detect_adapter_for_pe --fix_mgi_id"
    threads:
        8
    shell:
        """
        bin/fastp {params} -w {threads} -i {input.R1} -I {input.R2} -o {output.R1} \
        -O {output.R2} -j {output.json} -h {output.html} 2> {log}
        """

# Step 02: map to genome by BWA
rule DataCleanUp_02_Map2Genome:
    input:
        R1   = OUTPUTDIR + "./DataCleanUp/01.FastqFilter/{sample}/{sample}.R1.fq.gz",
        R2   = OUTPUTDIR + "./DataCleanUp/01.FastqFilter/{sample}/{sample}.R2.fq.gz",
        IDX  = BWA_INDEX_PREFIX,
    output:
        bam = OUTPUTDIR + "./DataCleanUp/02.Map2Genome/{sample}.bam"
    log:
        OUTPUTDIR + "./AllLogs/DataCleanUp/02.Map2Genome/{sample}.bwa_align.logs"
    benchmark:
        OUTPUTDIR + "./Benchmark/DataCleanUp/02.Map2Genome/{sample}.benchmark"
    threads:
        8 
    params:
        "-k 31 -P -M"
    shell:
        """
        bwa mem {params} -t 8 {input.IDX} {input.R1} {input.R2} 2>{log} \
             | samtools view -Sb -@ 4 -m 5G -o {output.bam}
        """

# Step 03: AddOrReplaceReadGroups
rule DataCleanUp_03_AddOrReplaceReadGroups:
    input:
        bam = OUTPUTDIR + "./DataCleanUp/02.Map2Genome/{sample}.bam",
    output:
        bam = OUTPUTDIR + "./DataCleanUp/03.AddOrReplaceReadGroups/{sample}.add_rg.bam",
    log:
        OUTPUTDIR + "./AllLogs/DataCleanUp/03.AddOrReplaceReadGroups/{sample}.add_rg.log"
    benchmark:
        OUTPUTDIR + "./Benchmark/DataCleanUp/03.AddOrReplaceReadGroups/{sample}.benchmark"
    threads:
        12
    params:
        "-SO coordinate -RGID {sample} -RGLB lib1 -RGPL MGISEQ -RGPU hiseq -RGSM {sample}"
    shell:
        """
        bin/gatk-4.2.6.1/gatk --java-options '-Xms5G -Xmx15G' AddOrReplaceReadGroups \
            {params} -I {input.bam} -O {output.bam} 2> {log}
        """

# Step 04: MarkDuplicates
rule DataCleanUp_04_MarkDuplicates:
    input:
        bam = OUTPUTDIR + "./DataCleanUp/03.AddOrReplaceReadGroups/{sample}.add_rg.bam",
    output:
        bam = OUTPUTDIR + "./DataCleanUp/04.MarkDuplicates/{sample}.dedupped.bam",
        mtx = OUTPUTDIR + "./DataCleanUp/04.MarkDuplicates/{sample}.dedupped.metrics",
    message:
        "Begin to MarkDuplicates by Picard!"
    log:
        OUTPUTDIR + "./AllLogs/DataCleanUp/04.MarkDuplicates/{sample}.dedup.log"
    benchmark:
        OUTPUTDIR + "./Benchmark/DataCleanUp/04.MarkDuplicates/{sample}.benchmark"
    threads:
        8
    params:
        "--CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"
    shell:
        """
        bin/gatk-4.2.6.1/gatk --java-options '-Xms5G -Xmx15G' MarkDuplicates \
                {params} --INPUT {input.bam} --OUTPUT {output.bam} --METRICS_FILE {output.mtx}
        """

# Step 05: BaseRecalibrator
rule DataCleanUp_05_BaseRecalibrator:
    input:
        dna = DNA,
        vcf = VCF,
        bam = OUTPUTDIR + "./DataCleanUp/04.MarkDuplicates/{sample}.dedupped.bam",
    output:
        tab = OUTPUTDIR + "./DataCleanUp/05.BaseRecalibrator/{sample}.recal.tab",
    log:
        OUTPUTDIR + "./AllLogs/DataCleanUp/05.BaseRecalibrator/{sample}.log"
    benchmark:
        OUTPUTDIR + "./Benchmark/DataCleanUp/05.BaseRecalibrator/{sample}.benchmark"
    message:
        "Begin to BaseRecalibrator by GATK !"
    threads:
        8
    params:
        "--use-original-qualities"
    shell:
        """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms5G -Xmx10G" BaseRecalibrator {params} \
            -R {input.dna} -I {input.bam} -O {output.tab} -known-sites {input.vcf}
        """

# Step 06: ApplyBQSR
rule DataCleanUp_06_ApplyBQSR:
    input:
        dna = DNA,
        bam = OUTPUTDIR + "./DataCleanUp/04.MarkDuplicates/{sample}.dedupped.bam",
        tab = OUTPUTDIR + "./DataCleanUp/05.BaseRecalibrator/{sample}.recal.tab",
    output:
        bam = OUTPUTDIR + "./DataCleanUp/06.ApplyBQSR/{sample}.BQSR.bam",
        bai = OUTPUTDIR + "./DataCleanUp/06.ApplyBQSR/{sample}.BQSR.bai",
    log:
        OUTPUTDIR + "./AllLogs/DataCleanUp/06.ApplyBQSR/{sample}.BQSR.log"
    benchmark:
        OUTPUTDIR + "./Benchmark/DataCleanUp/06.ApplyBQSR/{sample}.benchmark"
    resources:
        mem_mb=10000
    threads:
        8
    params:
        "--add-output-sam-program-record --use-original-qualities"
    shell:
        """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms5G -Xmx10G" ApplyBQSR {params} \
            -R {input.dna} -I {input.bam} -O {output.bam} --bqsr-recal-file {input.tab}
        """

# =========================================================================== #
# ``````````````````````````````````````````````````````````````````````````` #
#                      Part 02 Variants Discovery                             #
# ........................................................................... #
# =========================================================================== #

# Step 01: HaplotypeCaller
rule VarDiscover_01_HaplotypeCaller:
    input:
        dna = DNA,
        bam = OUTPUTDIR + "./DataCleanUp/06.ApplyBQSR/{sample}.BQSR.bam",
    output:
        vcf = OUTPUTDIR + "./VarDiscover/01.HaplotypeCaller/{sample}.g.vcf.gz"
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/01.HaplotypeCaller/{sample}.g.vcf.log"
    benchmark:
        OUTPUTDIR + "./Benchmark/VarDiscover/01.HaplotypeCaller/{sample}.benchmark"
    resources:
        mem_mb=10000
    threads:
        4
    params:
        "--emit-ref-confidence GVCF"
    shell:
        """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads={threads}" \
            HaplotypeCaller {params} -R {input.dna} -I {input.bam} -O {output.vcf}
        """

# Step 02: Consolidate GVCFs
rule VarDiscover_02_CombineGVCFs:
    input:
        dna = DNA,
        gvcf = expand(OUTPUTDIR + "./VarDiscover/01.HaplotypeCaller/{sample}.g.vcf.gz", sample = SAMPLES),
    output:
        lst = OUTPUTDIR + "./VarDiscover/02.CombineGVCFs/gvcf.list",
        src = OUTPUTDIR + "./VarDiscover/02.CombineGVCFs/cohort.script",
        vcf = OUTPUTDIR + "./VarDiscover/02.CombineGVCFs/cohort.g.vcf.gz",
        ok  = OUTPUTDIR + "./VarDiscover/02.CombineGVCFs/cohort.ok",
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/02.CombineGVCFs/cohort.g.vcf.log"
    benchmark:
        OUTPUTDIR + "./Benchmark/VarDiscover/02.CombineGVCFs/cohort.benchmark"
    resources:
        mem_mb=200000
    threads:
        16
    params:
        dir = OUTPUTDIR + "./VarDiscover/02.CombineGVCFs/"
    run:
        import os
        import subprocess

        if not os.path.exists(params.dir):
            os.makedirs(params.dir)
            
        with open(output.lst, 'w') as f:
            for file in input.gvcf:
                print(file, file=f)
        
        cmd = """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms50G -Xmx200G -XX:ParallelGCThreads={t}" \
            CombineGVCFs -R {r} -V {l} -O {o}
        """.format(t=threads, r=input.dna, l=output.lst, o=output.vcf)

        with open(output.src, "w") as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.vcf) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)
        
# Step 03: Joint-Call Cohort
rule VarDiscover_03_GenotypeGVCFs:
    input:
        dna = DNA,
        vcf = OUTPUTDIR + "./VarDiscover/02.CombineGVCFs/cohort.g.vcf.gz",
    output:
        src = OUTPUTDIR + "./VarDiscover/03.GenotypeGVCFs/RawGeno.script",
        vcf = OUTPUTDIR + "./VarDiscover/03.GenotypeGVCFs/RawGeno.vcf.gz",
        ok  = OUTPUTDIR + "./VarDiscover/03.GenotypeGVCFs/RawGeno.ok",
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/03.GenotypeGVCFs/GenotypeGVCFs.log"
    benchmark:
        OUTPUTDIR + "./Benchmark/VarDiscover/03.GenotypeGVCFs/GenotypeGVCFs.benchmark"
    resources:
        mem_mb=150000
    threads:
        16
    params:
        dir = OUTPUTDIR + "./VarDiscover/03.GenotypeGVCFs/"
    run:
        import os
        import subprocess
        
        if not os.path.exists(params.dir):
            os.makedirs(params.dir)

        cmd = """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms10G -Xmx150G -XX:ParallelGCThreads={t}" \
            GenotypeGVCFs -R {r} -V {v} -O {o} 2> {l}
        """.format(t=threads, r=input.dna, v=input.vcf, o=output.vcf, l=log)

        with open(output.src, "w") as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.vcf) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)

# Step 04: GATK SelectVariants -- SNP
rule VarDiscover_04_SelectVariantsSNP:
    input:
        vcf = OUTPUTDIR + "./VarDiscover/03.GenotypeGVCFs/RawGeno.vcf.gz",
    output:
        src = OUTPUTDIR + "./VarDiscover/04.SelectVariants/SNP/RawGeno.SNP.script",
        vcf = OUTPUTDIR + "./VarDiscover/04.SelectVariants/SNP/RawGeno.SNP.vcf.gz",
        ok  = OUTPUTDIR + "./VarDiscover/04.SelectVariants/SNP/RawGeno.SNP.ok",
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/04.SelectVariants/SNP/RawGeno.SNP.log"
    resources:
        mem_mb=20000
    threads:
        4
    params:
        opt = "--select-type-to-include SNP"
    run:
        import os
        import subprocess
        
        outputpath = os.path.dirname(output.vcf)
        if not os.path.exists(outputpath):
            os.makedirs(params.dir)

        cmd = """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms5G -Xmx20G -XX:ParallelGCThreads={t}" \
            SelectVariants {p} -V {v} -O {o} 2> {l}
        """.format(t=threads, p=params.opt, v=input.vcf, o=output.vcf, l=log)

        with open(output.src, "w") as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.vcf) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)

# Step 04: GATK SelectVariants -- INDEL
rule VarDiscover_04_SelectVariantsIndel:
    input:
        vcf = OUTPUTDIR + "./VarDiscover/03.GenotypeGVCFs/RawGeno.vcf.gz",
    output:
        src = OUTPUTDIR + "./VarDiscover/04.SelectVariants/INDEL/RawGeno.Indel.script",
        vcf = OUTPUTDIR + "./VarDiscover/04.SelectVariants/INDEL/RawGeno.Indel.vcf.gz",
        ok  = OUTPUTDIR + "./VarDiscover/04.SelectVariants/INDEL/RawGeno.Indel.ok",
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/04.SelectVariants/INDEL/RawGeno.Indel.log"
    resources:
        mem_mb=20000
    threads:
        4
    params:
        opt = "--select-type-to-include INDEL --max-indel-size 40"
    run:
        import os
        import subprocess
        
        outputpath = os.path.dirname(output.vcf)
        if not os.path.exists(outputpath):
            os.makedirs(outputpath)

        cmd = """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms5G -Xmx20G -XX:ParallelGCThreads={t}" \
            SelectVariants {p} -V {v} -O {o} 2> {l}
        """.format(t=threads, p=params.opt, v=input.vcf, o=output.vcf, l=log)

        with open(output.src, "w") as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.vcf) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)

# Step 05:  Hard-filter SNPs on multiple expressions using VariantFiltration
rule VarFilter_05_HardFilterSNPs:
    input:
        vcf = OUTPUTDIR + "./VarDiscover/04.SelectVariants/SNP/RawGeno.SNP.vcf.gz",
    output:
        src = OUTPUTDIR + "./VarDiscover/05.VariantFiltration/SNP/snps_filtered.script",
        vcf = OUTPUTDIR + "./VarDiscover/05.VariantFiltration/SNP/snps_filtered.vcf.gz",
        ok  = OUTPUTDIR + "./VarDiscover/05.VariantFiltration/SNP/snps_filtered.ok",
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/05.VariantFiltration/SNP/snp_filter.log"
    resources:
        mem_mb=20000
    threads:
        8
    params:
        opt = """ -cluster 3 -window 10 \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"   """
    run:
        import os
        import subprocess
        
        outputpath = os.path.dirname(output.vcf)
        if not os.path.exists(outputpath):
            os.makedirs(outputpath)

        cmd = """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms10G -Xmx20G -XX:ParallelGCThreads={t}" \
            VariantFiltration {p} -V {v} -O {o} 2> {l}
        """.format(t=threads, p=params.opt, v=input.vcf, o=output.vcf, l=log)

        with open(output.src, "w") as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.vcf) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)

# Step 05:  Hard-filter Indels -- VariantFiltration
rule VarFilter_05_HardFilterIndels:
    input:
        vcf = OUTPUTDIR + "./VarDiscover/04.SelectVariants/INDEL/RawGeno.Indel.vcf.gz",
    output:
        src = OUTPUTDIR + "./VarDiscover/05.VariantFiltration/INDEL/indels_filtered.script",
        vcf = OUTPUTDIR + "./VarDiscover/05.VariantFiltration/INDEL/indels_filtered.vcf.gz",
        ok  = OUTPUTDIR + "./VarDiscover/05.VariantFiltration/INDEL/indels_filtered.ok",
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/05.VariantFiltration/INDEL/indels_filter.log"
    resources:
        mem_mb=10000
    threads:
        8
    params:
        opt = """ -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 10.0" --filter-name "SOR10" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "InbreedingCoeff < -0.8" --filter-name "InbreedingCoeff-0.8" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" """
    run:
        import os
        import subprocess
        
        outputpath = os.path.dirname(output.vcf)
        if not os.path.exists(outputpath):
            os.makedirs(outputpath)

        cmd = """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms5G -Xmx10G -XX:ParallelGCThreads={t}" \
            VariantFiltration {p} -V {v} -O {o} 2> {l}
        """.format(t=threads, p=params.opt, v=input.vcf, o=output.vcf, l=log)

        with open(output.src, "w") as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.vcf) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)

# Step 06: Merge SNP and Indel
rule VarFilter_06_MergeVcfs:
    input:
        snp   = OUTPUTDIR + "./VarDiscover/05.VariantFiltration/SNP/snps_filtered.vcf.gz",
        indel = OUTPUTDIR + "./VarDiscover/05.VariantFiltration/INDEL/indels_filtered.vcf.gz",
    output:
        src = OUTPUTDIR + "./VarDiscover/06.MergeVcfs/SNP_Indel.Filtered.src",
        vcf = OUTPUTDIR + "./VarDiscover/06.MergeVcfs/SNP_Indel.Filtered.vcf.gz",
        ok  = OUTPUTDIR + "./VarDiscover/06.MergeVcfs/SNP_Indel.Filtered.ok"
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/01.MergeVcfs/snp_indel_merge.log"
    resources:
        mem_mb=10000
    threads:
        4
    run:
        import os
        import subprocess
        
        outputpath = os.path.dirname(output.vcf)
        if not os.path.exists(outputpath):
            os.makedirs(outputpath)
        
        cmd = """
        bin/gatk-4.2.6.1/gatk --java-options "-Xms5G -Xmx10G -XX:ParallelGCThreads={t}" \
            MergeVcfs -I {snp} -I {indel} -O {o} 2> {l}
        """.format(t=threads, snp=input.snp, indel=input.indel, o=output.vcf, l=log)

        with open(output.src, "w") as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.vcf) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)

# Step 07: Fetch Passed Variants
rule VarFilter_07_PurePassedVCF:
    input:
        vcf = OUTPUTDIR + "./VarDiscover/06.MergeVcfs/SNP_Indel.Filtered.vcf.gz",
    output:
        src = OUTPUTDIR + "./VarDiscover/07.PurePassedVCF/SNP_Indel.Passed.script",
        vcf = OUTPUTDIR + "./VarDiscover/07.PurePassedVCF/SNP_Indel.Passed.vcf.gz",
        tbi = OUTPUTDIR + "./VarDiscover/07.PurePassedVCF/SNP_Indel.Passed.vcf.gz.tbi",
        ok = OUTPUTDIR + "./VarDiscover/07.PurePassedVCF/SNP_Indel.Passed.vcf.gz.ok",
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/07.PurePassedVCF/SNP_Indel.Passed.log"
    threads:
        8
    params:
        opt_filter = " -g10 -G10 ",
        opt_view = " -f 'PASS,.' -O z -m2 -M2 "
    run:
        import os
        import subprocess

        outputpath = os.path.dirname(output.vcf)
        if not os.path.exists(outputpath):
            os.makedirs(outputpath)

        cmd = """
        bcftools filter --threads {t} {p1} {v1} | bcftools view --threads {t} {p2} - -o {v2} 2> {l} && \
        bcftools index -t --threads {t} -o {b} {v2}
        """.format(t=threads, p1=params.opt_filter, v1=input.vcf, p2=params.opt_view, 
            v2=output.vcf, b=output.tbi, l=log)

        with open(output.src, 'w') as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.tbi) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)

# Step 08: Variants Rename
rule VarPure_08_VarAddReName:
    input:
        vcf = OUTPUTDIR + "./VarDiscover/07.PurePassedVCF/SNP_Indel.Passed.vcf.gz",
    output:
        vcf = OUTPUTDIR + "./VarDiscover/08.VarAddReName/SNP_Indel_AddName.vcf.gz",
        tbi = OUTPUTDIR + "./VarDiscover/08.VarAddReName/SNP_Indel_AddName.vcf.gz.tbi",
    threads:
        8
    shell:
        """
        perl bin/VarAddName.pl {input.vcf} |bgzip > {output.vcf} && \
        bcftools index --tbi --threads {threads} -o {output.tbi} {output.vcf}
        """
# Step 09: Filter Twice And get final genotype
rule VarPure_09_CalculateHetRate:
    input:
        vcf = OUTPUTDIR + "./VarDiscover/08.VarAddReName/SNP_Indel_AddName.vcf.gz",
    output:
        src = OUTPUTDIR + "./VarDiscover/09.VariantStat/calc_het_rate.script",
        hwe = OUTPUTDIR + "./VarDiscover/09.VariantStat/variants.hwe.gz",
        het = OUTPUTDIR + "./VarDiscover/09.VariantStat/high_het_rate_var.list",
        ok = OUTPUTDIR + "./VarDiscover/09.VariantStat/high_het_rate_var.ok",
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/09.VariantStat/high_het_rate_var.log"
    params:
        opt = "--hardy gz --allow-extra-chr",
        pfx = OUTPUTDIR + "./VarDiscover/09.VariantStat/variants"
    threads:
        16
    run:
        import os
        import subprocess

        outputpath = os.path.dirname(output.het)

        if not os.path.exists(outputpath):
            os.makedirs(outputpath)

        cmd = """
        bin/plink {p} --threads {t} --vcf {v} --out {pfx} && \
        zcat {hwe} | perl bin/fetch_high_het_rate.pl - > {het}
        """.format(p=params.opt, t=threads, v=input.vcf, pfx=params.pfx, 
        hwe=output.hwe, het=output.het, l=log)

        with open(output.src, 'w') as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.het) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)

# Step 10: Fetch final genotype
rule VarPure_10_FinalGenotype:
    input:
        vcf = OUTPUTDIR + "./VarDiscover/08.VarAddReName/SNP_Indel_AddName.vcf.gz",
        het = OUTPUTDIR + "./VarDiscover/09.VariantStat/high_het_rate_var.list",
    output:
        src = OUTPUTDIR + "./VarDiscover/10.FinalGenotype/fetch_final.src",
        vcf = OUTPUTDIR + "./VarDiscover/10.FinalGenotype/Genotype.vcf.gz",
        ok  = OUTPUTDIR + "./VarDiscover/10.FinalGenotype/fetch_final.ok",
    log:
        OUTPUTDIR + "./AllLogs/VarDiscover/10.FinalGenotype/fetch_final.log"
    threads:
        24
    params:
        pfx = OUTPUTDIR + "./VarDiscover/10.FinalGenotype/Genotype",
        opt = "--geno 0.2 --maf 0.01 --recode vcf-iid bgz --allow-extra-chr"
    run:
        import os
        import subprocess

        outputpath = os.path.dirname(output.vcf)

        if not os.path.exists(outputpath):
            os.makedirs(outputpath)

        cmd = """
        bin/plink {opt} --threads {t} --vcf {v} --exclude {het} --out {p} 2> {l}
        """.format(opt=params.opt, t=threads, v=input.vcf, het=input.het, 
        p=params.pfx, l=log)

        with open(output.src, 'w') as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.vcf) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)


# =========================================================================== #
# ``````````````````````````````````````````````````````````````````````````` #
#                      Part 03 Variants Annotation                            #
# ........................................................................... #
# =========================================================================== #                    
# Step 01: Annotation Variations
rule Annotation_01_VepAnnotation:
    input:
        vep_db = VEP_DB,
        vcf = OUTPUTDIR + "./VarDiscover/08.VarAddReName/SNP_Indel_AddName.vcf.gz",
    output:
        src = OUTPUTDIR + "./Annotaion/01.VepAnnotation/Annotation.script",
        ann = OUTPUTDIR + "./Annotaion/01.VepAnnotation/VarWithAnno.vcf.gz",
        tbi = OUTPUTDIR + "./Annotaion/01.VepAnnotation/VarWithAnno.vcf.gz.tbi",
        sum = OUTPUTDIR + "./Annotaion/01.VepAnnotation/Annotation_Summary.html",
        war = OUTPUTDIR + "./Annotaion/01.VepAnnotation/Annotation_Waring.txt",
        ok  = OUTPUTDIR + "./Annotaion/01.VepAnnotation/Variant_Annotation.ok",
    log:
        OUTPUTDIR + "./AllLogs/Annotaion/01.VepAnnotation/Annotation.log"
    threads:
        8
    params:
        opt = "--species zea_mays --distance 2000,0 --offline --cache --vcf --force_overwrite",
    run:
        import os
        import subprocess

        outputpath = os.path.dirname(output.ann)
        if not os.path.exists(outputpath):
            os.makedirs(outputpath)

        # Annotation with VEP program
        vep_dir = os.path.dirname(input.vep_db)
        cmd = """
        source activate genome_tools && \
        vep {p} -i {i} --dir {d} --cache_version {v} --warning_file {w} \
            --stats_file {s} -o - 2> {l} | bgzip > {o} && \
        bcftools index --tbi --threads {t} -o {b} {o}
        """.format(p=params.opt, i=input.vcf, d=vep_dir, v=VEP_VER, o=output.ann,
        w=output.war, s=output.sum, l=log, t=threads, b=output.tbi)

        with open(output.src, 'w') as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.ann) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)

# =========================================================================== #
# ``````````````````````````````````````````````````````````````````````````` #
#                      Part 04 Variants Imputation                            #
# ........................................................................... #
# =========================================================================== #
# Step 01: Beagle Imputation
rule Imputation_02_BeagleImpute:
    input:
        vcf = OUTPUTDIR + "./VarDiscover/10.FinalGenotype/Genotype.vcf.gz",
    output:
        src = OUTPUTDIR + "./Imputation/01.BeagleImpute/BeagleImputed.src",
        vcf = OUTPUTDIR + "./Imputation/01.BeagleImpute/BeagleImputed.vcf.gz",
        ok =  OUTPUTDIR + "./Imputation/01.BeagleImpute/BeagleImputed.ok",
    log:
        OUTPUTDIR + "./AllLogs/Imputation/01.BeagleImpute/BeagleImpute.log"
    threads:
        24
    params:
        pfx = OUTPUTDIR + "./Imputation/01.BeagleImpute/BeagleImputed"
    run:
        import os
        import subprocess

        outputpath = os.path.dirname(output.vcf)
        if not os.path.exists(outputpath):
            os.makedirs(outputpath)
        
        cmd = """
        java -Xms5g -Xmx20g -jar bin/beagle.jar nthreads={t} gt={i} out={p} 2> {log}
        """.format(t=threads, i=input.vcf, o=params.pfx, l=log)

        with open(output.src, 'w') as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.vcf) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)