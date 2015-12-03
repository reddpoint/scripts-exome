### Takes folders containing FASTQ.gz and returns

"""
Needs file structure:
    ./GA800x/reads/ - folder containing all reads
    files should be structured GA800x_lane#_run#_R{1,2}.fastq.gz
"""


configfile: "config.yaml"
workdir: config["WORKINGDIR"]

SAMPLENAMES = set(glob_wildcards("{FOLDERS}/reads/{SAMPLES,[a-zA-Z0-9]+}_{LANE,[a-zA-Z0-9]+}_{RUN,[a-zA-Z0-9]+}_{R}.fastq.gz").SAMPLES)


rule all:
    input:
        expand("{S}/merged/{S}.bam", S=SAMPLENAMES)

rule bwa:
    input:
        forward = "{FOLDERS}/reads/{SAMPLE_RUN_LANE}_R1.fastq",
        reverse = "{FOLDERS}/reads/{SAMPLE_RUN_LANE}_R2.fastq"
    output:
        temp("{FOLDERS}/aligned/{SAMPLE_RUN_LANE}.sam")
    threads: config["threads"]
    params:
        ref = config["REFERENCE_SEQUENCE"]
    shell:
        """
        module load bwa/0.7.10
        bwa mem -Mp -t {params.ref} {input.forward} {input.reverse} > {output}
        """

rule AddOrReplaceReadGroups:
    input:
        "{FOLDERS}/aligned/{SAMPLE_RUN_LANE}.sam"
    output:
        temp("{FOLDERS}/aligned/{SAMPLE_RUN_LANE}.bam")
    shell:
        """
        module load picard
        java -Xmx8g -Djava.io.tmpdir=/scratch/cooperjam/temp -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups \
        INPUT={input} \
        OUTPUT={output} \
        ID= {input.SAMPLE_RUN_LANE} \
        SM= {input.FOLDERS} \
        LB=1 \
        PL=il \
        PU=1 \
        SORT_ORDER=coordinate
        """

rule MarkDuplicates:
    input:
        "{FOLDERS}/aligned/{SAMPLE_RUN_LANE}.bam"
    output:
        "{FOLDERS}/aligned_dup_mark_bam/{SAMPLE_RUN_LANE}.dup.mark.bam"
    log:
        "{FOLDERS}/aligned_dup_mark_bam/{SAMPLE_RUN_LANE}.dup.mark.log"
    shell:
        """
        module load picard
        module load samtools
        java -Xmx8g -Djava.io.tmpdir=/scratch/cooperjam/temp -jar $PICARDJARPATH/picard.jar MarkDuplicates \
        INPUT={input} \
        OUTPUT={output} \
        METRICS_FILE="{FOLDERS}/aligned_dup_mark_bam/{SAMPLE_RUN_LANE}.dup.mark.metrics.txt" \
        """

rule bam_index:
    input:
        "{FOLDERS}/aligned_dup_mark_bam/{SAMPLE_RUN_LANE}.dup.mark.bam"
    output:
        "{FOLDERS}/aligned_dup_mark_bam/{SAMPLE_RUN_LANE}.dup.mark.bam.bai"
    shell:
        """
        module load samtools
        samtools index {input}
        """

def _get_ref(wildcards):
    return config["REFERENCE_SEQUENCE"][wildcards.reference]


rule gatk_realign_info:
    input:
        "{FOLDERS}/aligned_dup_mark_bam/{SAMPLE_RUN_LANE}.dup.mark.bam.bai",
        ref=_get_ref,
        bam="{FOLDERS}/aligned_dup_mark_bam/{SAMPLE_RUN_LANE}.dup.mark.bam"
    output:
        temp("{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.realign.intervals")
    log:
        "{FOLDERS}/mapping/log/{SAMPLE_RUN_LANE}.realign_info.log"
    threads: config['THREADS']
    shell:
        """
        java -Xmx8g -Djava.io.tmpdir=/scratch/cooperjam/temp -jar $GATK_HOME/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R {input.ref} \
        -nt {threads} \
        -I {input.bam} \
        -L {config[INTERVAL]} \
        -known {config[known_variants][dbsnp]} \
        -o {output} >& {log}
        """

rule gatk_realign_bam:
    input:
        ref=_get_ref,
        bam="{FOLDERS}/aligned_dup_mark_bam/{SAMPLE_RUN_LANE}.dup.mark.bam",
        intervals="{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.realign.intervals""
    output:
        "{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.realigned.bam"
    params:
        custom=config.get("params_gatk", "")
    log:
        "{FOLDERS}/mapping/log/{SAMPLE_RUN_LANE}.realign.log"
    shell:
        """
        java -Xmx8g -Djava.io.tmpdir=/scratch/cooperjam/temp -jar $GATK_HOME/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R {input.ref}  \
        -I {input.bam} \
        -targetIntervals {input.intervals} \
        -known {config[known_variants][dbsnp]} \
        -o {output} >& {log}
        """

rule gatk_recalibrate_info:
    input:
        ref=_get_ref,
        bam="{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.realigned.bam"
    output:
        temp("{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.recal_data.table")
    log:
        "{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.recalibrate_info.log"
    threads: 8
    shell:
        """
        java -Xmx8g -Djava.io.tmpdir=/scratch/cooperjam/temp -jar $GATK_HOME/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -R {input.ref}  \
        -L {config[INTERVAL]} \
        -nct 4 \
        -I {input.bam} \
        -knownSites {config[known_variants][dbsnp]} \
        -knownSites {config[known_variants][mills]} \
        -knownSites {config[known_variants][cosmic]} \
        -filterMBQ \
        -o {output} >& {log}
        """

rule gatk_post_recalibrate_info:
    input:
        ref=_get_ref,
        bam="{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.realigned.bam"
    output:
        "{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.post_recal_data.table"
    log:
        "{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.post_recal.log"
    threads: 8
    shell:
        """
        java -Xmx8g -Djava.io.tmpdir=/scratch/cooperjam/temp -jar $GATK_HOME/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -R {input.ref}  \
        -L {config[INTERVAL]} \
        -nct 8 \
        -I {input.bam} \
        -BQSR "{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.recal_data.table" \
        -knownSites {config[known_variants][dbsnp]} \
        -knownSites {config[known_variants][mills]} \
        -knownSites {config[known_variants][cosmic]} \
        -filterMBQ \
        -o {output} >& {log}
        """


rule gatk_analyzecovar:
    input:
        ref=_get_ref,
        before="{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.recal_data.table",
        after="{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.post_recal_data.table"
    output:
        "{FOLDERS}/mapping/{SAMPLE_RUN_LANE}_recalibration_plots.pdf"
    shell:
        """
        java -Xmx8g -Djava.io.tmpdir=/scratch/cooperjam/temp -jar $GATK_HOME/GenomeAnalysisTK.jar \
        -T AnalyzeCovariates \
        -R {input.ref}  \
        -before {input.before} \
        -after {input.after} \
        -plots {output}
        """

rule gatk_recalibrate_bam:
    input:
        ref=_get_ref,
        bam="{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.realigned.bam",
        BQSR="{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.post_recal_data.table"
    output:
        "{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.recalibrate.bam"
    log:
        "{FOLDERS}/mapping/{SAMPLE_RUN_LANE}.recalibrate.bam.log"
    threads: 4
    shell:
        """
        java -Xmx8g -Djava.io.tmpdir=/scratch/cooperjam/temp -jar $GATK_HOME/GenomeAnalysisTK.jar \
        -T PrintReads \
        -R {input.ref}  \
        -I {input.bam} \
        -BQSR {input.BQSR} \
        -o {output} \
        -nct 4 \
        -filterMBQ
        """
