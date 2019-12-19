#!/usr/bin/python3
import os
import subprocess
import argparse


def get_filename(PDIR, subdirectory, extension):
    filename = subprocess.run(f"ls {PDIR}/{subdirectory} | grep '.{extension}$'",
                              shell=True, capture_output=True).stdout.strip().decode('utf-8')
    return f"{PDIR}/{subdirectory}/{filename}"


def exists(file, extension):
    return os.path.exists(f"{file}.{extension}")


def run_shell(command):
    subprocess.run(command, shell=True)
    return 0


def main():
    """
    A parser for GATK Somatic Variant Discovery Pipeline.
    This parser will run GATK tools to complete the Pipeline
    and will generate (gbs) of data.

    Directory structure requiered:

    <Your-Project-Directory>
        * reference:
          fasta: reference genome (GRCh37 or GRCh38)

        * fastq:
          fastq: paired-end Illumina reads

        * variants:
          VCF: dbSNP VCF database

        * annotation:
          custom folder structure for annotation.
          example: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php

    """

    description = "Wrapper for GATK Somatic Variant Discovery"
    parser = argparse.ArgumentParser(description=description)
    variables = [
        {
            'name': 'project',
            'type': str,
            'help': 'Project Directory'
        }
    ]

    for variable in variables:
        parser.add_argument(
            variable['name'],
            type=variable['type'],
            help=variable['help']
        )

    args = parser.parse_args()
    PDIR = args.project

    if PDIR[-1] == '/':
        PDIR = PDIR[:-1]

    # [WARNING] hard-coded fastq-sequences
    PNAME = args.project.split('/')[0]
    REFGENOME = get_filename(PDIR, 'reference', 'fasta')
    VARIANTS = f"{PNAME}/variants/common_variants.vcf.gz"
    SOURCES = f"{PDIR}/annotation"
    F1 = f'{PDIR}/fastq/R1.fastq.gz'
    F2 = f'{PDIR}/fastq/R2.fastq.gz'

    # [INFO] Files to be created
    uBAM = f"{PDIR}/alignment/{PNAME}_unmapped.bam"
    gBAM = f"{PDIR}/alignment/{PNAME}_grouped.bam"
    rBAM = f"{PDIR}/alignment/{PNAME}_ready.bam"
    mBAM = f"{PDIR}/alignment/{PNAME}_mapped.bam"
    gBAM = f"{PDIR}/alignment/{PNAME}_qsorted.bam"
    scBAM = f"{PDIR}/alignment/{PNAME}_csorted.bam"
    mergedBAM = f'{PDIR}/alignment/{PNAME}_merged.bam'
    mutectVCF = f"{PDIR}/alignment/{PNAME}_mutect.vcf"
    filteredVCF = f"{PDIR}/alignment/{PNAME}_filtered.vcf"
    aVARIANTS = f"{PDIR}/alignment/{PNAME}_annotated.vcf"

    # Create missing directories if needed
    if not os.path.exists(f"{PDIR}/alignment"):
        make_dir = f'mkdir {PDIR}/alignment'
        run_shell(make_dir)

    # ReadFilesPipeline implementation
    # uBAM
    # [WARNING] This implementation requieres fixed names and should be modified for production.
    get_uBAM = f"gatk FastqToSam -F1 {F1} -F2 {F2} -O {uBAM} -SM {PNAME}"
    # run_shell(get_uBAM)

    # AddorReplaceReadGroups with default data
    set_RG = f"gatk AddOrReplaceReadGroups -I {uBAM} -O {gBAM} -LB trueSeqCancer -PL ILLUMINA -PU 01 -SM {PNAME} -SO queryname"
    # run_shell(set_RG)

    # Get IMG, DICT and indices for Reference Genome
    if not exists(REFGENOME, 'img'):
        get_img = f"gatk BwaMemIndexImageCreator -I {REFGENOME}"
        # run_shell(get_img)

    if not exists(f"{PDIR}/reference/{PNAME}", 'dict'):
        get_dict = f"gatk CreateSequenceDictionary -R {REFGENOME}"
        # run_shell(get_dict)

    if not exists(VARIANTS, 'tbi'):
        get_features = f"gatk IndexFeatureFile -I {VARIANTS}"
        # run_shell(get_features)

    if not exists(REFGENOME, 'fai'):
        get_fai = f"samtools faidx {REFGENOME}"
        # run_shell(get_fai)

    # Align BAM to Reference
    BAM_alignment = f"gatk BwaSpark \
        -I {gBAM} \
        -O {mBAM} \
        -R {REFGENOME}"
    # run_shell(BAM_alignment)

    # Merge BAM alignments
    merge_bams = f"gatk MergeBamAlignment \
    -O {mergedBAM} \
    -R {REFGENOME} \
    -UNMAPPED {gBAM} \
    -ALIGNED {mBAM}"
    # run_shell(merge_bams)

    # ReadsPipelineSpark
    run_pipeline = f"gatk ReadsPipelineSpark -I {mergedBAM} -O {VARIANTS} -R {REFGENOME} \
    -output-bam {rBAM} \
    -annotation AlleleFraction \
    -annotation DepthPerAlleleBySample \
    -annotation Coverage \
    -pairHMM AVX_LOGLESS_CACHING \
    --known-sites {VARIANTS}"
    # run_shell(run_pipeline)

    # Build Index for BAM visualization
    sort_bam_coordinate = f'gatk SortSam -I {rBAM} -O {scBAM} -SO coordinate'
    # run_shell(sort_bam_coordinate)

    build_bamidx = f'gatk BuildBamIndex -I {scBAM}'
    # run_shell(build_bamidx)

    # Mutect2
    run_mutect = f"gatk Mutect2 -I {rBAM} -O {mutectVCF} -R {REFGENOME}"
    # run_shell(run_mutect)

    # Filter Mutect Variants
    run_variant_filter = f"gatk FilterMutectCalls -O {filteredVCF} -R {REFGENOME} -V {mutectVCF}"
    # run_shell(run_variant_filter)

    # Annotate Variants
    funcotator = f"gatk Funcotator --data-sources-path {SOURCES} -O {aVARIANTS} --output-file-format VCF --ref-version hg19 -R {REFGENOME} -V {VARIANTS}"
    run_shell(funcotator)

    return 0


if __name__ == "__main__":
    main()
