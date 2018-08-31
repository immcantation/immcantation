#!/usr/bin/env bats

# Run parametes
DATE=$(date +"%Y.%m.%d")
DATA_DIR="data"
DATA_DIR=$(readlink -f ${DATA_DIR})
RUN_DIR="run/presto-${DATE}"

# Run parameters
NPROC=2
OUTDIR=false
GERMLINES="${HOME}/share/imgt/human/vdj"
V_GERMLINES="${HOME}/share/igblast/fasta/imgt_human_ig_v.fasta"
#GERMLINES="/usr/local/share/imgt/human/vdj"
#FORMAT=""
#FORMAT="--format airr"
#FASTA="--fasta"
CLUSTER="cd-hit-est"
ALIGNER="blastn"
#FAILED=""
FAILED="--failed"

# Create output parent
mkdir -p ${RUN_DIR}/logs ${RUN_DIR}/console ${RUN_DIR}/output
RUN_DIR=$(readlink -f ${RUN_DIR})

# Get output argument block
#${BATS_TEST_DESCRIPTION} is the description of the current test case.
get_output() {
    if $1; then
        echo "--outdir ${RUN_DIR}/output --outname ${2} ${FAILED}"
    else
        echo "-o ${RUN_DIR}/output/${2}.fastq"
    fi
}

@test "AlignSets-muscle" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run AlignSets.py muscle -s $READS --bf BARCODE --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "AlignSets-table" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	PRIMERS="${DATA_DIR}/primers/AbSeq_R2_TS.fasta"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run AlignSets.py table -p $PRIMERS $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "AssemblePairs-align" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS_1="${DATA_DIR}/sequences/HD13M-R2_consensus-pass_pair-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R1_consensus-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run AssemblePairs.py align -1 $READS_1 -2 $READS_2 --coord presto --rc tail --scanrev \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "AssemblePairs-join" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS_1="${DATA_DIR}/sequences/HD13M-R2_consensus-pass_pair-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R1_consensus-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run AssemblePairs.py join -1 $READS_1 -2 $READS_2 --coord presto --rc tail --gap 10 \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "AssemblePairs-reference" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS_1="${DATA_DIR}/sequences/HD13M-R2_consensus-pass_pair-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R1_consensus-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run AssemblePairs.py reference -1 $READS_1 -2 $READS_2 --coord presto --rc tail \
        -r $V_GERMLINES --minident 0.5 --evalue 1e-5 --maxhits 100 --aligner $ALIGNER \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "AssemblePairs-sequential" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS_1="${DATA_DIR}/sequences/HD13M-R2_consensus-pass_pair-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R1_consensus-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run AssemblePairs.py sequential -1 $READS_1 -2 $READS_2 --coord presto --rc tail --scanrev \
        -r $V_GERMLINES --minident 0.5 --evalue 1e-5 --maxhits 100 --aligner $ALIGNER \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "BuildConsensus" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run BuildConsensus.py -s $READS -q 20 --cf PRIMER --act majority --maxerr 0.1 --maxgap 0.7 \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "ClusterSets-all" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run ClusterSets.py all -s $READS --id 0.80 --prefix A --cluster $CLUSTER \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "ClusterSets-barcode" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run ClusterSets.py barcode -s $READS -f BARCODE --id 0.80 --prefix B --cluster $CLUSTER \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "ClusterSets-set" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run ClusterSets.py set -s $READS -f BARCODE --id 0.80 --prefix S \
        --start 25 --end 275 --cluster $CLUSTER \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "CollapseSeq" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/HD13M-final_total.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run CollapseSeq.py set -s $READS -n 20 --cf CONSCOUNT PRCONS --act sum set --inner \
        --log $LOG $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "FilterSeq-maskqual" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run FilterSeq.py maskqual -s $READS -q 20 \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "FilterSeq-quality" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run FilterSeq.py quality -s $READS -q 20 \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "FilterSeq-trimqual" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run FilterSeq.py trimqual -s $READS -q 20 \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "MaskPrimers-align" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq"
	PRIMERS="${DATA_DIR}/primers/AbSeq_R1_Human_IG_Primers.fasta"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run MaskPrimers.py align -s $READS -p $PRIMERS --maxlen 50 --maxerror 0.20 --barcode --mode cut \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "MaskPrimers-extract" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq"
	PRIMERS="${DATA_DIR}/primers/AbSeq_R1_Human_IG_Primers.fasta"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run MaskPrimers.py extract -s $READS --start 17 --len 20 --barcode --mode cut \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "MaskPrimers-score" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS="${DATA_DIR}/sequences/MG2v3_HD13M_BC13_AGTCAA_L001_R1_001.fastq"
	PRIMERS="${DATA_DIR}/primers/AbSeq_R1_Human_IG_Primers.fasta"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run MaskPrimers.py score -s $READS -p $PRIMERS --start 17 --maxerror 0.20 --barcode --mode cut \
        --log $LOG --nproc $NPROC $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "PairSeq" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	READS_1="${DATA_DIR}/sequences/HD13M-R1_primers-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R2_primers-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run PairSeq.py -1 $READS_1 -2 $READS_2 --coord illumina --1f PRIMER BARCODE --2f PRIMER \
        $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

@test "ParseLog" {
    TEST="${BATS_TEST_DESCRIPTION}-${BATS_TEST_NUMBER}"
	LOGFILE="${DATA_DIR}/logs/primers-1.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${OUTDIR} ${TEST})

    run ParseLog.py -l $LOGFILE -f ID PRIMER PRSTART BARCODE ERROR \
        $OUTPUT > $CONSOLE

	[ "$status" -eq 0 ]
}

# Presto
#ConvertHeaders 454 -s 454_headers.fasta $OUTPUT
#ConvertHeaders genbank -s genbank_headers.fasta $OUTPUT
#ConvertHeaders generic -s header_formats.fasta $OUTPUT
#ConvertHeaders illumina -s illumina_headers.fasta $OUTPUT
#ConvertHeaders imgt -s IMGT_Human_IGHC.fasta $OUTPUT
#ConvertHeaders migec -s migec_headers.fastq $OUTPUT
#ConvertHeaders sra -s sra_headers.fasta $OUTPUT
#EstimateError -s MS12_R1_primers-pass.fastq -n 100 $OUTPUT
#ParseHeaders collapse -s Assembled.fastq -f CONSCOUNT --act cat $OUTPUT
#ParseHeaders copy -s R2_100K_primers-mask_paired_align-pass_consensus-pass.fastq -f ID PRCOUNT PRIMER -k UID COUNT COUNT --act set cat set $OUTPUT
#ParseHeaders merge -s R2_100K_primers-mask_paired_align-pass_consensus-pass.fastq -f PRIMER PRCOUNT -k COUNT --delete --act cat $OUTPUT
#ParseHeaders rename -s R2_100K_primers-mask_paired_align-pass_consensus-pass.fastq -f PRCOUNT PRIMER -k COUNT NEWPRIMER $OUTPUT
#ParseHeaders table -s R2_100K_primers-mask_paired_align-pass_consensus-pass.fastq -f ID PRCOUNT PRIMER $OUTPUT
#SplitSeq count -s R1_100K_primers-mask_paired.fastq -n 20000 $OUTPUT
#SplitSeq group -s R1_100K_primers-mask_paired.fastq -f PRIMER $OUTPUT
#SplitSeq sample -s R1_100K_primers-mask_paired.fastq -n 3 10 -f BARCODE $OUTPUT
#SplitSeq samplepair -1 R1_100K_primers-mask_paired.fastq -2 R2_100K_primers-mask_paired.fastq -n 10 2000 -f BARCODE $OUTPUT
#SplitSeq select -s head.fasta -f PRIMER -u P3 P2 --not $OUTPUT
#UnifyHeaders consensus -s join_test.fastq -f BARCODE -k SAMPLE --nproc 2 --log $LOG $OUTPUT
#UnifyHeaders delete -s join_test.fastq -f BARCODE -k SAMPLE --nproc 2 --log $LOG $OUTPUT
