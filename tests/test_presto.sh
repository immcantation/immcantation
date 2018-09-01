#!/usr/bin/env bats

# Run parametes
DATE=$(date +"%Y.%m.%d")
DATA_DIR="data"
DATA_DIR=$(readlink -f ${DATA_DIR})
RUN_DIR="run/presto-${DATE}"

# Run parameters
NPROC=2
OUTDIR=false
FAILED=true
GERMLINES="${HOME}/share/imgt/human/vdj"
V_GERMLINES="${HOME}/share/igblast/fasta/imgt_human_ig_v.fasta"
#GERMLINES="/usr/local/share/imgt/human/vdj"
#FORMAT=""
#FORMAT="--format airr"
#FASTA=true
CLUSTER="cd-hit-est"
ALIGNER="blastn"

# Create output parent
mkdir -p ${RUN_DIR}/logs ${RUN_DIR}/console ${RUN_DIR}/output
RUN_DIR=$(readlink -f ${RUN_DIR})

# Get output argument block
#${BATS_TEST_DESCRIPTION} is the description of the current test case.
get_output() {
    # $1 : output name for --outname or -o argument
    # $2 : whether to set --outdir/--outname (true) or -o (false)
    # $3 : whether to set --failed
    if [[ -e $2 && $3 ]]; then
        echo "--outdir ${RUN_DIR}/output --outname ${1} --failed"
    elif [[ $2 ]]; then
        echo "--outdir ${RUN_DIR}/output --outname ${1}"
    else
        echo "-o ${RUN_DIR}/output/${1}.fastq"
    fi
}

@test "AlignSets-muscle" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AlignSets.py muscle -s $READS --bf BARCODE --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "AlignSets-table" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	PRIMERS="${DATA_DIR}/primers/AbSeq_R2_TS.fasta"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AlignSets.py table -p $PRIMERS $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "AssemblePairs-align" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS_1="${DATA_DIR}/sequences/HD13M-R2_consensus-pass_pair-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R1_consensus-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AssemblePairs.py align -1 $READS_1 -2 $READS_2 --coord presto --rc tail --scanrev \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "AssemblePairs-join" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS_1="${DATA_DIR}/sequences/HD13M-R2_consensus-pass_pair-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R1_consensus-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AssemblePairs.py join -1 $READS_1 -2 $READS_2 --coord presto --rc tail --gap 10 \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "AssemblePairs-reference" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS_1="${DATA_DIR}/sequences/HD13M-R2_consensus-pass_pair-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R1_consensus-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AssemblePairs.py reference -1 $READS_1 -2 $READS_2 --coord presto --rc tail \
        -r $V_GERMLINES --minident 0.5 --evalue 1e-5 --maxhits 100 --aligner $ALIGNER \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "AssemblePairs-sequential" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS_1="${DATA_DIR}/sequences/HD13M-R2_consensus-pass_pair-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R1_consensus-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AssemblePairs.py sequential -1 $READS_1 -2 $READS_2 --coord presto --rc tail --scanrev \
        -r $V_GERMLINES --minident 0.5 --evalue 1e-5 --maxhits 100 --aligner $ALIGNER \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "BuildConsensus" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run BuildConsensus.py -s $READS -q 20 --cf PRIMER --act majority --maxerr 0.1 --maxgap 0.7 \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ClusterSets-all" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run ClusterSets.py all -s $READS --id 0.80 --prefix A --cluster $CLUSTER \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ClusterSets-barcode" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run ClusterSets.py barcode -s $READS -f BARCODE --id 0.80 --prefix B --cluster $CLUSTER \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ClusterSets-set" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run ClusterSets.py set -s $READS -f BARCODE --id 0.80 --prefix S \
        --start 25 --end 275 --cluster $CLUSTER \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "CollapseSeq" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M-final_total.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run CollapseSeq.py -s $READS -n 20 --cf CONSCOUNT PRCONS --act sum set --inner $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "FilterSeq-maskqual" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_L001_R1_001.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run FilterSeq.py maskqual -s $READS -q 20 --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "FilterSeq-quality" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_L001_R1_001.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run FilterSeq.py quality -s $READS -q 20 \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "FilterSeq-trimqual" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_L001_R1_001.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run FilterSeq.py trimqual -s $READS -q 20 --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "MaskPrimers-align" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_L001_R1_001.fastq"
	PRIMERS="${DATA_DIR}/primers/AbSeq_R1_Human_IG_Primers.fasta"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run MaskPrimers.py align -s $READS -p $PRIMERS --maxlen 50 --maxerror 0.20 --barcode --mode cut \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "MaskPrimers-extract" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_L001_R1_001.fastq"
	PRIMERS="${DATA_DIR}/primers/AbSeq_R1_Human_IG_Primers.fasta"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run MaskPrimers.py extract -s $READS --start 17 --len 20 --barcode --mode cut \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "MaskPrimers-score" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_L001_R1_001.fastq"
	PRIMERS="${DATA_DIR}/primers/AbSeq_R1_Human_IG_Primers.fasta"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run MaskPrimers.py score -s $READS -p $PRIMERS --start 17 --maxerror 0.20 --barcode --mode cut \
        --log $LOG --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "PairSeq" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS_1="${DATA_DIR}/sequences/HD13M-R1_primers-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R2_primers-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
	OUTPUT=$(get_output ${TEST} true ${FAILED})

    run PairSeq.py -1 $READS_1 -2 $READS_2 --coord illumina --1f PRIMER BARCODE --2f PRIMER \
        $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseHeaders-collapse" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_assemble-pass.fastq"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run ParseHeaders.py collapse -s $READS -f CONSCOUNT --act cat $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseHeaders-collapse" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_assemble-pass.fastq"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run ParseHeaders.py copy -s $READS -f ID PRCONS CONSCOUNT -k UMI C_CALL COUNT --act set set sum $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseHeaders-merge" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run ParseHeaders.py merge -s $READS -f PRIMER BARCODE -k UMI --act cat --delete $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseHeaders-rename" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_assemble-pass.fastq"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run ParseHeaders.py rename -s $READS -f CONSCOUNT PRCONS -k COUNT C_CALL $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseHeaders-table" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M-R1_primers-pass_pair-pass.fastq"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run ParseHeaders.py table -s $READS -f ID PRIMER BARCODE $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseLog" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	LOGFILE="${DATA_DIR}/logs/primers-1.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run ParseLog.py -l $LOGFILE -f ID PRIMER PRSTART BARCODE ERROR \
        $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "SplitSeq-count" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_assemble-pass.fastq"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} true false)

    run SplitSeq.py count -s $READS -n 300 $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "SplitSeq-group" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_assemble-pass.fastq"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} true false)

    run SplitSeq.py group -s $READS -f PRCONS $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "SplitSeq-sample" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M_assemble-pass.fastq"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} true false)

    run SplitSeq.py sample -s $READS -n 10 100 -f PRCONS $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}


@test "SplitSeq-samplepair" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS_1="${DATA_DIR}/sequences/HD13M-R1_primers-pass.fastq"
	READS_2="${DATA_DIR}/sequences/HD13M-R2_primers-pass.fastq"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} true false)

    run SplitSeq.py samplepair -1 $READS_1 -2 $READS_2 -n 10 100 $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "SplitSeq-select" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/sequences/HD13M-R1_consensus-pass_pair-pass.fastq"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run SplitSeq.py select -s $READS -f PRCONS -u Human-IGHA Human-IGHG $OUTPUT

    echo "$output" > $CONSOLE
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
#UnifyHeaders consensus -s join_test.fastq -f BARCODE -k SAMPLE --nproc 2 --log $LOG $OUTPUT
#UnifyHeaders delete -s join_test.fastq -f BARCODE -k SAMPLE --nproc 2 --log $LOG $OUTPUT
