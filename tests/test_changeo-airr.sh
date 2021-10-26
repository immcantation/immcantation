#!/usr/bin/env bats

# Paths
DATE=$(date +"%Y.%m.%d")
DATA_DIR="data"
DATA_DIR=$(realpath ${DATA_DIR})
RUN_DIR="run/changeo-airr-${DATE}"

# Run parameters
NPROC=2
OUTDIR=true
FAILED=true
GERMLINES="${HOME}/local/share/germlines/human/vdj"
V_GERMLINES="${HOME}/local/share/igblast/fasta/imgt_human_ig_v.fasta"
IGBLAST_DATA="${HOME}/local/share/igblast"
FORMAT="airr"

# Create output parent
mkdir -p ${RUN_DIR}/logs ${RUN_DIR}/console ${RUN_DIR}/output
RUN_DIR=$(realpath ${RUN_DIR})

# Get output arguments
get_output() {
    # $1 : output name for --outname or -o argument
    # $2 : whether to set --outdir/--outname (true) or -o (false)
    # $3 : whether to set --failed
    if $2 && $3; then
        echo "--outdir ${RUN_DIR}/output --outname ${1} --failed"
    elif $2; then
        echo "--outdir ${RUN_DIR}/output --outname ${1}"
    else
        echo "-o ${RUN_DIR}/output/${1}.txt"
    fi
}

@test "AssignGenes-igblast-airr" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    READS="${DATA_DIR}/db/HD13M-VH.fasta"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run AssignGenes.py igblast -s $READS -b $IGBLAST_DATA --organism human --loci ig \
        --format airr $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "AssignGenes-igblast-blast" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    READS="${DATA_DIR}/db/HD13M-VH.fasta"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run AssignGenes.py igblast -s $READS -b $IGBLAST_DATA --organism human --loci ig \
        --format blast $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "MakeDb-igblast" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    READS="${DATA_DIR}/db/HD13M-VH.fasta"
    ALIGNMENT="${DATA_DIR}/db/HD13M-VH.fmt7"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run MakeDb.py igblast -i $ALIGNMENT -s $READS -r $GERMLINES \
        --extended --log $LOG --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "MakeDb-ihmm" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    READS="${DATA_DIR}/db/sim_ihmm.fasta"
    ALIGNMENT="${DATA_DIR}/db/sim_ihmm.txt"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run MakeDb.py ihmm -i $ALIGNMENT -s $READS -r $GERMLINES \
        --extended --log $LOG --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "MakeDb-imgt" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    READS="${DATA_DIR}/db/HD13M-VH.fasta"
    ALIGNMENT="${DATA_DIR}/db/HD13M-VH.txz"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run MakeDb.py imgt -i $ALIGNMENT -s $READS -r $GERMLINES \
        --extended --log $LOG --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "CreateGermlines" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run CreateGermlines.py -d $DB -r $GERMLINES -g vonly dmask full regions \
        --log $LOG --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "CreateGermlines-cloned" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run CreateGermlines.py -d $DB -r $GERMLINES -g vonly dmask full regions --cloned \
        --log $LOG --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "DefineClones" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run DefineClones.py -d $DB --model ham --dist 0.15 --mode gene --maxmiss 0 --act set \
        --log $LOG --format $FORMAT --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "AlignRecords-across" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_germ-pass.tsv"
    SEQ_FIELD="sequence_alignment germline_alignment_d_mask"
    GROUP_FIELD="junction_length"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AlignRecords.py across -d $DB --sf $SEQ_FIELD --gf $GROUP_FIELD \
        --calls v d j --log $LOG --nproc $NPROC --format $FORMAT $OUTPUT
        
    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "AlignRecords-block" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_germ-pass.tsv"
    SEQ_FIELD="sequence_alignment germline_alignment_d_mask"
    GROUP_FIELD="clone_id"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AlignRecords.py block -d $DB --sf $SEQ_FIELD --gf $GROUP_FIELD \
        --calls v d j --log $LOG --nproc $NPROC --format $FORMAT $OUTPUT
        
    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "AlignRecords-within" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_germ-pass.tsv"
    SEQ_FIELD="sequence germline_alignment_d_mask"
    GROUP_FIELD="clone_id"
    LOG="${RUN_DIR}/logs/${TEST}.log"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AlignRecords.py within -d $DB --sf $SEQ_FIELD \
        --log $LOG --nproc $NPROC --format $FORMAT $OUTPUT
        
    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ConvertDb-airr" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_germ-pass.tsv"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py airr -d $DB $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ConvertDb-changeo" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_germ-pass.tsv"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py changeo -d $DB $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ConvertDb-baseline" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_germ-pass.tsv"
    ID_FIELD="sequence_id"
    SEQ_FIELD="sequence_alignment"
    GERM_FIELD="germline_alignment_d_mask"
    META_FIELD="v_call j_call"
    CLONE_FIELD="clone_id"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py baseline -d $DB --if $ID_FIELD --sf $SEQ_FIELD \
        --gf $GERM_FIELD --mf $META_FIELD --cf $CLONE_FIELD $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ConvertDb-fasta" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    ID_FIELD="sequence_id"
    SEQ_FIELD="sequence_alignment"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py fasta -d $DB --if $ID_FIELD --sf $SEQ_FIELD \
        --mf V_CALL J_CALL $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ConvertDb-genbank" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_germ-pass.tsv"
    SBT="${DATA_DIR}/db/template.sbt"
    YAML="${DATA_DIR}/db/genbank.yaml"
    CREGION_FIELD="c_call"
    COUNT_FIELD="duplicate_count"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py genbank -d $DB --inf "IgBLAST:1.14.0" --organism "Homo sapiens" \
        --sex Male --tissue "Peripheral blood" --cf $CREGION_FIELD --nf $COUNT_FIELD \
        --asis-id --asn --sbt $SBT -y $YAML --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ParseDb-add" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py add -d $DB -f add_1 add_2 -u 1 2 $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ParseDb-delete" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    PARSE_FIELD="c_call"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py delete -d $DB -f $PARSE_FIELD -u "IGHA|IGHG" --regex $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ParseDb-drop" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    PARSE_FIELD="v_call j_call"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py drop -d $DB -f $PARSE_FIELD $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ParseDb-index" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py index -d $DB -f index $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ParseDb-rename" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    PARSE_FIELD="v_call j_call"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py rename -d $DB -f $PARSE_FIELD -k rename_1 rename_2 $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ParseDb-select" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    PARSE_FIELD="v_call j_call"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py select -d $DB -f $PARSE_FIELD -u IGH --logic all --regex $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ParseDb-sort" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    PARSE_FIELD="duplicate_count"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py sort -d $DB -f $PARSE_FIELD --num --descend $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ParseDb-split" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    PARSE_FIELD="c_call"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} true false)

    run ParseDb.py split -d $DB -f $PARSE_FIELD $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}

@test "ParseDb-update" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
    DB="${DATA_DIR}/db/HD13M-VH_clone-pass.tsv"
    PARSE_FIELD="c_call"
    CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py update -d $DB -f $PARSE_FIELD -u IGHA IGHG -t IgA IgG $OUTPUT

    echo "$output" > $CONSOLE
    [ "$status" -eq 0 ]
}
