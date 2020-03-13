#!/usr/bin/env bats

# Paths
DATE=$(date +"%Y.%m.%d")
DATA_DIR="data"
DATA_DIR=$(realpath ${DATA_DIR})
RUN_DIR="run/changeo-${DATE}"

# Run parameters
NPROC=2
OUTDIR=true
FAILED=true
GERMLINES="${HOME}/local/share/germlines/human/vdj"
V_GERMLINES="${HOME}/local/share/igblast/fasta/imgt_human_ig_v.fasta"
IGBLAST_DATA="${HOME}/local/share/igblast"
FORMAT="changeo"
#FORMAT="airr"

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


@test "AlignRecords-across" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AlignRecords.py across -d $DB --sf SEQUENCE_INPUT GERMLINE_IMGT_D_MASK --gf JUNCTION_LENGTH \
        --calls v d j --log $LOG --nproc $NPROC --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "AlignRecords-block" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AlignRecords.py block -d $DB --sf SEQUENCE_INPUT GERMLINE_IMGT_D_MASK --gf CLONE \
        --calls v d j --log $LOG --nproc $NPROC --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "AlignRecords-within" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run AlignRecords.py within -d $DB --sf SEQUENCE_INPUT GERMLINE_IMGT_D_MASK \
        --log $LOG --nproc $NPROC --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "AssignGenes-igblast-airr" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/db/HD13M-final_collapse-unique_atleast-2.fasta"
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
	READS="${DATA_DIR}/db/HD13M-final_collapse-unique_atleast-2.fasta"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
	OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run AssignGenes.py igblast -s $READS -b $IGBLAST_DATA --organism human --loci ig \
        --format blast $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "CreateGermlines" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_clone-pass.tab"
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
	DB="${DATA_DIR}/db/HD13M_clone-pass.tab"
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
	DB="${DATA_DIR}/db/HD13M_db-pass.tab"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run DefineClones.py -d $DB --model ham --dist 0.15 --mode gene --maxmiss 0 --act set \
        --log $LOG --format $FORMAT --nproc $NPROC $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "MakeDb-igblast" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	READS="${DATA_DIR}/db/HD13M-final_collapse-unique_atleast-2.fasta"
	ALIGNMENT="${DATA_DIR}/db/HD13M-final_collapse-unique_atleast-2.fmt7"
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
	READS="${DATA_DIR}/db/S43_atleast-2.fasta"
	ALIGNMENT="${DATA_DIR}/db/S43_atleast-2.txz"
	LOG="${RUN_DIR}/logs/${TEST}.log"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} ${FAILED})

    run MakeDb.py imgt -i $ALIGNMENT -s $READS -r $GERMLINES \
        --extended --log $LOG --format $FORMAT $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ConvertDb-airr" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py airr -d $DB $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ConvertDb-changeo" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/S43_airr.tsv"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py changeo -d $DB $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ConvertDb-baseline" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py baseline -d $DB --if SEQUENCE_ID --sf SEQUENCE_IMGT \
        --gf GERMLINE_IMGT_D_MASK --mf V_CALL J_CALL --cf CLONE $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ConvertDb-fasta" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py fasta -d $DB --if SEQUENCE_ID --sf SEQUENCE_IMGT \
        --mf V_CALL J_CALL $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ConvertDb-genbank-airr" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/S43_airr.tsv"
	SBT="${DATA_DIR}/db/template.sbt"
	YAML="${DATA_DIR}/db/genbank.yaml"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py genbank -d $DB --inf "IgBLAST:1.7.0" --organism "Homo sapiens" \
        --sex Male --tissue "Peripheral blood" --cf CPRIMER --nf DUPCOUNT --if MID \
        --asis-id --asn --sbt $SBT -y $YAML --format airr $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ConvertDb-genbank-changeo" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/S43_changeo.tab"
	SBT="${DATA_DIR}/db/template.sbt"
	YAML="${DATA_DIR}/db/genbank.yaml"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ConvertDb.py genbank -d $DB --inf "IgBLAST:1.7.0" --organism "Homo sapiens" \
        --sex Male --tissue "Peripheral blood" --cf CPRIMER --nf DUPCOUNT --if MID \
        --asis-id --asn --sbt $SBT -y $YAML --format changeo $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseDb-add" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py add -d $DB -f ADD_1 ADD_2 -u 1 2 $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseDb-delete" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py delete -d $DB -f PRCONS -u "IGHA|IGHG" --regex $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseDb-drop" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py drop -d $DB -f V_CALL J_CALL $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseDb-index" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py index -d $DB -f INDEX $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseDb-rename" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py rename -d $DB -f V_CALL J_CALL -k V_RENAME J_RENAME $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseDb-select" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py select -d $DB -f V_CALL J_CALL -u IGH --logic all --regex $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseDb-sort" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py sort -d $DB -f DUPCOUNT --num --descend $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseDb-split" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} true false)

    run ParseDb.py split -d $DB -f PRCONS $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}

@test "ParseDb-update" {
    TEST="${BATS_TEST_NUMBER}-${BATS_TEST_DESCRIPTION}"
	DB="${DATA_DIR}/db/HD13M_germ-pass.tab"
	CONSOLE="${RUN_DIR}/console/${TEST}.out"
    OUTPUT=$(get_output ${TEST} ${OUTDIR} false)

    run ParseDb.py update -d $DB -f PRCONS -u Human-IGHA Human-IGHG -t IGHA IGHG $OUTPUT

    echo "$output" > $CONSOLE
	[ "$status" -eq 0 ]
}


# AIRR
#AlignRecords across -d S43_IGBLAST_AIRR_db-pass.tsv --sf sequence junction_nt --calls v d j --format airr $OUTPUT $LOG $FAILED $NPROC
#AlignRecords block -d S43_IGBLAST_AIRR_db-pass.tsv --sf sequence junction_nt --gf junction_length --calls v d j --format airr \
#    $OUTPUT $LOG $FAILED $NPROC
#AlignRecords within -d S43_IGBLAST_AIRR_db-pass.tsv --sf sequence junction --format airr $OUTPUT $LOG $FAILED $NPROC
#CreateGermlines -d S43_IGBLAST_AIRR_db-pass_clone-pass.tsv -r germlines_gapped -g vonly dmask full regions --format airr $LOG $OUTPUT
#CreateGermlines -d S43_IGBLAST_AIRR_db-pass_clone-pass.tsv -r germlines_gapped -g vonly dmask full regions --cloned --format airr $LOG $OUTPUT
#DefineClones -d S43_IGBLAST_AIRR_db-pass.tsv --model ham --dist 0.15 --mode gene --maxmiss 5 --act first --format airr $FAILED $LOG $OUTPUT
#MakeDb igblast -i igblast/v1.7/S43_atleast-2.fmt7 -s igblast/v1.7/S43_atleast-2.fasta -r $PROJECT_DIR$/../germlines/IMGT/human/vdj \
#    --scores --partial --regions --format airr $OUTPUT
#MakeDb imgt -i imgt/S43_atleast-2.txz -s imgt/S43_atleast-2.fasta --regions --scores --junction --format airr $FAILED $OUTPUT
