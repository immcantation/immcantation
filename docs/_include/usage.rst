.. Start preprocess-phix

Usage: preprocess-phix [OPTIONS]
  -s   FASTQ sequence file.
  -r   Directory containing phiX174 reference db.
       Defaults to /usr/local/share/phix.
  -n   Sample identifier which will be used as the output file prefix.
       Defaults to a truncated version of the input filename.
  -o  Output directory. Will be created if it does not exist.
      Defaults to a directory matching the sample identifier in the current working directory.
  -p   Number of subprocesses for multiprocessing tools.
       Defaults to the available cores.
  -h   This message.

.. End preprocess-phix

.. Start presto-abseq

Usage: presto-abseq [OPTIONS]
  -1  Read 1 FASTQ sequence file.
      Sequence beginning with the C-region or J-segment).
  -2  Read 2 FASTQ sequence file.
      Sequence beginning with the leader or V-segment).
  -j  Read 1 FASTA primer sequences.
      Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R1_Human_IG_Primers.fasta.
  -v  Read 2 FASTA primer or template switch sequences.
      Defaults to /usr/local/share/protocols/AbSeq/AbSeq_R2_TS.fasta.
  -c  C-region FASTA sequences for the C-region internal to the primer.
      If unspecified internal C-region alignment is not performed.
  -r  V-segment reference file.
      Defaults to /usr/local/share/igblast/fasta/imgt_human_ig_v.fasta.
  -y  YAML file providing description fields for report generation.
  -n  Sample identifier which will be used as the output file prefix.
      Defaults to a truncated version of the read 1 filename.
  -o  Output directory. Will be created if it does not exist.
      Defaults to a directory matching the sample identifier in the current working directory.
  -x  The mate-pair coordinate format of the raw data.
      Defaults to illumina.
  -p  Number of subprocesses for multiprocessing tools.
      Defaults to the available cores.
  -h  This message.

.. End presto-abseq

.. Start presto-clontech

Usage: presto-clontech [OPTIONS]
  -1  Read 1 FASTQ sequence file.
      Sequence beginning with the C-region.
  -2  Read 2 FASTQ sequence file.
      Sequence beginning with the leader.
  -j  C-region reference sequences (reverse complemented).
      Defaults to /usr/local/share/protocols/Universal/Mouse_IG_CRegion_RC.fasta.
  -r  V-segment reference file.
      Defaults to /usr/local/share/igblast/fasta/imgt_mouse_ig_v.fasta.
  -y  YAML file providing description fields for report generation.
  -n  Sample identifier which will be used as the output file prefix.
      Defaults to a truncated version of the read 1 filename.
  -o  Output directory. Will be created if it does not exist.
      Defaults to a directory matching the sample identifier in the current working directory.
  -x  The mate-pair coordinate format of the raw data.
      Defaults to illumina.
  -p  Number of subprocesses for multiprocessing tools.
      Defaults to the available cores.
  -h  This message.

.. End presto-clontech

.. Start changeo-10x

Usage: changeo-10x [OPTIONS]
  -s  FASTA or FASTQ sequence file.
  -a  10x Genomics cellranger-vdj contig annotation CSV file.
      Must corresponding with the FASTA/FASTQ input file (all, filtered or consensus).
  -r  Directory containing IMGT-gapped reference germlines.
      Defaults to /usr/local/share/germlines/imgt/human/vdj when species is human.
      Defaults to /usr/local/share/germlines/imgt/mouse/vdj when species is mouse.
  -g  Species name. One of human or mouse. Defaults to human.
  -t  Receptor type. One of ig or tr. Defaults to ig.
  -x  Distance threshold for clonal assignment. Specify "auto" for automatic detection.
      If unspecified, clonal assignment is not performed.
  -m  Distance model for clonal assignment.
      Defaults to the nucleotide Hamming distance model (ham).
  -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
      and auxillary_data directories. Defaults to /usr/local/share/igblast.
  -n  Sample identifier which will be used as the output file prefix.
      Defaults to a truncated version of the sequence filename.
  -o  Output directory. Will be created if it does not exist.
      Defaults to a directory matching the sample identifier in the current working directory.
  -f  Output format. One of changeo or airr. Defaults to airr.
  -p  Number of subprocesses for multiprocessing tools.
      Defaults to the available cores.
  -i  Specify to allow partial alignments.
  -h  This message.

.. End changeo-10x

.. Start changeo-igblast

Usage: changeo-igblast [OPTIONS]
  -s  FASTA or FASTQ sequence file.
  -r  Directory containing IMGT-gapped reference germlines.
      Defaults to /usr/local/share/germlines/imgt/human/vdj when species is human.
      Defaults to /usr/local/share/germlines/imgt/mouse/vdj when species is mouse.
  -g  Species name. One of human or mouse. Defaults to human.
  -t  Receptor type. One of ig or tr. Defaults to ig.
  -b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
      and auxillary_data directories. Defaults to /usr/local/share/igblast.
  -n  Sample identifier which will be used as the output file prefix.
      Defaults to a truncated version of the sequence filename.
  -o  Output directory. Will be created if it does not exist.
      Defaults to a directory matching the sample identifier in the current working directory.
  -f  Output format. One of airr (default) or changeo.
  -p  Number of subprocesses for multiprocessing tools.
      Defaults to the available cores.
  -k  Specify to filter the output to only productive/functional sequences.
  -i  Specify to allow partial alignments.
  -h  This message.

.. End changeo-igblast

.. Start changeo-clone

Usage: changeo-clone [OPTIONS]
  -d  Change-O formatted TSV (TAB) file.
  -x  Distance threshold for clonal assignment.
  -m  Distance model for clonal assignment.
      Defaults to the nucleotide Hamming distance model (ham).
  -r  Directory containing IMGT-gapped reference germlines.
      Defaults to /usr/local/share/germlines/imgt/human/vdj.
  -n  Sample identifier which will be used as the output file prefix.
      Defaults to a truncated version of the input filename.
  -o  Output directory. Will be created if it does not exist.
      Defaults to a directory matching the sample identifier in the current working directory.
  -f  Output format. One of airr (default) or changeo.
  -p  Number of subprocesses for multiprocessing tools.
      Defaults to the available cores.
  -a  Specify to clone the full data set.
      By default the data will be filtering to only productive/functional sequences.
  -h  This message.

.. End changeo-clone

.. Start shazam-threshold

Usage: shazam-threshold [options]
	-d DB, --db=DB
		Tabulated data file, in Change-O (TAB) or AIRR format (TSV).
	-m METHOD, --method=METHOD
		Threshold inferrence to use. One of gmm, density, or none. 
		If none, the distance-to-nearest distribution is plotted without threshold detection. 
		Defaults to density.
	-n NAME, --name=NAME
		Sample name or run identifier which will be used as the output file prefix. 
		Defaults to a truncated version of the input filename.
	-o OUTDIR, --outdir=OUTDIR
		Output directory. Will be created if it does not exist. 
		Defaults to the current working directory.
	-f FORMAT, --format=FORMAT
		File format. One of 'airr' (default) or 'changeo'.
	-p NPROC, --nproc=NPROC
		Number of subprocesses for multiprocessing tools. 
		Defaults to the available processing units.
	--model=MODEL
		Model to use for the gmm model. 
		One of gamma-gamma, gamma-norm, norm-norm or norm-gamma. 
		Defaults to gamma-gamma.
	--subsample=SUBSAMPLE
		Number of distances to downsample the data to before threshold calculation. 
		By default, subsampling is not performed.
	--repeats=REPEATS
		Number of times to recalculate. 
		Defaults to 1.
	-h, --help
		Show this help message and exit


.. End shazam-threshold

.. Start tigger-genotype

Usage: tigger-genotype [options]
	-d DB, --db=DB
		Change-O formatted TSV (TAB) file.
	-r REF, --ref=REF
		FASTA file containing IMGT-gapped V segment reference germlines. 
		Defaults to /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta.
	-v VFIELD, --vfield=VFIELD
		Name of the output field containing genotyped V assignments. 
		Defaults to V_CALL_GENOTYPED.
	-x MINSEQ, --minseq=MINSEQ
		Minimum number of sequences in the mutation/coordinate range. 
		Samples with insufficient sequences will be excluded. 
		Defaults to 50.
	-y MINGERM, --mingerm=MINGERM
		Minimum number of sequences required to analyze a germline allele. 
		Defaults to 200.
	-n NAME, --name=NAME
		Sample name or run identifier which will be used as the output file prefix. 
		Defaults to a truncated version of the input filename.
	-o OUTDIR, --outdir=OUTDIR
		Output directory. Will be created if it does not exist. 
		Defaults to the current working directory.
	-f FORMAT, --format=FORMAT
		File format. One of 'airr' (default) or 'changeo'.
	-p NPROC, --nproc=NPROC
		Number of subprocesses for multiprocessing tools. 
		Defaults to the available processing units.
	-h, --help
		Show this help message and exit


.. End tigger-genotype

