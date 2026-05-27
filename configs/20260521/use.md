nextflow run ~/work/pipelines/vbcf-pipseq --input $INPUTCSV --index=$INDEX --outdir $OUTDIR --project_id $PID --pipspeaker $PIPSPEAKERYAML \
                                                 --umilen 12 --umi_offset 10 --cblen 24 --clip3pNbases 45 \
                                                 --cellfilter $CELLFILTER --cellparams " $CELLPARAMS " \
                                                 --vitessce --cloupe \
                                                 --tar \
                                                 -resume -profile cbe

script:
    def pipspeak = params.pipspeak
    def (forward, reverse) = reads.collate(2).transpose()
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_pipspeak"
    def r1m = "${meta.id}_R1.fq.gz"
    def r2m = "${meta.id}_R2.fq.gz"

    def to_merge = true //reverse.length > 1

    def merge = to_merge ? 
    """
    cat ${forward.join( " " )} > ${r1m}
    cat ${reverse.join( " " )} > ${r2m}
    """ : ""

 
    def r1 = to_merge ? r1m : reads[0]
    def r2 = to_merge ? r2m : reads[1]
    def offset = params.offset ?: 8
    def umi_offset = params.umi_offset ?: 0



```bash
pipspeak \
        --loglevel info \
        --config $config \
        --prefix $prefix \
        --r1 ${r1} \
        --r2 ${r2} \
        --offset ${offset} \
        --umi-len ${params.umilen} \
        --umi-offset ${umi_offset} \
        --threads $task.cpus
```

  def cb_len = params.cblen 
    def umi_len = params.umilen
    def clip3pNbases = params.clip3pNbases //was hardcoded to 45 because of PE150RL -> 105!!!

    if( meta.protocol == "10x" ) {
        cb_len = 16
        umi_len = 12
    }


    def cellfilter = params.cellfilter ?: "EmptyDrops_CR"
    def cellparams = params.cellparams ?: "5000 0.99 10 45000 90000 500 0.01 20000 0.01 10000"

    //pipseq v2
    // CBstart: 1 , CBlen: 24 , UMIstart: 24 , UMIlen: 12  cb_len: 24

    //pipseq v3
    // CBstart: 1 , CBlen: 28 , UMIstart: 29 , UMIlen: 12
    //SRR19180490: length R1: 55, length R2: 66 -> no trimming

    //fluent v4:
    // CBstart: 1 , CBlen: 28 , UMIstart: 29 , UMIlen: 12
    //pipseq 20240422
    // CBstart: 1 , CBlen: 30 , UMIstart: 31 , UMIlen: 8

    //pipseq 20241008 (removed spacer 3 and barcode 4)
    // CBstart: 1 , CBlen: 24 , UMIstart: 25 , UMIlen: 8

    //pipseq 20241101 (removed spacer 3 and barcode 4)
    // CBstart: 1 , CBlen: 24 , UMIstart: 25 , UMIlen: 12    

    //now clipping in cutadapt which makes much more sense 20250724
    //clip3pNbases : 45 from 3' end => 105 for PE 150
    // --clip3pNbases ${clip3pNbases} \\ was shown to be necessary for PIPseq batch R19026

    """
    STAR \\
        --genomeDir $index \\
        --readFilesCommand zcat \\
        --outFileNamePrefix $prefix. \\
        --readFilesIn ${reads[1]} ${reads[0]} \\
        --runThreadN $task.cpus \\
        --clip3pNbases ${clip3pNbases} \\
        --soloCBstart 1 \\
        --soloCBlen ${cb_len} \\
        --soloUMIstart ${cb_len + 1} \\
        --soloUMIlen ${umi_len}  \\
        --soloType CB_UMI_Simple \\
        --soloCBwhitelist None \\
        --soloStrand Forward \\
        --soloCellFilter ${cellfilter} ${cellparams.replaceAll("'","")} \\
        --soloUMIdedup 1MM_CR \\
        --soloUMIfiltering  MultiGeneUMI_CR \\
        --outSAMattributes CB UB \\
        --outSAMtype BAM SortedByCoordinate \\
        --soloFeatures Gene GeneFull SJ \\
