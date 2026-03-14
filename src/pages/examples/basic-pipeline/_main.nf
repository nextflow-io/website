#!/usr/bin/env nextflow


params.in = "${baseDir}/data/sample.fa"


process splitSequences {
    input:
    path 'input.fa'

    output:
    path 'seq_*'

    script:
    """
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    """
}


process reverse {
    input:
    path x

    output:
    stdout

    script:
    """
    cat ${x} | rev
    """
}


workflow {

    splitSequences(params.in)
        | reverse
        | view
}
