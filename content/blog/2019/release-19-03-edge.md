title=Release 19.03: The Sequence Read Archive & more!
date=2019-03-18
type=post
tags=nextflow,release
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

It's time for the monthly Nextflow release for March, edge version 19.03. This is another great release with some cool new features, bug fixes and improvements.


### SRA channel factory
This sees the introduction of the long-awaited sequence read archive (SRA) channel factory.
The SRA is a key public repository for sequencing data and run in coordination between The National Center for
Biotechnology Information (NCBI), The European Bioinformatics Institute (EBI) and the DNA Data Bank of Japan (DDBJ).

This feature originates all the way back in [2015](https://github.com/nextflow-io/nextflow/issues/89) and was worked on during a 2018 Nextflow hackathon. It was brought to fore again with the release of Phil Ewels [SRA Explorer](https://ewels.github.io/sra-explorer/). The SRA channel factory allows users to pull read data in FASTQ format directly from SRA by referencing a study, accession ID or even a keyword. It works in a similar way to [`fromFilePairs`](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs), returning a sample id and files (or pairs of files) for each sample.

The code snippet below creates a channel containing 24 samples from a chromatin dynamics study and runs FASTQC on the resulting files.

```
Channel
    .fromSRA('SRP043510')
    .set {reads}

process fastqc {
    input:
    set sample_id, file(reads_file) from reads

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}
```

See the [documentation](https://www.nextflow.io/docs/edge/channel.html#fromsra) for more details. When combined with downstream processes, you can quickly open a firehose of data on your workflow!

### Edge release
Note that this is a monthly edge release. To use it simply execute the following command prior to running Nextflow:
```
export NXF_VER=19.03.0-edge
```

### If you need help…
…please don’t hesitate to use our very active [Gitter](https://gitter.im/nextflow-io/nextflow) channel or create a thread in the [forum](https://groups.google.com/forum/#!forum/nextflow).

### Reporting Issues
Experiencing issues introduced by this release? Please report them in our [issue tracker](https://github.com/nextflow-io/nextflow/issues). Make sure to fill in the fields of the issue template.

### Contributions
Special thanks to the contributors of this release:
* Paolo Di Tommaso - [pditommaso](https://github.com/pditommaso)
* Jon Haitz Legarreta Gorroño - [jhlegarreta](https://github.com/jhlegarreta)
* Toni Hermoso Pulido - [toniher](https://github.com/toniher)
* Lukas Jelonek - [lukasjelonek](https://github.com/lukasjelonek)
* Jonathan Leitschuh - [JLLeitschuh](https://github.com/JLLeitschuh)
* [phue](https://github.com/phue)
* Philippe Hupé [phupe](https://github.com/phupe)
* Kevin Sayers - [KevinSayers](https://github.com/KevinSayers)
* Akira Sekiguchi - [pachiras](https://github.com/pachiras)


### Complete changes
- Fix Nextflow hangs submitting jobs to AWS batch #1024
- Fix process builder incomplete output [2fe1052c]
- Fix Grid executor reports invalid queue status #1045
- Fix Script execute permission is lost in container #1060
- Fix K8s serviceAccount is not honoured #1049
- Fix K8s kuberun login path #1072
- Fix K8s imagePullSecret and imagePullPolicy #1062
- Fix Google Storage docs #1023
- Fix Env variable NXF_CONDA_CACHEDIR is ignored #1051
- Fix failing task due to legacy sleep command [3e150b56]
- Fix SplitText operator should accept a closure parameter #1021
- Add Channel.fromSRA factory method #1070
- Add voluntary/involuntary context switches to metrics #1047
- Add noHttps option to singularity config #1041)
- Add docker-daemon Singularity support #1043 [dfef1391]
- Use peak_vmem and peak_rss as default output in the trace file instead of rss and vmem #1020
- Improve ansi log rendering #996 [33038a18]

### Breaking changes:
None known.
