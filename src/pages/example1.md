---
title: Basic pipeline
layout: @layouts/Page.astro
---

<div class="blg-summary example">
<h3><a href="javascript:void(0)">Basic pipeline</a></h3>

<p class="text-muted" >
    This example shows how to write a pipeline with two simple Bash processes, so that the results produced by the first process are consumed by the second process.
</p>


<script type="syntaxhighlighter" class="brush: groovy">
<![CDATA[
#!/usr/bin/env nextflow

params.in = "$baseDir/data/sample.fa"

/*
 * Split a fasta file into multiple files
 */
process splitSequences {

    input:
    path 'input.fa'

    output:
    path 'seq_*'

    """
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    """
}

/*
 * Reverse the sequences
 */
process reverse {

    input:
    path x

    output:
    stdout

    """
    cat $x | rev
    """
}

/*
 * Define the workflow
 */
workflow {
    splitSequences(params.in) \
      | reverse \
      | view
}
]]>
</script>
</div>


### Synopsis

* __Line 1__ The script starts with a shebang declaration. This allows you to launch your pipeline just like any other Bash script.

* __Line 3__: Declares a pipeline parameter named `params.in` that is initialized with the value `$HOME/sample.fa`. This value can be overridden when launching the pipeline, by simply adding the option `--in <value>` to the script command line.

* __Lines 8-19__: The process that splits the provided file.

  * __Line 10__: Opens the input declaration block. The lines following this clause are interpreted as input definitions.

  * __Line 11__: Declares the process input file, which will be named `input.fa` in the process script.

  * __Line 13__: Opens the output declaration block. The lines following this clause are interpreted as output declarations.

  * __Line 14__: Files whose names match the pattern `seq_*` are declared as the output of this process.

  * __Lines 16-18__: The actual script executed by the process to split the input file.

* __Lines 24-35__: The second process, which receives the splits produced by the
previous process and reverses their content.

  * __Line 26__: Opens the input declaration block. Lines following this clause are
interpreted as input declarations.

  * __Line 27__: Defines the process input file.

  * __Line 29__: Opens the output declaration block. Lines following this clause are
interpreted as output declarations.

  * __Line 30__: The standard output of the executed script is declared as the process
output.

  * __Lines 32-34__: The actual script executed by the process to reverse the content of the input files.

* __Lines 40-44__: The workflow that connects everything together!

  * __Line 41__: First, the input file specified by `params.in` is passed to the `splitSequences` process.

  * __Line 42__: The outputs of `splitSequences` are passed as inputs to the `reverse` process, which processes each split file in parallel.

  * __Line 43__: Finally, each output emitted by `reverse` is printed.
