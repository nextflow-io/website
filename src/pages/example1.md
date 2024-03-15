---
title: Basic pipeline
layout: "@layouts/MarkdownPage.astro"
---

<div class="blg-summary example">
<h3>Basic pipeline</h3>

<p class="text-muted" >
    This example shows how to write a pipeline with two simple Bash processes, so that the results produced by the first process are consumed by the second process.
</p>

```groovy
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
```

</div>

### Synopsis

- **Line 1** The script starts with a shebang declaration. This allows you to launch your pipeline just like any other Bash script.

- **Line 3**: Declares a pipeline parameter named `params.in` that is initialized with the value `$HOME/sample.fa`. This value can be overridden when launching the pipeline, by simply adding the option `--in <value>` to the script command line.

- **Lines 8-19**: The process that splits the provided file.

  - **Line 10**: Opens the input declaration block. The lines following this clause are interpreted as input definitions.

  - **Line 11**: Declares the process input file, which will be named `input.fa` in the process script.

  - **Line 13**: Opens the output declaration block. The lines following this clause are interpreted as output declarations.

  - **Line 14**: Files whose names match the pattern `seq_*` are declared as the output of this process.

  - **Lines 16-18**: The actual script executed by the process to split the input file.

- **Lines 24-35**: The second process, which receives the splits produced by the
  previous process and reverses their content.

  - **Line 26**: Opens the input declaration block. Lines following this clause are
    interpreted as input declarations.

  - **Line 27**: Defines the process input file.

  - **Line 29**: Opens the output declaration block. Lines following this clause are
    interpreted as output declarations.

  - **Line 30**: The standard output of the executed script is declared as the process
    output.

  - **Lines 32-34**: The actual script executed by the process to reverse the content of the input files.

- **Lines 40-44**: The workflow that connects everything together!

  - **Line 41**: First, the input file specified by `params.in` is passed to the `splitSequences` process.

  - **Line 42**: The outputs of `splitSequences` are passed as inputs to the `reverse` process, which processes each split file in parallel.

  - **Line 43**: Finally, each output emitted by `reverse` is printed.
