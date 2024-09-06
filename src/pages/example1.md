---
title: Basic pipeline
layout: "@layouts/MarkdownPage.astro"
---

<div class="blg-summary example">
<h3>Basic pipeline</h3>

<p class="text-muted" >
    This example shows a simple Nextflow pipeline consisting of two Bash processes.
</p>

```groovy
#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input
params.greeting = "Hello World!"

/*
 * Redirect a string to a text file
 */
process sayHello {

    input:
    val x

    output:
    path 'output.txt'

    script:
    """
    echo '$x' > output.txt
    """
}

/*
 * Convert lowercase letters to uppercase letters
 */
process convertToUpper {

    input:
    path y

    output:
    stdout

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}

/*
 * Workflow definition
 */
workflow {

    // Creates channel using the Channel.of() channel factory
    greeting_ch = Channel.of(params.greeting)

    // Redirects a string to a text file
    sayHello(greeting_ch)

    // Concatenates a text file and transforms lowercase letters to uppercase letters
    convertToUpper(sayHello.out)

    // View convertToUpper output
    convertToUpper.out.view()
}
```

</div>

### Script synopsis

This example shows a simple Nextflow pipeline consisting of two Bash processes. The `sayHello` process takes a string as input and redirects it to an output text file. The `convertToUpper` process takes the output text file from `sayHello` as input, concatenates the text, and converts all of the lowercase letters to uppercase letters. The output from the `convertToUpper` process is then printed to screen.

### Try it

To try this pipeline:

1. Follow the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html#install-nextflow) to install Nextflow (if not already available).
2. Copy the script above and save it as `hello-world.nf`.
3. Launch the pipeline:

    nextflow run hello-world.nf

**NOTE**: To run this example with versions of Nextflow older than 22.04.0, you must include the `-dsl2` flag with `nextflow run`.
