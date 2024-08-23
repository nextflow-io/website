---
title: Basic pipeline
layout: "@layouts/MarkdownPage.astro"
---

<div class="blg-summary example">
<h3>Basic pipeline</h3>

<p class="text-muted" >
    This example shows a simple Nextflow pipeline consisting of two Bash processes, where the output from the first process is used as input for the second process.
</p>

```groovy
#!/usr/bin/env nextflow

params.greeting = "Hello World!"

/*
 * Redirect a string to a text file
 */
process sayHello {

    input:
    val x

    output:
    path 'output.txt'

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

    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}

/*
 * Workflow definition
 */
workflow {
    sayHello(params.greeting)
        | convertToUpper
        | view
}
```

</div>

### Try it

To try this pipeline:

1. Follow the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html#install-nextflow) to install Nextflow.
2. Copy the script above and save it as `hello-world.nf`.
3. Launch the pipeline:

    nextflow run hell-world.nf

4. Launch the pipeline again with a custom greeting:

    nextflow run hello-world.nf --greeting "Bonjour le monde!"

### Script synopsis

- **Line 1** Declares Nextflow as the interpreter.

- **Line 3**: Declares a pipeline parameter named `greeting` that is initialized with the value `"Hello World!"`.

- **Lines 8-19**: Declares a process named `sayHello` that redirects a string to a text file.

  - **Line 10**: Opens the input declaration block.

  - **Line 11**: Defines the process input `x`.

  - **Line 13**: Opens the output declaration block.

  - **Line 14**: Defines the process output `'output.txt'`.

  - **Lines 16-18**: Defines a script that redirects the string `x` to a text file named `output.txt`.

- **Lines 24-35**: Declares a process named `convertToUpper` that concatenates a file and transforms all of the lowercase letters to uppercase letters.

  - **Line 26**: Opens the input declaration block.

  - **Line 27**: Defines the process input `y`.

  - **Line 29**: Opens the output declaration block.

  - **Line 30**: Defines standard output (`stdout`) as the output.

  - **Lines 32-34**: Defines a script that concatenates the variable `y` and transforms all of the lowercase letters to uppercase letters.

- **Lines 40-44**: Declares the workflow that connects everything together!

  - **Line 41**: Passes the string specified by `params.greeting` to the `sayHello` process.

  - **Line 42**: Passes of the output from `sayHello` to the `convertToUpper` process.

  - **Line 43**: Prints the standard output from `convertToUpper`.
