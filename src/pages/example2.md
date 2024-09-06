---
title: Mixing scripting languages
layout: "@layouts/MarkdownPage.astro"
---

<div class="blg-summary example">
<h3>Mixing scripting languages</h3>

<p class="text-muted">
    You are not limited to Bash scripts with Nextflow -- you can use any scripting language that can be executed by the Linux platform.
</p>

```groovy
#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Range
params.range = 100

/*
 * A trivial Perl script that produces a list of number pairs
 */
process perlTask {

    input:
    val x

    output:
    stdout

    shell:
    '''
    #!/usr/bin/env perl
    use strict;
    use warnings;

    my $count;
    my $range = !{x};
    for ($count = 0; $count < 10; $count++) {
        print rand($range) . ', ' . rand($range) . "\n";
    }
    '''
}

/*
 * A Python script which parses the output of the previous script
 */
process pyTask {

    input:
    stdin

    output:
    stdout

    script:
    """
    #!/usr/bin/env python
    import sys

    x = 0
    y = 0
    lines = 0
    for line in sys.stdin:
        items = line.strip().split(",")
        x += float(items[0])
        y += float(items[1])
        lines += 1

    print("avg: %s - %s" % ( x/lines, y/lines ))
    """
}

workflow {

    // A Perl script that produces a list of number pairs
    perlTask(params.range)

    // A Python script which parses the output of the previous script
    pyTask(perlTask.out)

    // View pyTask output
    pyTask.out.view()
}
```

</div>

### Synopsis

This example shows a simple Nextflow pipeline consisting of two processes written in different languages. The `perlTask` process starts with a Perl _shebang_ declaration and executes a Perl script that produces pairs of numbers. Since Perl uses the `$` character for variables, the special `shell` block is used instead of the normal `script` block to distinguish the Perl variables from Nextflow variables. Similarly, the `pyTask` process starts with a Python _shebang_ declaration. It takes the output from the Perl script and executes a Python script that averages the number pairs. The output from the `pyTask` process is then printed to screen.

### Try it

To try this pipeline:

1. Follow the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html#install-nextflow) to install Nextflow.
2. Copy the script above and save as `mixed-languages.nf`.
3. Launch the pipeline:

    nextflow run mixed-languages.nf

**NOTE**: To run this example with versions of Nextflow older than 22.04.0, you must include the `-dsl2` flag with `nextflow run`.
