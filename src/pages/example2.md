---
title: Mixing scripting languages
layout: "@layouts/Page.astro"
---

<div class="blg-summary example">
<h3>Mixing scripting languages</h3>

<p class="text-muted">
    With Nextflow, you are not limited to Bash scripts -- you can use any scripting language! In other words, for each <i>process</i> you can use the language that best fits the specific task or that you simply prefer.
</p>

```groovy
#!/usr/bin/env nextflow

params.range = 100

/*
 * A trivial Perl script that produces a list of number pairs
 */
process perlTask {
    output:
    stdout

    shell:
    '''
    #!/usr/bin/env perl
    use strict;
    use warnings;

    my $count;
    my $range = !{params.range};
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
    perlTask | pyTask | view
}
```

</div>

### Synopsis

In the above example we define a simple pipeline with two processes.

The first process executes a Perl script, because the script block definition starts
with a Perl _shebang_ declaration (line 14). Since Perl uses the `$` character for variables, we use the special `shell` block instead of the normal `script` block to easily distinguish the Perl variables from the Nextflow variables.

In the same way, the second process will execute a Python script, because the script block starts with a Python shebang (line 36).
