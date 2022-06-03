title=Basic pipeline
date=2014-07-23
type=page
status=published
syntaxhighlighter=yes
~~~~~~

<div class="blg-summary example">
<h3><a href="javascript:void(0)">Basic pipeline</a></h3>

<p class="text-muted" >
	This example shows how write a pipeline made up of two simple BASH processes, so
	that the results produced by the first are consumed by the second process.
</p>


<script type="syntaxhighlighter" class="brush: groovy">
<![CDATA[
#!/usr/bin/env nextflow

params.in = "$baseDir/data/sample.fa"

/* 
 * split a fasta file in multiple files 
 */
process splitSequences {

    input:
    path 'input.fa' from params.in

    output:
    path 'seq_*' into records

    """
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    """
}

/* 
 * Simple reverse the sequences 
 */
process reverse {

    input:
    path x from records
    
    output:
    stdout into result

    """
    cat $x | rev
    """
}

/* 
 * print the channel content 
 */
result.subscribe { println it }


]]>
</script>
</div>


### Synopsis

* Line 1: The script starts with a shebang declaration. This allows you to launch your pipeline, 
 as any other BASH script.

* Line 3: Declares a pipeline parameter named <code>params.in</code> that is initialized with the value
 <code>$HOME/sample.fa</code>.This value can be overridden when launching the pipeline, 
 by simply adding the option <code>--in &lt;value&gt;</code> to the script command line.
 
* Lines 8-20: The process that splits the provided file.
 
* Line 10: Opens the input declaration block. The lines following this clause are interpreted 
as input definitions.

* Line 11: Declares the process input file. This file is taken from the `params.in` parameter 
and named <code>input.fa</code>.

* Line 13: Opens the output declaration block. Lines following this clause are interpreted 
as output declarations.

* Line 14: Defines that the process outputs files whose names match the pattern <code>seq_*</code>.
These files are sent over the channel `records`.

* Lines 16-18: The actual script executed by the process to split the provided file.
  
* Lines 24-35: Defines the second process, that receives the splits produced by the
previous process and reverses their content.

* Line 26: Opens the input declaration block. Lines following this clause are
interpreted as input declarations.

* Line 27: Defines the process input file. This file is received through the channel `records`.

* Line 29: Opens the output declaration block. Lines following this clause are
interpreted as output declarations.

* Line 30: The standard output of the executed script is declared as the process
output. This output is sent over the channel `result`.

* Lines 32-34: The actual script executed by the process to reverse the content of the
received files.

* Line 40: Prints a result each time a new item is received on the `result` channel.

### Try it in your computer 

Then launch the pipeline execution using this command: 

    $ ./nextflow run example-1.nf
