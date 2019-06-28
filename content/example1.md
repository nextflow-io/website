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
sequences = file(params.in)

/* 
 * split a fasta file in multiple files 
 */
process splitSequences {

    input:
    file 'input.fa' from sequences

    output:
    file 'seq_*' into records

    """
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    """

}

/* 
 * Simple reverse the sequences 
 */
process reverse {

    input:
    file x from records
    
    output:
    stdout result

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

* Line 4: Defines a variable sequences holding a reference for the file whose name is
 specified by the <code>params.in</code> parameter.
 
* Lines 9-21: The process that splits the provided file.
 
* Line 11: Opens the input declaration block. The lines following this clause are interpreted 
as input definitions.

* Line 12: Defines the process input file. This file is received from the variable 
sequences and will be named <code>input.fa</code>.

* Line 14: Opens the output declaration block. Lines following this clause are interpreted 
as output definitions.

* Line 15: Defines that the process outputs files whose names match the pattern <code>seq_*</code>.
These files are sent over the channel records.

* Lines 17-19: The actual script executed by the process to split the provided file.
  
* Lines 26-37: Defines the second process, that receives the splits produced by the
previous process and reverses their content.

* Line 28: Opens the input declaration block. Lines following this clause are
interpreted as input definitions.

* Line 29: Defines the process input file. This file is received through the channel records.

* Line 31: Opens the output declaration block. Lines following this clause are
interpreted as output definitions.

* Line 32: The standard output of the executed script is declared as the process
output. This output is sent over the channel result.

* Lines 34-36: The actual script executed by the process to reverse the content of the
received files.

* Line 42: Prints a result each time a new item is received on the result channel.

### Try it in your computer 

Then launch the pipeline execution using this command: 

    $ ./nextflow run example-1.nf
