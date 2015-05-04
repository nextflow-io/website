title=Introducing Nextflow REPL Console
date=2015-04-14
type=post
tags=data-analysis,pipelines,repl,groovy
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

The latest version of Nextflow introduces a new *console* graphical interface. 

The Nextflow console is a REPL ([read-eval-print loop](http://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop)) 
environment that allows one to quickly test part of a script or pieces of Nextflow code 
in an interactive manner. 

It is a handy tool that allows one to evaluate fragments of Nextflow/Groovy code 
or fast prototype a complete pipeline script. 


### Getting started 

The console application is included in the latest version of Nextflow 
([0.13.1](https://github.com/nextflow-io/nextflow/releases) or higher). 

You can try this feature out, having Nextflow installed on your computer, by entering the 
following command in your shell terminal: ``nextflow console ``.
    
When you execute it for the first time, Nextflow will spend a few seconds downloading 
the required runtime dependencies. When complete the console window will appear as shown in 
the picture below. 

<img src='/img/nextflow-console1.png' alt="Nextflow console" style='width: 100%; padding: 2em 1em;' />

It contains a text editor (the top white box) that allows you to enter and modify code snippets. 
The results area (the bottom yellow box) will show the executed code's output. 

At the top you will find the menu bar (not shown in this picture) and the actions 
toolbar that allows you to open, save, execute (etc.) the code been tested.

As a practical execution example, simply copy and paste the following piece of code in the 
console editor box:

    echo true 

    process sayHello {
  
     """
     echo Hello world
     """ 
 
    }
    
    
Then, in order to evaluate it, open the ``Script`` menu in the top menu bar and select the ``Run`` 
command. Alternatively you can use the ``CTRL+R`` keyboard shortcut to run it (``⌘+R`` on the Mac). 
In the result box an output similar to the following will appear: 

    [warm up] executor > local
    [00/d78a0f] Submitted process > sayHello (1)
    Hello world
 
Now you can try to modify the entered process script, execute it again and check that 
the printed result has changed. 

If the output doesn't appear, open the ``View`` menu and make sure that the entry ``Capture Standard
Output`` is selected (it must have a tick on the left).

It is worth noting that the global script context is maintained across script executions. 
This means that variables declared in the global script scope are not lost when the 
script run is complete, and they can be accessed in further executions of the same or another 
piece of code.

In order to reset the global context you can use the command ``Clear Script Context``
available in the ``Script`` menu.  


### Conclusion

The Nextflow console is a REPL environment which allows you to experiment and get used 
to the Nextflow programming environment. By using it you can prototype or test your code 
without the need to create/edit script files.

Note: the Nextflow console is implemented by sub-classing the [Groovy console](http://groovy-lang.org/groovyconsole.html) tool. 
For this reason you may find some labels that refer to the Groovy programming environment
in this program. 

 
