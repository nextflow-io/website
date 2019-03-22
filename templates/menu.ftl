<div class="navbar navbar-inverse navbar-fixed-top" role="navigation">
<#if (content.type="index")>
    <a href="https://github.com/nextflow-io/nextflow" target='_blank' class="hidden-xs"><img
            style="position: absolute; top: 0; right: 0; border: 0;"
            src="https://s3.amazonaws.com/github/ribbons/forkme_right_darkblue_121621.png" alt="Fork me on GitHub"></a>
</#if>

    <div class="container">
        <div class="navbar-header">
            <button type="button" class="navbar-toggle" data-toggle="collapse" data-target=".navbar-collapse">
                <span class="sr-only">Toggle navigation</span>
                <span class="icon-bar"></span>
                <span class="icon-bar"></span>
                <span class="icon-bar"></span>
            </button>
            <a class="navbar-brand" href="/index.html" style='padding: 8px'><img src='/img/nextflow2014_no-bg.png' height='40px' style='border: 0;'></a>
        </div>
        <div class="navbar-collapse collapse">

            <ul class="nav navbar-nav">
                <li class="show animated flipInX" >
                    <a href="/index.html#Features" class="<#if (content.type="index")> scroll</#if>">Features</a>
                </li>
                <li class="show animated flipInX" >
                    <a href="/index.html#GetStarted" class="<#if (content.type="index")> scroll</#if>">Quick start</a>
                </li>

                <li class="dropdown show animated flipInX">
                    <a href="#" class="dropdown-toggle" data-toggle="dropdown">Examples <b class="caret"></b></a>
                    <ul class="dropdown-menu">
                        <li><a href="/example1.html">Basic pipeline</a></li>
                        <li><a href="/example2.html">Mixing scripting languages</a></li>
                        <li><a href="/example3.html">BLAST pipeline</a></li>
                        <li><a href="/example4.html">RNA-Seq pipeline</a></li>
                        <li><a href="https://github.com/CRG-CNAG/CalliNGS-NF/" target="_blank">Variant Calling pipeline</a></li>
                    </ul>
                </li>

                <li class="dropdown show animated flipInX">
                    <a href="#" class="dropdown-toggle" data-toggle="dropdown">Developers<b class="caret"></b></a>
                    <ul class="dropdown-menu">
                        <li><a href="/docs/latest/index.html">Reference documentation</a></li>
                        <li><a href="/docs/edge/index.html">Edge release documentation</a></li>
                        <li><a href="http://nextflow-io.github.io/patterns/index.html">Implementation patterns</a></li>
                        <li><a href="https://github.com/nextflow-io/nextflow">GitHub repository</a></li>
                        <li><a href="https://groups.google.com/forum/#!forum/nextflow">Discussion forum</a></li>
                        <li><a href="https://gitter.im/nextflow-io/nextflow">Gitter chat</a></li>
                    </ul>
                </li>

                <li class="show animated flipInX">
                    <a href="/blog.html">Blog</a>
                </li>

                <li class="show animated flipInX">
                    <a href="/about-us.html">About Us</a>
                </li>
                   
            </ul>
            <ul class="list-inline pull-right hidden-xs">
                <li style='position: relative; top: 10px; right: 75px' >
                	<a href='http://www.crg.eu/' target='_blank' style='padding: 0'><img src='/img/crg_logo.png' height='45'/></a>
                </li>
            </ul>
        </div>
    </div>
</div>
