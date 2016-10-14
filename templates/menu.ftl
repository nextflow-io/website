<div class="navbar navbar-inverse navbar-fixed-top" role="navigation">
<#if (content.type="index")>
    <a href="https://github.com/nextflow-io/nextflow" target='_blank'><img
            style="position: absolute; top: 0; right: 0; border: 0;"
            src="https://s3.amazonaws.com/github/ribbons/forkme_right_darkblue_121621.png" alt="Fork me on GitHub"></a>
</#if>

    <div class="container">
        <div class="navbar-header">
            <a class="navbar-brand" href="/index.html" style='padding: 8px'><img src='/img/nextflow2014_no-bg.png' height='40px' style='border: 0;'></a>
        </div>
        <div class="navbar-collapse collapse">

            <ul class="nav navbar-nav navbar-left hidden-xs">
                <li class="show animated flipInX" >
                    <a href="/index.html#Features" class="pull-right<#if (content.type="index")> scroll</#if>">Features</a>
                </li>
                <li class="show animated flipInX" >
                    <a href="/index.html#GetStarted" class="pull-right<#if (content.type="index")> scroll</#if>">Quick start</a>
                </li>

                <li class="dropdown show animated flipInX">
                    <a href="#" class="dropdown-toggle" data-toggle="dropdown">Examples <b class="caret"></b></a>
                    <ul class="dropdown-menu">
                        <li><a href="/example1.html">Basic pipeline</a></li>
                        <li><a href="/example2.html">Mixing scripting languages</a></li>
                        <li><a href="/example3.html">BLAST pipeline</a></li>
                        <li><a href="/example4.html">RNA-Seq pipeline</a></li>
                        <li><a href="https://github.com/nextflow-io/examples" target="_blank">More examples..</a></li>
                    </ul>
                </li>

                <li class="show animated flipInX">
                    <a href="/docs/latest/index.html" class="pull-right" target="_blank">Documentation</a>
                </li>

                <li class="show animated flipInX">
                    <a href="/blog.html">Blog</a>
                </li>

                <li class="show animated flipInX">
                    <a href="/about-us.html">About Us</a>
                </li>
                
                <li class="show">
                	<a href='http://www.crg.eu' target='_blank' style='padding: 10px 20px 0px 30px'><img src='/img/crg_logo.png' height='45'/></a>
                </li>
            </ul>
        </div>
    </div>
</div>
