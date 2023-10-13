<div class="navbar navbar-inverse navbar-fixed-top" role="navigation">
<#if (content.type="index")>
    <a href="https://github.com/nextflow-io/nextflow" target='_blank' class="hidden-xs"><img
            style="position: absolute; top: 0; right: 0; border: 0;"
            src="/img/forkme_right_darkblue_121621.png" alt="Fork me on GitHub"></a>
</#if>

    <div class="container">
        <div class="navbar-header">
            <button type="button" class="navbar-toggle" data-toggle="collapse" data-target=".navbar-collapse">
                <span class="sr-only">Toggle navigation</span>
                <span class="icon-bar"></span>
                <span class="icon-bar"></span>
                <span class="icon-bar"></span>
            </button>
            <a class="navbar-brand" href="/index.html" style='padding: 15px'><img src='/img/nextflow.svg' height='35px' style='border: 0;'></a>
        </div>
        <div class="navbar-collapse collapse">

            <ul class="nav navbar-nav">
                <li class="show animated flipInX" >
                    <a href="/index.html#Features" class="<#if (content.type="index")> scroll</#if>">Features</a>
                </li>
                <li class="show animated flipInX" >
                    <a href="/index.html#GetStarted" class="<#if (content.type="index")> scroll</#if>">Quick start</a>
                </li>

                <li class="show animated flipInX">
                    <a href="/blog/2023/learn-nextflow-in-2023.html">Learn Nextflow</a>
                </li>

                <li class="dropdown show animated flipInX">
                    <a href="#" class="dropdown-toggle" data-toggle="dropdown">Examples <b class="caret"></b></a>
                    <ul class="dropdown-menu">
                        <li><a href="/example1.html">Basic pipeline</a></li>
                        <li><a href="/example2.html">Mixing scripting languages</a></li>
                        <li><a href="/example3.html">BLAST pipeline</a></li>
                        <li><a href="/example4.html">RNA-Seq pipeline</a></li>
                        <li><a href="/example5.html">Machine Learning pipeline</a></li>
                        <li><a href="https://github.com/CRG-CNAG/CalliNGS-NF/" target="_blank">Variant Calling pipeline</a></li>
                    </ul>
                </li>

                <li class="dropdown show animated flipInX">
                    <a href="#" class="dropdown-toggle" data-toggle="dropdown">Community<b class="caret"></b></a>
                    <ul class="dropdown-menu">
                        <li><a href="/docs/latest/index.html">Reference documentation</a></li>
                        <li><a href="/docs/edge/index.html">Edge release documentation</a></li>
                        <li><a href="http://training.nextflow.io">Nextflow training</a></li>
                        <li><a href="http://nextflow-io.github.io/patterns/index.html">Implementation patterns</a></li>
                        <li><a href="https://github.com/nextflow-io/nextflow">GitHub repository</a></li>
                        <li><a href="https://community.seqera.io/c/nextflow/5">Community forum</a></li>
                        <li><a href="https://github.com/nextflow-io/nextflow/discussions">GitHub discussions</a></li>
                        <li><a href="https://www.nextflow.io/slack-invite.html">Slack community chat</a></li>
                        <li><a href="https://nf-co.re">nf-core Community pipelines</a></li>
                    </ul>
                </li>

                <li class="show animated flipInX">
                    <a href="/blog.html">Blog</a>
                </li>

                <li class="show animated flipInX">
                    <a href="/podcasts.html">Podcast</a>
                </li>

                <li class="show animated flipInX">
                    <a href="/about-us.html">About</a>
                </li>

            </ul>

        </div>
    </div>
</div>
