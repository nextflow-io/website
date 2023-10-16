<!DOCTYPE html>
<html lang="en">
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    <meta charset="utf-8">
    <!--<meta http-equiv="X-UA-Compatible" content="IE=edge">-->
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="author" content="">
    <link rel="shortcut icon" href="/img/favicon.png">
    <link rel="alternate" type="application/rss+xml" title="RSS" href="https://www.nextflow.io/feed.xml">
    <title><#if (content.title)??><#if (content.uri!='error-page.html')><#escape x as x?xml>${content.title} | Nextflow</#escape><#else>${content.title}</#if><#else>Nextflow</#if></title>
    <#if (content.description)??>
    <meta name="description" content="<#escape x as x?xml>${content.description}</#escape>">
    </#if>
    <#if (content.title)??>
    <meta property="og:title" content="<#escape x as x?xml>${content.title} | Nextflow</#escape>" />
    </#if>
    <#if (content.description)??>
    <meta property="og:description" content="<#escape x as x?xml>${content.description}</#escape>" />
    </#if>
    <#if (content.image)??>
    <meta property="og:image" content="https://nextflow.io/${content.image}" />
    <#else>
    <meta property="og:image" content="https://nextflow.io/img/share.png" />
    </#if>
    <meta property="og:type" content="website" />
    <meta data-react-helmet="true" name="twitter:card" content="summary_large_image" />
    <meta data-react-helmet="true" name="twitter:creator" content="@nextflowio" />
    <#if (content.title)??>
    <meta property="twitter:title" content="<#escape x as x?xml>${content.title} | Nextflow</#escape>" />
    </#if>
    <#if (content.description)??>
    <meta property="twitter:description" content="<#escape x as x?xml>${content.description}</#escape>" />
    </#if>
    <#if (content.image)??>
    <meta property="twitter:image" content="https://nextflow.io/${content.image}" />
    <#else>
    <meta property="twitter:image" content="https://nextflow.io/img/share.png" />
    </#if>

    <!-- Bootstrap core CSS -->
    <link href="/css/bootstrap.css" rel="stylesheet">
    <!-- HighlightJS -->
    <link href="/css/xcode.min.css" rel="stylesheet">

    <!-- Custom styles for this template -->
    <link href="/css/color-styles.css" rel="stylesheet">
    <link href="/css/ui-elements.css" rel="stylesheet">
    <link href="/css/custom.css?${.now?long?c}" rel="stylesheet">

    <!-- Resources -->
    <link href="/css/animate.css" rel="stylesheet">
    <link href="https://netdna.bootstrapcdn.com/font-awesome/4.0.3/css/font-awesome.css" rel="stylesheet">
    <link href='https://fonts.googleapis.com/css?family=Lobster' rel='stylesheet' type='text/css'>
		<link href="/fonts/degular/degular.css" rel="stylesheet">

    <!-- HTML5 shim and Respond.js IE8 support of HTML5 elements and media queries -->
    <!--[if lt IE 9]>
    <script src="https://oss.maxcdn.com/libs/html5shiv/3.7.0/html5shiv.js"></script>
    <script src="https://oss.maxcdn.com/libs/respond.js/1.3.0/respond.min.js"></script>
    <![endif]-->

    <!-- Google Tag Manager -->
    <script>(function(w,d,s,l,i){w[l]=w[l]||[];w[l].push({'gtm.start':
    new Date().getTime(),event:'gtm.js'});var f=d.getElementsByTagName(s)[0],
    j=d.createElement(s),dl=l!='dataLayer'?'&l='+l:'';j.async=true;j.src=
    'https://www.googletagmanager.com/gtm.js?id='+i+dl;f.parentNode.insertBefore(j,f);
    })(window,document,'script','dataLayer','GTM-TNCXSWG');</script>
    <!-- End Google Tag Manager -->
</head>

<body class="body-green">
    <!-- Google Tag Manager (noscript) -->
    <noscript><iframe src="https://www.googletagmanager.com/ns.html?id=GTM-TNCXSWG"
    height="0" width="0" style="display:none;visibility:hidden"></iframe></noscript>
    <!-- End Google Tag Manager (noscript) -->
