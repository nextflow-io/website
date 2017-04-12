<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9"> 
<#list published_content as content>
<#if (content.uri = 'index.html')>
<url><loc>${config.site_host}</loc><lastmod>${content.date?string("yyyy-MM-dd")}</lastmod><priority>1.0</priority></url>
<url><loc>${config.site_host}/${content.uri}</loc><lastmod>${content.date?string("yyyy-MM-dd")}</lastmod><priority>0.85</priority></url>
<#elseif (content.uri != 'error-page.html')>
<url><loc>${config.site_host}/${content.uri}</loc><lastmod>${content.date?string("yyyy-MM-dd")}</lastmod><priority>0.85</priority></url>
</#if>
</#list>
<url><loc>${config.site_host}/blog.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><priority>0.85</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/index.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.85</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/getstarted.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/basic.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/script.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/process.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/channel.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/operator.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/executor.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/config.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/amazons3.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/docker.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.85</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/dnanexus.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/ignite.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/tracing.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.85</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/sharing.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.85</priority></url>
<url><loc>https://www.nextflow.io/docs/latest/metadata.html</loc><lastmod>${published_date?string("yyyy-MM-dd")}</lastmod><changefreq>weekly</changefreq><priority>0.69</priority></url>
</urlset>