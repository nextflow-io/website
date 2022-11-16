<#include "header.ftl">
<#include "menu.ftl">


  <!-- Body -->
  <div class="wrapper"> <!-- wrapper -->
    <div class="container">
      <div class="row">
      <div class="col-xs-12">
        <h4 class='text-muted'>Talking about news in the Nextflow ecosytem and speaking with pioneers in the field.</h4>
        <hr>
      <div class="row">
           <div class="col-sm-8">
           <div class="timeline">
           <#assign count = 0>
           <#list podcasts as post>
             <div class="blg-summary">
               <h3 ><a href="${post.uri}"><#escape x as x?xml>${post.title}</#escape></a></h3>
               <h4 class="text-muted blg-description"><#escape x as x?xml>${post.description}</#escape></h4>
               <div class="timeline-info hidden-xs">
				        <img src="/img/${post.icon}" class="blg-author" alt="${post.author}">
				       </div>
               <ul class="text-muted list-inline blg-header">
               <li><i class="fa fa-info-circle"></i> ${post.subtype}</li>
               <li><i class="fa fa-calendar"></i> ${post.date?string("dd MMMM yyyy")} </li>
               <!--<li><i class="fa fa-comments-o"></i> 21 comments</li> -->
               </ul>
               <hr>
               <#if post.image??>
                <a href="${post.uri}"><img src="${post.image}" class="podcast-thumb-img"></a>
               </#if>
               <p class="blg-text">
                 <#assign MAX = 150>

                 <#assign words = post.body?word_list>
                 <#if words?size gt 150 >
                   <#assign body = words[0..MAX-1]?join(' ') >
                   ${body} .. (<a href="${post.uri}">click here to read more</a>)
                 <#else>
                   <#assign body = words?join(' ') >
                   ${body}
                 </#if>

               </p>

             </div>
          </#list>
           </div>

           <div class="clearfix"></div>
           </div>
           <div class="col-sm-4">
           <h3>Stay Tuned <small>Where to listen</small></h3>
           <hr>
           <ul class="blg-social">
             <li>
               <a href="https://open.spotify.com/show/1slEz7EL46cHa9vdRmPLY4" target='_blank'>
                 <span class="icon spotify">
                   <img src="/img/spotify.png" alt="Spotify">
                 </span>
                 <span class="text-muted">Listen on Spotify</span>
               </a>
             </li>
             <li>
               <a href="https://podcasts.apple.com/ca/podcast/the-next-bio-informatics-podcast/id1554921146" target='_blank'>
                <span class="icon applepodcasts">
                   <img src="/img/apple_podcasts.png" alt="Apple Podcasts">
                 </span>
                 <span class="text-muted">Listen on Apple Podcasts</span>
               </a>
             </li>
             <li>
               <a href="https://www.google.com/podcasts?feed=aHR0cHM6Ly9hbmNob3IuZm0vcy80ZTAzZGMxOC9wb2RjYXN0L3Jzcw==" target='_blank'>
                <span class="icon googlepodcasts">
                   <img src="/img/google_podcasts.png" alt="Google Podcasts">
                 </span>
                 <span class="text-muted">Listen on Google Podcasts</span>
               </a>
             </li>
             <li>
               <a href="https://www.youtube.com/playlist?list=PLPZ8WHdZGxmUAV23hZ9lcZFtt3MAa9IOj" target='_blank'>
                 <span class="icon googlepodcasts">
                   <img src="/img/youtube.png" alt="YoyTube">
                </span>
                <span class="text-muted">Watch on YouTube</span>
             </a>
             </li>
             <li>
               <a href="https://twitter.com/nextflowio" target='_blank'>
                 <span class="icon twitter">
                 <i class="fa fa-twitter fa-2x"></i>
               </span>
               <span class="text-muted">Watch on Twitter spaces</span>
             </a>
             </li>
            <li>
              <a href="https://anchor.fm/s/4e03dc18/podcast/rss">
                <span class="icon rss">
                  <i class="fa fa-rss fa-2x"></i>
                </span>
                <span class="text-muted">Podcast RSS Feed</span>
              </a>
            </li>
           </ul>
           </div>
         </div>
        <hr>
      </div>
    </div>

    </div>
  </div> <!-- / wrapper -->



<#include "footer.ftl">
