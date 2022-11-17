<#include "header.ftl">
<#include "menu.ftl">

    <!-- Body -->
  <div class="wrapper"> <!-- wrapper -->
    <div class="container">
      <div class="row">
      <div class="col-sm-8">
      <div class="blg-summary">
        <h2><a href="${post.uri}"><span class="label label-success">Episode ${post.episode}</span> ${post.subtype?cap_first}: <#escape x as x?xml>${post.title}</#escape></a></h2>
        <h4 class="text-muted blg-description"><#escape x as x?xml>${content.description}</#escape></h4>
        <ul class="text-muted list-inline blg-header">
        <li><i class="fa fa-info-circle"></i> ${content.subtype}</li>
        <li><i class="fa fa-calendar"></i> ${content.date?string("dd MMMM yyyy")}</li>
        <!--<li><i class="fa fa-comments-o"></i> 21 comments</li> -->
        </ul>
        <hr>
        <div class="blg-text">
          <#if content.youtubeid?length gt 6 >
            <iframe width="560" height="315" src="https://www.youtube.com/embed/${content.youtubeid}" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
            <hr>
          </#if>
          ${content.body}
        </div>

        <p class="tags">
                <#list content.tags as x >
        <a href="/tags/${x}.html" class="background-color bg-hover-color">${x}</a>
        </#list>
        </p>

      </div>
      <!-- COMMENTS -->
      <div class="comments">
        <!--<h4>3 comments</h4>-->

          <div id="disqus_thread"></div>
                <script type="text/javascript">
                    if( window.location.hostname == 'localhost' || window.location.hostname == '127.0.0.1' ) { throw new Error('Skip Disqus on localhost') };
                    /* * * CONFIGURATION VARIABLES: EDIT BEFORE PASTING INTO YOUR WEBPAGE * * */
                    var disqus_shortname = 'nextflow'; // required: replace example with your forum shortname

                    /* * * DON'T EDIT BELOW THIS LINE * * */
                    (function() {
                        var dsq = document.createElement('script'); dsq.type = 'text/javascript'; dsq.async = true;
                        dsq.src = '//' + disqus_shortname + '.disqus.com/embed.js';
                        (document.getElementsByTagName('head')[0] || document.getElementsByTagName('body')[0]).appendChild(dsq);
                    })();
                </script>
                <noscript>Please enable JavaScript to view the <a href="http://disqus.com/?ref_noscript">comments powered by Disqus.</a></noscript>
                <a href="http://disqus.com" class="dsq-brlink">comments powered by <span class="logo-disqus">Disqus</span></a>

      </div>

      <!-- // Comments -->
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
            <!--
      <h3>Popular Stories</h3>
      <hr>
      <div class="block-inverse">
        <ul>
        <li><a href="#">Lorem ipsum dolor sit amet, consectetur adipiscing elit.</a></li>
        <li><a href="#">Nullam id ipsum varius, tincidunt odio nec, placerat.</a></li>
        <li><a href="#">Sed sit amet auctor augue, nec dignissim ligula.</a></li>
        <li><a href="#">Lorem ipsum dolor sit amet, consectetur adipiscing elit.</a></li>
        <li><a href="#">Lorem ipsum dolor sit amet, consectetur adipiscing elit.</a></li>
        </ul>
      </div>
            -->
      </div>
    </div>
    </div>
  </div> <!-- / wrapper -->


<#include "footer.ftl">
