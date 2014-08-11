<#include "header.ftl">
<#include "menu.ftl">
	
    <!-- Body -->
	<div class="wrapper"> <!-- wrapper -->
	  <div class="container">
	    <div class="row">
            <div class="col-sm-8">
                <h3>Tag: ${tag}</h3>
                <hr>
                <ul>
                    <#list tag_posts as post>
                    <#if (last_month)??>
                        <#if post.date?string("MMMM yyyy") != last_month>
                            </ul>
                            <h4>${post.date?string("MMMM yyyy")}</h4>
                            <ul>
                        </#if>
                    <#else>
                        <h4>${post.date?string("MMMM yyyy")}</h4>
                        <ul>
                    </#if>
        
                    <li>${post.date?string("dd")} - <a href="/${post.uri}">${post.title}</a></li>
                    <#assign last_month = post.date?string("MMMM yyyy")>
                    </#list>
                </ul>              

            </div>
		  <div class="col-sm-4">
			<h3>Stay Tuned <small>Social Links</small></h3>
			<hr>
			<ul class="blg-social">
			  <li>
			    <a href="/feed.xml">
			      <span class="icon rss">
				    <i class="fa fa-rss fa-2x"></i>
				  </span>
				  <span class="text-muted">Subscribe to our RSS Feed.</span>
				</a>
			  </li>
			  <li>
			    <a href="https://twitter.com/nextflowio">
			      <span class="icon twitter">
				    <i class="fa fa-twitter fa-2x"></i>
				  </span>
				  <span class="text-muted">Subscribe to our Twitter Feed.</span>
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