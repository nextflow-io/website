<#include "header.ftl">
<#include "menu.ftl">


	<!-- Body -->
	<div class="wrapper"> <!-- wrapper -->
	  <div class="container">
	    <div class="row">
		  <div class="col-xs-12">
		    <h4 class='text-muted'>Blogging about Nextflow, computational workflows, containers and cloud computing</h4>
		    <hr>
		  </div>
		</div>
	    <div class="row">
		  <div class="col-sm-8">
			<div class="timeline">
			<#assign count = 0>
			<#list published_posts as post>
			  <div class="blg-summary">
			    <h3 ><a href="${post.uri}"><#escape x as x?xml>${post.title}</#escape></a></h3>
				<div class="timeline-info hidden-xs">
				  <img src="/img/${post.icon}" class="blg-author" alt="${post.author}">
				</div>
			    <ul class="text-muted list-inline blg-header">
				  <li><i class="fa fa-user"></i> ${post.author}</li>
				  <li><i class="fa fa-calendar"></i> ${post.date?string("dd MMMM yyyy")} </li>
				  <!--<li><i class="fa fa-comments-o"></i> 21 comments</li> -->
			    </ul>
			    <hr>
			    <p class="blg-text">
						<#if post.title == "Nextflow workshop is coming!" >
							<#assign MAX = 50>
						<#elseif post.title == "More syntax sugar for Nextflow developers!" >
							<#assign MAX = 70>
						<#elseif post.title == "Easy provenance reporting" >
							<#assign MAX = 70>
						<#elseif post.title == "Bringing Nextflow to Google Cloud Platform with WuXi NextCODE" >
							<#assign MAX = 108>
						<#elseif post.title == "Nextflow Summit 2022" >
							<#assign MAX = 140>
						<#elseif post.title == "Nextflow Summit 2022 Recap" >
							<#assign MAX = 140>
						<#else>
							<#assign MAX = 150>
						</#if>

						<#assign post_intro = post.body?split("<!-- end-archive-description -->")[0]>
						<#assign words = post_intro?word_list>
						<#if words?size gt 150 >
							<#assign body = words[0..MAX-1]?join(' ') >
							${body} .. (<a href="${post.uri}">click here to read more</a>)
						<#else>
							<#assign body = words?join(' ') >
							${body}
						</#if>

			    </p>

			    <!--
			    <p class="tags">
				  <a href="#" class="background-color bg-hover-color">Bootstrap</a>
				  <a href="#" class="background-color bg-hover-color">HTML5</a>
				  <a href="#" class="background-color bg-hover-color">CSS</a>
				  <a href="#" class="background-color bg-hover-color">jQuery</a>
			    </p>
			    -->
			  </div>
            </#list>
			</div>

		    <!-- Pagination
		    <ul class="pagination pull-right">
			  <li class="disabled"><a href="#">&laquo;</a></li>
			  <li class="active"><a href="#">1 <span class="sr-only">(current)</span></a></li>
			  <li><a href="#">2</a></li>
			  <li><a href="#">3</a></li>
			  <li><a href="#">4</a></li>
			  <li><a href="#">5</a></li>
			  <li><a href="#">&raquo;</a></li>
		    </ul>
		    -->

			<div class="clearfix"></div>
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
				  <span class="text-muted">Subscribe to our blog RSS Feed.</span>
				</a>
			  </li>
			  <li>
			    <a href="https://twitter.com/nextflowio" target='_blank'>
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
			<h3>Tags</h3>
			<hr>
			<p class="tags">
			  <a href="#" class="background-color bg-hover-color">Bootstrap</a>
			  <a href="#" class="background-color bg-hover-color">HTML5</a>
			  <a href="#" class="background-color bg-hover-color">CSS</a>
			  <a href="#" class="background-color bg-hover-color">jQuery</a>
			  <a href="#" class="background-color bg-hover-color">Lorem</a>
			  <a href="#" class="background-color bg-hover-color">ipsum</a>
			  <a href="#" class="background-color bg-hover-color">dolor</a>
			  <a href="#" class="background-color bg-hover-color">sit</a>
			  <a href="#" class="background-color bg-hover-color">amet</a>
			  <a href="#" class="background-color bg-hover-color">consectetur</a>
			  <a href="#" class="background-color bg-hover-color">adipiscing</a>
			  <a href="#" class="background-color bg-hover-color">elit</a>
			  <a href="#" class="background-color bg-hover-color">you</a>
			  <a href="#" class="background-color bg-hover-color">should</a>
			  <a href="#" class="background-color bg-hover-color">like</a>
			  <a href="#" class="background-color bg-hover-color">this</a>
			</p>
    -->
		  </div>
		</div>
	  </div>
	</div> <!-- / wrapper -->



<#include "footer.ftl">
