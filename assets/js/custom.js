/*------------------------------------------------------------------
Project:	Heaven
Version:	1.1
Created: 		27/11/2013
Last change:	12/01/2014
-------------------------------------------------------------------*/

// Recent works thumbnail image height resize
//===========================================

$('.recent-works .thumbnail > .image').on( 'resize', function () {
    $('.recent-works .thumbnail > .image').height( $('.recent-works .thumbnail > .image').width() / 1.6 );
}).resize();

// Sign In & Sign Out
// ==================

$('#sign-in').on('click', function() {
	$("#user-bar").toggleClass("show hidden");
	$("#user-bar").toggleClass("animated flipInX");
	$("#sign-in").toggleClass("hidden show");
	$("#sign-up").toggleClass("hidden show");
	$("#sign-in").removeClass("animated flipInX");
	$("#sign-up").removeClass("animated flipInX");
	return false;
});

$('#sign-out').on('click', function() {
	$("#user-bar").toggleClass("show hidden");
	$("#user-bar").toggleClass("animated flipInX");
	$("#sign-in").toggleClass("hidden show");
	$("#sign-in").addClass("animated flipInX");
	$("#sign-up").toggleClass("hidden show");
	$("#sign-up").addClass("animated flipInX");
	return false;
});

// Style Toggle
// ============

$('.style-toggle-btn').on('click', function() {
	$(".style-toggle-btn").toggleClass("show hidden");
	$(".style-toggle").toggleClass("hidden show");
	return false;
});

$('.style-toggle-close').on('click', function() {
	$(".style-toggle-btn").toggleClass("show hidden");
	$(".style-toggle").toggleClass("hidden show");
	return false;
});

// Body Style Change
// =================

$('.style-toggle ul > li.green').on('click', function() {
	$("body").addClass("body-green");
	$("body").removeClass("body-blue");
	$("body").removeClass("body-red");
	$("body").removeClass("body-amethyst");
	return false;
});

$('.style-toggle ul > li.blue').on('click', function() {
	$("body").addClass("body-blue");
	$("body").removeClass("body-green");
	$("body").removeClass("body-red");
	$("body").removeClass("body-amethyst");
	return false;
});

$('.style-toggle ul > li.red').on('click', function() {
	$("body").addClass("body-red");
	$("body").removeClass("body-green");
	$("body").removeClass("body-blue");
	$("body").removeClass("body-amethyst");
	return false;
});

$('.style-toggle ul > li.amethyst').on('click', function() {
	$("body").addClass("body-amethyst");
	$("body").removeClass("body-blue");
	$("body").removeClass("body-red");
	$("body").removeClass("body-green");
	return false;
});

// Lost password
// =============

$('#lost-btn').on('click', function() {
	$("#lost-form").toggleClass("show hidden");
	return false;
});

// Contact Us
// ==========

$('#signed-in').on('click', function() {
	$(".form-white > .contact-avatar > span").toggleClass("show hidden");
	$(".form-white > .contact-avatar > img").toggleClass("show hidden animated flipInX");
	return false;
});

$('#signed-in').on('click', function() {
	$("#email-contact").toggleClass("show hidden");
	$("#email-contact-disabled").toggleClass("show hidden");
	$("#name-1").toggleClass("show hidden");
	$("#name-1-disabled").toggleClass("show hidden");
	$("#name-2").toggleClass("show hidden");
	$("#name-2-disabled").toggleClass("show hidden");
	return false;
});

$('#signed-in').on('click', function() {
	$("#signed-in > i").toggleClass("fa-circle-o fa-dot-circle-o");
	return false;
});

// Search box toggle
// =================
$('#search-btn').on('click', function() {
	$("#search-icon").toggleClass('fa-search fa-times margin-2');
	$("#search-box").toggleClass('hidden show animated flipInX');
	return false;
});

// Error page
// ==========

var divs = $("i.random").get().sort(function(){
  return Math.round(Math.random())-0.5; //random so we get the right +/- combo
}).slice(0,1)
$(divs).show();

var divs = $("i.random2").get().sort(function(){
  return Math.round(Math.random())-0.5; //random so we get the right +/- combo
}).slice(0,1)
$(divs).show();

var divs = $("i.random3").get().sort(function(){
  return Math.round(Math.random())-0.5; //random so we get the right +/- combo
}).slice(0,1)
$(divs).show();

// Corporate Index Features
// ========================

$('.crp-ft').hover (function() {
	$(this).children("a").toggleClass("show hidden");
	$(this).children("a").toggleClass("animated flipInX");
	return false;
});

// Smooth scroll
// =================

$(".scroll").click(function (event) {
    event.preventDefault();
    //calculate destination place
    var dest = 0;
    if ($(this.hash).offset().top > $(document).height() - $(window).height()) {
        dest = $(document).height() - $(window).height();
    } else {
        dest = $(this.hash).offset().top;
    }
    //go to destination
    $('html,body').animate({
        scrollTop: dest
    }, 1000, 'swing');
});

// Podcast - styling for future events
// ====================================
$('.podcast').each(function () {
    var golive = parseFloat($(this).data("golive"));
	var now = new Date().getTime();
	if (golive > now) {
		$(this).addClass('podcast-future')
    }
});
