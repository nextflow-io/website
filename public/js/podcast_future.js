var podcasts = {};
$('.podcast').each(function () {
  var golive = parseFloat($(this).data("golive"));
	var now = new Date().getTime();
	if (golive > now) {
		$(this).addClass('podcast-future')
    podcasts[golive] = $(this);
  }
});

// Sort the podcasts by date and hide all by the first one
var sortedPodcasts = Object.keys(podcasts).sort();
// NB: Starts at i=1
for (var i = 1; i < sortedPodcasts.length; i++) {
  podcasts[sortedPodcasts[i]][0].style.display = "none";
}
