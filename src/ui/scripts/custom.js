window.$ = window.jQuery = require('JQUERY');

$(document).ready(function() {
  // executes when HTML-Document is loaded and DOM is ready
  $('.navbar-nav>li>a').on('click', function(){
    $('.navbar-collapse').collapse('hide');
  });
  console.log("document is ready");
});
$(window).ready(function(){
  console.log("window is ready");
});
