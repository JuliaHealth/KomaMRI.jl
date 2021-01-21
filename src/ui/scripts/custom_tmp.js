window.$ = window.jQuery = require('/home/ccp/Documents/MRIsim.jl/src/ui/scripts/jquery-3.4.1.slim.min.js');

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
