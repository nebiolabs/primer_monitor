import jQuery from "jquery";
//if jquery has not been loaded (window.jQuery is not null or undefined)
if(window.jQuery == null) {
    window.jQuery = jQuery
    window.$ = jQuery
}
