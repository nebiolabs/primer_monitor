// datatables.net-responsive-bm@4.1.1 downloaded from https://ga.jspm.io/npm:datatables.net-responsive-bm@4.1.1/js/responsive.bulma.mjs

import e from"datatables.net-bm";export{default}from"datatables.net-bm";import"datatables.net-responsive";
/*! Responsive Bulma styling 4.1.1 for DataTables
* Copyright (c) SpryMedia Ltd - datatables.net/license
*/
var t=e.Dom,n=e.Responsive.display,r;function i(){return r||(r=t.c(`div`).classAdd(`modal DTED`).append(t.c(`div`).classAdd(`modal-background`)).append(t.c(`div`).classAdd(`modal-content`).append(t.c(`div`).classAdd(`modal-header`)).append(t.c(`div`).classAdd(`modal-body`))).append(t.c(`button`).attr(`type`,`button`).attr(`aria-label`,`Close`).classAdd(`modal-close is-large`))),r}n.modal=function(e){return function(n,r,a,o){var s=a(),c=i();if(s===!1)return!1;if(!r){if(e&&e.header){var l=c.find(`div.modal-header`);l.find(`button`).detach(),l.empty().append(t.c(`h4`).classAdd(`modal-title subtitle`).html(e.header(n)))}c.find(`div.modal-body`).empty().append(s),c.attr(`data-dtr-index`,n.index()).appendTo(`body`),c.classAdd(`is-active is-clipped`),t.s(`.modal-close`).one(`click`,function(){c.classRemove(`is-active is-clipped`),o()}),t.s(`.modal-background`).one(`click`,function(){c.classRemove(`is-active is-clipped`),o()})}else if(c.isAttached()&&n.index()===c.attr(`data-dtr-index`))c.find(`div.modal-body`).empty().append(s);else return null;return!0}};

