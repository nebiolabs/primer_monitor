// datatables.net-buttons-bm@4.1.2 downloaded from https://ga.jspm.io/npm:datatables.net-buttons-bm@4.1.2/js/buttons.bulma.mjs

import e from"datatables.net-bm";export{default}from"datatables.net-bm";import"datatables.net-buttons";
/*! Buttons Bulma styling 4.1.2 for DataTables
* Copyright (c) SpryMedia Ltd - datatables.net/license
*/
var t=e.Dom;e.util.object.assignDeep(e.Buttons.defaults,{dom:{container:{className:`dt-buttons field is-grouped`},button:{className:`button`,active:`is-active`,disabled:`is-disabled`,dropHtml:`<span class="icon is-small"><i class="fa fa-angle-down" aria-hidden="true"></i></span>`,dropClass:``},collection:{action:{tag:`div`,className:`dropdown-content`},button:{tag:`a`,className:`dt-button dropdown-item`,active:`dt-button-active`,disabled:`is-disabled`,spacer:{className:`dropdown-divider`,tag:`hr`}},closeButton:!1,container:{className:`dt-button-collection dropdown dropdown-menu`,content:{className:`dropdown-content`}}},split:{action:{tag:`button`,className:`dt-button-split-drop-button button`,closeButton:!1},dropdown:{tag:`button`,className:`button`,closeButton:!1,align:`split-left`,splitAlignClass:`dt-button-split-left`},wrapper:{tag:`div`,className:`dt-button-split dropdown-trigger buttons has-addons`,closeButton:!1}}},buttonCreated:function(e,n){return e.buttons&&(e._collection=t.c(`div`).classAdd(`dropdown-menu`).append(e._collection)),n}});

