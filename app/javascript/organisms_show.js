import 'init_jquery';
import "igv";

let config = {};

function initBrowser() {

    const browserConfig =
        {
            reference: {
                "id": config['reference_accession'],
                "name": config['organism_name'] + " (" + config['reference_accession'] + ")",
                "fastaURL": config['data_server'] + "/" + config['organism_slug'] + "/ref/" + config['organism_slug'] + ".fasta",
                "indexURL": config['data_server'] + "/" + config['organism_slug'] + "/ref/" + config['organism_slug'] + ".fasta.fai",
                tracks: [
                    {
                        "name": "Genes",
                        "type": "annotation",
                        "url": config['data_server'] + "/" + config['organism_slug'] + "/ref/" + config['organism_slug'] + ".gff3",
                        "format": "gff3",
                        "filterTypes": ['CDS', 'mature_protein_region_of_CDS', 'region', 'stem_loop', 'five_prime_UTR', 'three_prime_UTR'],
                        "displayMode": "EXPANDED",
                        "colorBy": "gbkey",
                        "colorTable": {
                            "Gene": "rgb(0,190,0)",
                        }
                    }
                ]
            }
        };

    $('.igv_div').children('.igv-container, .igv-message').remove();

    let igvBrowser = null;

    const browser_div = document.getElementById("igv");
    igv.createBrowser(browser_div, browserConfig).then(function (igvBrowser) {
        if (config['variants_exist']) {
            const variantsTrack = {
                "name": "Variants",
                "url": config['data_server'] + "/" + config['organism_slug'] + "/lineage_variants/all.bed",
                "format": "bed",
                "color": "#575757",
                "displayMode": "COLLAPSED",
                "autoHeight": true
            }
            igvBrowser.loadTrack(variantsTrack);
        }
    });



}

$(document).ready(function(){
    config = JSON.parse($('#config')[0].innerHTML);
    initBrowser();
});