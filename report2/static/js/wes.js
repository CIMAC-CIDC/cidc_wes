/* Len Taing 2020 (TGBTG) */

/* Function to toggle sections in the main panel */
//var last_section = 'meta';
function toggler(divId) {
    if (last_section != null) {
        $("#" + last_section).toggle();
    }
    $("#" + divId).toggle();
    last_section = divId;
}

//HANDLE table cells image clicks
$('.wes-image-modal').on('click',function(){
    $('#wesImageModal').modal({show:true});
    $('#wesImageModal_img').attr("src", this.src);
});

$('.data-coloured.Total_Mutations').on('click', function() {
    var table_id = this.closest('table').id;
    var sample_id = $(this.closest('tr')).find('th').attr('data-original-sn');
    console.log(table_id);
    console.log(sample_id);
    var tmp = JSON.parse($("#wes_resources").text());
    console.log(tmp['somatic_summary'][sample_id]);
});

var wes_resources = JSON.parse($("#wes_resources").text());

//Should generalize this fn- ColName, resource, handler
function handlerFactory(colName, resource, handler) {
    $(".data-coloured."+colName).on("click", function() {
	var table_id = this.closest('table').id;
	var sample_id = $(this.closest('tr')).find('th').attr('data-original-sn');
	var modal = $('#wesSubModal');
	var data = new Object();
	data[sample_id] = wes_resources[resource][sample_id];
	handler(data, modal);
    });
    $(".data-coloured."+colName).css('cursor', 'pointer');   
}
handlerFactory('Total_Mutations', 'somatic_summary', somaticSummary_submodal);
handlerFactory('TiTv', 'ti_tv', tiTv_submodal);
handlerFactory('TMB', 'tmb', somaticSummary_submodal);

//Builds a table in the modal and displays it
function somaticSummary_submodal(data, modal) {
   //console.log(data);
   var modal_body = $(modal).find('.modal-body');
   //Clear contents
   modal_body.empty()

   //BUILD table
   var content = "<table class=\"table table-condensed mqc_table\">";
   //add header
   content += "<thead><tr>";
   //Sample names is first col
   content += "<th>&nbsp;</th>";
   var first_elm = Object.keys(data)[0];
   var hdr = Object.keys(data[first_elm]);
   $.each(hdr, function(i,val) { content += "<th>"+val+"</th>";})
   content += "</tr></thead>";
   //Add content
   const samples = Object.keys(data);
   for (s of samples) {
      content += "<tr><td>"+s+"</td>";
      var items = Object.values(data[s]);
      for (v of items) {
         content += "<td>"+v+"</td>";
      }
   }
   content += "</tr>";
   content += "</table>";
   modal_body.append(content);
   modal.modal({show:true});
}

function tiTv_submodal(data, modal) {
    console.log(data);
    var modal_body = $(modal).find('.modal-body');
    //Clear contents
    modal_body.empty()

    //BUILD table
    var content = "<table class=\"table table-condensed mqc_table\">";
    //add header
    content += "<thead><tr>";
    content += "<th>&nbsp;</th>";
    const ACGT = ["A","C","G","T"];
    for (b of ACGT) { content += "<th>"+b+"</th>";}
    content += "</tr></thead>";
    //Add content
   var first_elm = Object.values(data)[0];
    for (b of ACGT) {
        content += "<tr><th>"+b+"</th>";
        console.log(b);
        console.log(first_elm[b]);
        var tmp = Object.values(first_elm[b]);
        for (v of tmp) { content += "<td>"+v+"</td>";}
    }
    content += "</tr>";
    content += "</table>";
    modal_body.append(content);
    modal.modal({show:true});
}
