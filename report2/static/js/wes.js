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

//HANDLE somatic_summaries
$('.data-coloured.Total_Mutations').on('click', function() {
    var table_id = this.closest('table').id;
    var sample_id = $(this.closest('tr')).find('th').attr('data-original-sn');
    var tmp = JSON.parse($("#wes_resources").text());
    //console.log(table_id);
    //console.log(sample_id);
    //console.log(tmp['somatic_summary'][sample_id]);

    var modal = $('#wesSubModal');
    var data = new Object();
    data[sample_id] = tmp['somatic_summary'][sample_id]
    table_submodal(data, modal);

});
//Change roll-over cursor
$('.data-coloured.Total_Mutations').css('cursor', 'pointer');

//Builds a table in the modal and displays it
function table_submodal(data, modal) {
   console.log(data);
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
