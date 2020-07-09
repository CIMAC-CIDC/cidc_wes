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

//DELETE the following
/* function to handle the dropdown select --not fully implemented */
$("#menu-toggle").click(function(e) {
    e.preventDefault();
    $("#wrapper").toggleClass("toggled");
    currVal = $("#menu-toggle").text();
    nextVal = (currVal=='X') ? '>' : 'X';
    $("#menu-toggle").text(nextVal);
});

