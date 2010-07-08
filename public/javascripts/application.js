// Place your application-specific JavaScript functions and classes here
// This file is automatically included by javascript_include_tag :defaults

function optionalButton(id,name){
    var vis = document.getElementById(id).style.display;
    if (vis == 'block'){
        document.getElementById(id).style.display = 'none';
        document.getElementsByName(name)[0].value='show more options';  
    }
    else{
        document.getElementById(id).style.display = 'block';
        document.getElementsByName(name)[0].value='hide more options';
    }
                  
}

function optionalCheckboxElements(id1,name1,id2,name2){
    var vis = document.getElementById(name1);
    if (vis.checked){
        document.getElementById(id1).style.display = 'block';
        document.getElementById(id2).style.display = 'none';
        document.getElementById(name2).checked = false;
    }
    else{
        document.getElementById(id1).style.display = 'none';
    }
}

function showHideCheckboxElements(checkboxID, elementsID){
    var vis = document.getElementById(checkboxID);
    if (vis.checked){
        document.getElementById(elementsID).style.display = 'block';
    }
    else {
        document.getElementById(elementsID).style.display = 'none';
    }
}

function showTooltip(id, event)
{
    var yoffset = parseInt(document.body.scrollTop);
    var xoffset = parseInt(document.body.scrollLeft);
    var mouseY = (event.clientY) ? event.clientY : event.pageY;
    var mouseX = (event.clientX) ? event.clientX : event.pageX;
    document.getElementById(id).style.top  = mouseY + 20 + yoffset+ "px";
    document.getElementById(id).style.left = mouseX + 10 + xoffset + "px";
    document.getElementById(id).style.visibility = "visible";
   }
function hideTooltip(id)
   {
      document.getElementById(id).style.visibility = "hidden";
   }
        

