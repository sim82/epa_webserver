// Place your application-specific JavaScript functions and classes here
// This file is automatically included by javascript_include_tag :defaults

function optionalButton(id,name){
    var vis = document.getElementById(id).style.visibility;
    if (vis == 'visible'){
        document.getElementById(id).style.visibility = 'hidden';
        document.getElementsByName(name)[0].value='show more options';  
    }
    else{
        document.getElementById(id).style.visibility = 'visible';
        document.getElementsByName(name)[0].value='hide more options';
    }
                  
}

function optionalCheckboxElements(id1,name1,id2,name2){
    var vis = document.getElementById(name1);
    if (vis.checked){
        document.getElementById(id1).style.visibility = 'visible';
        document.getElementById(id2).style.visibility = 'hidden';
        document.getElementById(name2).checked = false;
    }
    else{
        document.getElementById(id1).style.visibility = 'hidden';
    }
}

function showHideCheckboxElements(checkboxID, elementsID){
    var vis = document.getElementById(checkboxID);
    if (vis.checked){
        document.getElementById(elementsID).style.visibility = 'visible';
    }
    else {
        document.getElementById(elementsID).style.visibility = 'hidden';
    }
}
        

        

