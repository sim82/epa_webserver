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
        

        

