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

