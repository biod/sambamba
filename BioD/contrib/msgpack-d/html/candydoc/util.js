/* This file is a part of CanDyDOC fileset.
   File is written by Victor Nakoryakov and placed into the public domain.

   This file is javascript with cross-browser utility functions. */

function getLeft(elem)
{
    var ret = 0;
	while (elem.offsetParent)
	{
		ret += elem.offsetLeft;
		elem = elem.offsetParent;
	}

	return ret;
}

function getTop(elem)
{
    var ret = 0;
	while (elem.offsetParent)
	{
		ret += elem.offsetTop;
		elem = elem.offsetParent;
	}

	return ret;
}

function getWindowHeight()
{
    var ret = 0;
    if (typeof(window.innerHeight) == "number")
        ret = window.innerHeight;
    else if (document.documentElement && document.documentElement.clientHeight)
        ret = document.documentElement.clientHeight;
    else if (document.body && document.body.clientHeight)
        ret = document.body.clientHeight;
    
    return ret;
}
