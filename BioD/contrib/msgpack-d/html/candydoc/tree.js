/* This file is a part of CanDyDOC fileset.
   File is written by Victor Nakoryakov and placed into the public domain.

   This file is javascript with classes that represents native style tree control. */
   
var pmNone = 0;
var pmPlus = 1;
var pmMinus = 2;

var hlNone = 0;
var hlGrey = 1;
var hlSelected = 2;

function TreeView(hrefMode)
{
    this.domEntry = document.createElement("div");
    this.children = new Array();
    this.selection = null;
    this.hrefMode = hrefMode;
    
    this.createBranch = function(text, iconSrc)
    {
        var root = new TreeNode(text, iconSrc, this.hrefMode);
        root.owner = this;
        this.children[ this.children.length ] = root;
        this.domEntry.appendChild( root.domEntry );
        return root;
    }
    
    this.branch = function(text)
    {
        var ret = null;
        for (var i = 0; i < this.children.length; ++i)
            if (this.children[i].textElement.data == text)
            {
                ret = this.children[i];
                break;
            }
            
        return ret;
    }
    
    this.domEntry.style.fontSize = "10px";
    this.domEntry.style.cursor = "default";
    this.domEntry.style.whiteSpace = "nowrap";
}

var idCounter = 0;
function TreeNode(text, iconSrc, hrefMode)
{
    this.id             = idCounter++;
    this.parentNode     = null;
    this.children       = new Array();
    this.domEntry       = document.createElement("div");
    this.icon           = document.createElement("img");
    this.textElement    = document.createTextNode(text);
    this.textSpan       = document.createElement("span");
    this.lineDiv        = document.createElement("div");
    this.hierarchyImgs  = new Array();
    this.onclick        = null;
    
    function createIcon()
    {
        var img = document.createElement("img");
        img.style.verticalAlign = "middle";
        img.style.position = "relative";
        img.style.top = "-1px";
        img.width = 16;
        img.height = 16;
        return img;
    }
    
    function createHierarchyImage()
    {
        var img = createIcon();
        img.pointsTop = false;
        img.pointsBottom = false;
        img.pointsRight = false;
        img.pmState = pmNone;
        return img;
    }
    
    function genHierarchyImageSrc(hierarchyImg)
    {
        var name = "";
        if (hierarchyImg.pointsTop)
            name += "t";
            
        if (hierarchyImg.pointsBottom)
            name += "b";
            
        if (hierarchyImg.pointsRight)
            name += "r";
            
        if (hierarchyImg.pmState == pmPlus)
            name += "p";
        else if (hierarchyImg.pmState == pmMinus)
            name += "m";
        
        if (name == "")
            name = "shim";
        
        return "candydoc/img/tree/" + name + ".gif";
    }
    
    function setSrc(icon, src)
    {
        icon.src = src;
        // After src change width and height are reseted in IE.
        // Bug workaround:
        icon.width = 16;
        icon.height = 16;
    }
    
    this.createChild = function(text, iconSrc)
    {
        var child = new TreeNode(text, iconSrc, this.owner.hrefMode);
        this.children[ this.children.length ] = child;
        this.domEntry.appendChild( child.domEntry );
        child.parentNode = this;
        child.owner = this.owner;
        
        // insert hierarchy images according to deepness level
        // of created child.
        
        if (this.children.length > 1)
        {
            // there were already added child before. So copy `level-1`
            // hierarchy images from it.
            
            var prevAddedChild = this.children[ this.children.length - 2 ];
            
            for (var i = 0; i < prevAddedChild.hierarchyImgs.length - 1; ++i)
            {
                var prevAddedChildImg = prevAddedChild.hierarchyImgs[i];
                var img = createHierarchyImage();
                setSrc(img, prevAddedChildImg.src);
                img.pointsTop = prevAddedChildImg.pointsTop;
                img.pointsBottom = prevAddedChildImg.pointsBottom;
                img.pointsRight = prevAddedChildImg.pointsRight;
                img.pmState = prevAddedChildImg.pmState;
                
                child.hierarchyImgs[ child.hierarchyImgs.length ] = img;
                child.lineDiv.insertBefore(img, child.icon);
            }
            
            // change last hierarchy image of prevAddedChild from |_ to |-
            var lastHierarchyImg = prevAddedChild.hierarchyImgs[ prevAddedChild.hierarchyImgs.length - 1 ];
            lastHierarchyImg.pointsBottom = true;
            setSrc(lastHierarchyImg, genHierarchyImageSrc(lastHierarchyImg));
                        
            // change hierarchy images of prevAddedChild's children on it's last
            // level to |
            prevAddedChild.addHierarchyTBLine(prevAddedChild.hierarchyImgs.length - 1);
        }
        else
        {
            // this is a first child. So copy `level-2`
            // hierarchy images from parent, i.e. this.
            
            for (var i = 0; i < this.hierarchyImgs.length - 1; ++i)
            {
                var parentImg = this.hierarchyImgs[i];
                var img = createHierarchyImage();
                setSrc(img, parentImg.src);
                img.pointsTop = parentImg.pointsTop;
                img.pointsBottom = parentImg.pointsBottom;
                img.pointsRight = parentImg.pointsRight;
                img.pmState = parentImg.pmState;
                
                child.hierarchyImgs[ child.hierarchyImgs.length ] = img;
                child.lineDiv.insertBefore(img, child.icon);
            }
            
            if (this.hierarchyImgs.length > 0) // we are not root
            {
                // change last hierarchy image of parent (i.e. this): add minus to it
                var lastHierarchyImg = this.hierarchyImgs[ this.hierarchyImgs.length - 1];
                lastHierarchyImg.pmState = pmMinus;
                setSrc(lastHierarchyImg, genHierarchyImageSrc(lastHierarchyImg));
                lastHierarchyImg.owner = this;
                lastHierarchyImg.onclick = new Function("e", "this.owner.processPMClick(e);");
                
                // make decision on image on `level-1`. It depends on parent's (ie this)
                // image on same level.
                var parentL1HierarchyImg = lastHierarchyImg;
                var l1HierarchyImg = createHierarchyImage();
                if (parentL1HierarchyImg.pointsBottom)
                {
                    l1HierarchyImg.pointsTop = true;
                    l1HierarchyImg.pointsBottom = true;
                }
                setSrc(l1HierarchyImg, genHierarchyImageSrc(l1HierarchyImg));
                child.hierarchyImgs[ child.hierarchyImgs.length ] = l1HierarchyImg;
                child.lineDiv.insertBefore(l1HierarchyImg, child.icon);
            }
        }
        
        // in any case on last level our child will have icon |_
        var img = createHierarchyImage();
        img.pointsTop = true;
        img.pointsRight = true;
        setSrc(img, genHierarchyImageSrc(img));
        
        child.hierarchyImgs[ child.hierarchyImgs.length ] = img;
        child.lineDiv.insertBefore(img, child.icon);
        
        return child;
    }
    
    this.lastChild = function()
    {
        return this.children[ this.children.length - 1 ];
    }
    
    this.child = function(text)
    {
        var ret = null;
        for (var i = 0; i < this.children.length; ++i)
            if (this.children[i].textElement.data == text)
            {
                ret = this.children[i];
                break;
            }
            
        return ret;
    }
    
    this.addHierarchyTBLine = function(level)
    {
        for (var i = 0; i < this.children.length; ++i)
        {
            var img = this.children[i].hierarchyImgs[level];
            img.pointsTop = true;
            img.pointsBottom = true;
            setSrc(img, genHierarchyImageSrc(img));
            this.children[i].addHierarchyTBLine(level);
        }
    }
    
    this.expand = function()
    {
        var img = this.hierarchyImgs[ this.hierarchyImgs.length - 1 ];
        
        if (img.pmState == pmPlus)
        {
            img.pmState = pmMinus;
            setSrc(img, genHierarchyImageSrc(img));
            
            for (var i = 0; i < this.children.length; ++i)
                this.children[i].domEntry.style.display = "";
        }
    }
    
    this.collapse = function()
    {
        var img = this.hierarchyImgs[ this.hierarchyImgs.length - 1 ];
        
        if (img.pmState == pmMinus)
        {
            img.pmState = pmPlus;
            setSrc(img, genHierarchyImageSrc(img));
            
            for (var i = 0; i < this.children.length; ++i)
                this.children[i].domEntry.style.display = "none";
        }
    }
    
    this.toggle = function()
    {
        var img = this.hierarchyImgs[ this.hierarchyImgs.length - 1 ];
        if (img.pmState == pmMinus)
            this.collapse();
        else
            this.expand();
    }
    
    this.select = function()
    {
        if (this.owner.selection != this)
        {
            if (this.owner.selection)
                this.owner.selection.setHighlight(hlNone);
                
            this.owner.selection = this;
            this.setHighlight(hlSelected);
        }
    }
    
    this.setHighlight = function(mode)
    {
        if (mode == hlNone)
        {
            this.textSpan.style.backgroundColor = "";
            this.textSpan.style.color = "";
            this.textSpan.style.border = "";
        }
        else if (mode == hlGrey)
        {
            this.textSpan.style.backgroundColor = "#aaaaaa";
            this.textSpan.style.color = "";
            this.textSpan.style.border = "";
        }
        else if (mode == hlSelected)
        {
            this.textSpan.style.backgroundColor = "3399cc";
            this.textSpan.style.color = "white";
            this.textSpan.style.border = "dotted 1px red";
        }
    }
    
    this.setOnclick = function(proc)
    {
        this.onclick = proc;
    }
    
    this.setRef = function(url)
    {
        if (this.anchor)
            this.anchor.href = url;
    }
    
    this.processPMClick = function(e)
    {
        this.toggle();
        
        // prevent this line selection, stop bubbling
        if (e)
            e.stopPropagation(); // Mozilla way
        if (window.event)
            window.event.cancelBubble = true; // IE way
    }
    
    this.processOnclick = function()
    {
        this.select();
        if (this.onclick instanceof Function)
            this.onclick();
    }
    
    ///////////////////////////////////////////////////////////////////////////
    if (iconSrc)
        this.icon.src = iconSrc;
    else
    {
        this.icon.width = 0;
        this.icon.height = 0;
    }
    
    this.icon.style.verticalAlign = "middle";
    this.icon.style.position = "relative";
    this.icon.style.top = "-1px";
    this.icon.style.paddingRight = "2px";
    
    if (!hrefMode)
    {
        this.textSpan.appendChild( this.textElement );
    }
    else
    {
        this.anchor = document.createElement("a");
        this.anchor.appendChild( this.textElement );
        this.textSpan.appendChild( this.anchor );
    }
    
    this.lineDiv.appendChild( this.icon );
    this.lineDiv.appendChild( this.textSpan );
    this.domEntry.appendChild( this.lineDiv );
    
    this.lineDiv.owner = this;
    
    if (!hrefMode)
        this.lineDiv.onclick = new Function("this.owner.processOnclick();");
}
