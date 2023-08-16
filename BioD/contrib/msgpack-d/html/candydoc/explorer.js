/* This file is a part of CanDyDOC fileset.
   File is written by Victor Nakoryakov and placed into the public domain.

   This file is javascript with classes that represents explorer window.
   And things related to navigation. */
   
var explorer = new Explorer();

///////////////////////////////////////////////////////////////////////////////
// Current symbol marker class constructor
///////////////////////////////////////////////////////////////////////////////
function Marker()
{
    this.top    = document.createElement("div");
    this.middle = document.createElement("div");
    this.bottom = document.createElement("div");
    this.container = document.createElement("div");
    
    this.setTo = function(term)
    {
        // find definition related to `term`
        var def = term.nextSibling;
        while (def && def.nodeName != "DD")
            def = def.nextSibling;
            
        var defHeight = 0;
        var childrenHeight = 0; // children of current declaration
        if (def)
        {
            defHeight = def.offsetHeight;
            var child = def.firstChild;
            
            // traverse until DL tag, until children definition
            while (child && child.nodeName != "DL")
                child = child.nextSibling;
                
            if (child)
                childrenHeight = child.offsetHeight;
        }
        
        this.top.style.height = term.offsetHeight;
        this.middle.style.height = defHeight - childrenHeight;
        this.bottom.style.height = childrenHeight;
        
        if (childrenHeight == 0)
            this.bottom.style.display = "none";
        else
            this.bottom.style.display = "";
        
        this.container.style.left = getLeft(term) - 8;
        this.container.style.top = getTop(term);
        this.container.style.display = "";
    }
        
    ///////////////////////////////////////////////////////////////////////////
    this.container.style.position = "absolute";
    this.container.style.display = "none";
    
    this.top.className = "markertop";
    this.middle.className = "markermiddle";
    this.bottom.className = "markerbottom";
    
    this.container.appendChild(this.top);
    this.container.appendChild(this.middle);
    this.container.appendChild(this.bottom);
    
    //document.body.appendChild( this.container );
    
    // Workaround bug in IE 5/6. We can not append anything to document body until
    // full page load.
    window.marker = this;
    if (window.addEventListener)
        window.addEventListener("load", new Function("document.body.appendChild( window.marker.container );"), false);
    else if (window.attachEvent)
        window.attachEvent("onload", new Function("document.body.appendChild( window.marker.container );"));
}

///////////////////////////////////////////////////////////////////////////////
// Outline class constructor
///////////////////////////////////////////////////////////////////////////////
function Outline()
{
    this.tree           = new TreeView();
    this.mountPoint     = null;    
    this.writeEnabled   = false;
    this.marker         = new Marker();        
    this.classRegExp    = new RegExp;
    this.structRegExp   = new RegExp;
    this.enumRegExp     = new RegExp;
    this.templateRegExp = new RegExp;
    this.aliasRegExp    = new RegExp;
    this.funcRegExp     = new RegExp;
    
    this.incSymbolLevel = function()
    {
        if (this.mountPoint == null)
            this.mountPoint = this.tree.children[ 0 ];
        else
            this.mountPoint = this.mountPoint.lastChild();
    }

    this.decSymbolLevel = function()
    {
        // place icons near items according to extracted below type
        for (var i = 0; i < this.mountPoint.children.length; ++i)
        {
            child = this.mountPoint.children[i];
            var term = child.termRef;
            
            // find first span node
            var n = term.firstChild;
            while (n && n.nodeName != "SPAN")
                n = n.nextSibling;
            
            if (!n) // shouldn't happen
                continue;
            
            var iconSrc;
            if (n.firstChild.nodeName == "#text")
            {
                var text = n.firstChild.data; // text before declaration
                
                if ( this.classRegExp.test(text) )
                    iconSrc = "candydoc/img/outline/class.gif";
                else if ( this.structRegExp.test(text) )
                    iconSrc = "candydoc/img/outline/struct.gif";
                else if ( this.enumRegExp.test(text) )
                    iconSrc = "candydoc/img/outline/enum.gif";
                else if ( this.templateRegExp.test(text) )
                    iconSrc = "candydoc/img/outline/template.gif";
                else if ( this.aliasRegExp.test(text) )
                    iconSrc = "candydoc/img/outline/alias.gif";
                else // function or variable? check whether '(' ')' exists on the right
                {
                    var np = n.firstChild;
                    while (np && np.nodeName != "SCRIPT") // find our script "onDecl"
                        np = np.nextSibling;
                        
                    if (np && np.nextSibling && np.nextSibling.nodeName == "#text" &&
                        this.funcRegExp.test(np.nextSibling.data))
                    {
                        iconSrc = "candydoc/img/outline/func.gif";
                    }
                    else
                        iconSrc = "candydoc/img/outline/var.gif";
                }
            }
            else // enum member ?
                iconSrc = "candydoc/img/outline/var.gif";
                    
            child.icon.src = iconSrc;
            child.icon.width = 16;
            child.icon.height = 16;
        }
        
        this.mountPoint = this.mountPoint.parentNode;
    }

    this.addDecl = function(decl)
    {
        function getLastLeaf(elem)
        {
            if (elem.childNodes.length > 0)
                return getLastLeaf(elem.lastChild);
            else
                return elem;
        }

        function getCurrentTerm()
        {
            var ret = getLastLeaf( document.getElementById("content") );
            while (ret && ret.nodeName != "DT")
                ret = ret.parentNode;
            
            return ret;
        }
        
        if (this.writeEnabled)
        {
            var node = this.mountPoint.createChild(decl);
            node.termRef = getCurrentTerm();
            node.setOnclick( new Function("explorer.outline.mark(this.termRef);") );
        }
    }

    this.mark = function(term)
    {
        this.marker.setTo(term);
        window.scrollTo(0, getTop(term) - getWindowHeight() / 6);    
    }
        
    
    this.classRegExp.compile("(.*\b)?class(\b.*)?");
    this.structRegExp.compile("(.*\b)?struct(\b.*)?");
    this.enumRegExp.compile("(.*\b)?enum(\b.*)?");
    this.templateRegExp.compile("(.*\b)?template(\b.*)?");
    this.aliasRegExp.compile("(.*\b)?alias(\b.*)?");
    this.funcRegExp.compile(/.*\(.*/);    
}




///////////////////////////////////////////////////////////////////////////////
// Package explorer class constructor
///////////////////////////////////////////////////////////////////////////////
function PackageExplorer()
{
    this.tree = new TreeView(true);    
    
    this.addModule = function(mod)
    {
        var moduleIco = "candydoc/img/outline/module.gif";
        var packageIco = "candydoc/img/outline/package.gif";

        var path = mod.split("\.");
        var node = this.tree.branch(path[0]);
        if ( !node )
            node = this.tree.createBranch(path[0], (path.length == 1) ? moduleIco : packageIco);
        
        for (var i = 1; i < path.length; ++i)
        {
            var prev = node;
            node = node.child(path[i]);
            if (!node)
                node = prev.createChild(path[i], (path.length == i + 1) ? moduleIco : packageIco);
                
            if (path.length == i + 1)
                node.setRef(path[i] + ".html");
        }
    }
}



///////////////////////////////////////////////////////////////////////////////
// Explorer class constructor
///////////////////////////////////////////////////////////////////////////////
function Explorer()
{
    this.outline         = new Outline();
    this.packageExplorer = new PackageExplorer();
    this.tabs            = new Array();
    this.tabCount        = 0;
    
    this.initialize = function(moduleName)
    {
        this.tabArea = document.getElementById("tabarea");
        this.clientArea = document.getElementById("explorerclient");
        
        // prevent text selection
        this.tabArea.onmousedown = new Function("return false;");
        this.tabArea.onclick = new Function("return true;");
        this.tabArea.onselectstart = new Function("return false;");
        this.clientArea.onmousedown = new Function("return false;");
        this.clientArea.onclick = new Function("return true;");
        this.clientArea.onselectstart = new Function("return false;");
        
        this.outline.tree.createBranch( moduleName, "candydoc/img/outline/module.gif" );
        
        // create tabs
        this.createTab("Outline", this.outline.tree.domEntry);
        this.createTab("Package", this.packageExplorer.tree.domEntry);
    }
    
    this.createTab = function(name, domEntry)
    {
        var tab = new Object();
        this.tabs[name] = tab;
        this.tabCount++;
        
        tab.domEntry = domEntry;
        tab.labelSpan = document.createElement("span");
        
        if (this.tabCount > 1)
        {
            tab.labelSpan.className = "inactivetab";
            tab.domEntry.style.display = "none";
        }
        else
        {
            tab.labelSpan.className = "activetab";
            tab.domEntry.style.display = "";
        }
        
        tab.labelSpan.appendChild( document.createTextNode(name) );
        tab.labelSpan.owner = this;
        tab.labelSpan.onclick = new Function("this.owner.setSelection('" + name + "');");
        
        this.tabArea.appendChild( tab.labelSpan );
        this.clientArea.appendChild( domEntry );
    }
    
    this.setSelection = function(tabName)
    {
        for (name in this.tabs)
        {
            this.tabs[name].labelSpan.className = "inactivetab";
            this.tabs[name].domEntry.style.display = "none";
        }
        
        this.tabs[tabName].labelSpan.className = "activetab";
        this.tabs[tabName].domEntry.style.display = "";
    }
}
