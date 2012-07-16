import sambamba.view;
import sambamba.index;
import sambamba.merge;
import sambamba.sort;
import sambamba.flagstat;

import std.stdio;

void printUsage() {
    writeln("Usage: sambamba [command] [args...]");
    writeln();
    writeln("    Available commands: 'view', 'index', 'merge', 'sort', 'flagstat'.");
    writeln("    To get help on a particular command, just call it without args.");
    // TODO: man pages
}

int main(string[] args) {
    if (args.length == 1) {
        printUsage();
        return 1;
    }

    auto _args = args[0] ~ args[2 .. $];

    switch (args[1]) {
        case "view":     return view_main(_args);
        case "index":    return index_main(_args);
        case "merge":    return merge_main(_args);
        case "sort":     return sort_main(_args);
        case "flagstat": return flagstat_main(_args);
        default: 
            printUsage();
            return 1;
    }
}
