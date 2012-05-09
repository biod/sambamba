import samheader;
import reference;
mixin(import("utils/scaffold.d"));

import std.stdio;

void createScaffold(T)() {
    File("bindings/scaffolds/" ~ T.stringof ~ ".rb", "w")
        .write(toRuby!T);
}

void main() {
    createScaffold!SqLine();
    createScaffold!RgLine();
    createScaffold!PgLine();
    createScaffold!SamHeader();
    createScaffold!ReferenceSequenceInfo();
}
