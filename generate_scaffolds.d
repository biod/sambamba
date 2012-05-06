import samheader;
mixin(import("utils/scaffold.d"));

import std.stdio;

void main() {
    File("bindings/scaffolds/SqLine.rb", "w").write(toRuby!SqLine);
    File("bindings/scaffolds/RgLine.rb", "w").write(toRuby!RgLine);
    File("bindings/scaffolds/PgLine.rb", "w").write(toRuby!PgLine);
    File("bindings/scaffolds/SamHeader.rb", "w").write(toRuby!SamHeader);
}
