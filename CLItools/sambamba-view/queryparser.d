import filter;
import utils.pratt_parser;

import std.array;
import std.algorithm;
import std.conv;
import std.ascii;
import std.typetuple;
import std.exception;

interface Node {
    string toString() const;
}

class TagNameNode : Node {
    string tagname;
    this(string tagname) {
        this.tagname = tagname;
    }

    string toString() const {
        return "[" ~ tagname ~ "]";
    }
}

class IntegerNode : Node {
    long value;
    this(long value) {
        this.value = value;
    }
    string toString() const {
        return to!string(value);
    }
}

class IntegerFieldNode : Node {
    string fieldname;
    this(string fieldname) {
        this.fieldname = fieldname;
    }
    string toString() const {
        return fieldname;
    }
}

class StringNode : Node {
    string value;
    this(string value) {
        this.value = value;
    }
    string toString() const {
        return "'" ~ replace(value, "'", "\\'") ~ "'";
    }
}

class StringFieldNode : Node {
    string fieldname;
    this(string fieldname) {
        this.fieldname = fieldname;
    }
    string toString() const {
        return fieldname;
    }
}

abstract class ConditionNode : Node {
    Filter condition;
    abstract string toString() const;
}

/// LeftNodeType refers to TagNameNode, IntegerFieldNode, or StringFieldNode
/// RightNodeType refers to IntegerNode or StringNode
class ComparisonNode(T, LeftNodeType, RightNodeType, alias Filter) : ConditionNode {
    private {
        LeftNodeType _node;
        RightNodeType _valuenode;
        string _op;
    }
    this(string op, LeftNodeType node, RightNodeType valuenode) {
        _node = node;
        _valuenode = valuenode;
        _op = op;
        static if (is(LeftNodeType == TagNameNode)) {
            auto n = node.tagname;
        } else {
            auto n = node.fieldname;
        }
        auto v = valuenode.value;
        switch (op) {
            case ">": condition = new Filter!">"(n, v); break;
            case "<": condition = new Filter!"<"(n, v); break;
            case ">=": condition = new Filter!">="(n, v); break;
            case "<=": condition = new Filter!"<="(n, v); break;
            case "!=": condition = new Filter!"!="(n, v); break;
            case "==": condition = new Filter!"=="(n, v); break;
            default: assert(0); // must be catched before calling constructor
        }
    }

    override string toString() const {
        return _node.toString() ~ " " ~ _op ~ " " ~ _valuenode.toString();
    }
}

alias ComparisonNode!(long, TagNameNode, IntegerNode, IntegerTagFilter) IntegerTagConditionNode;
alias ComparisonNode!(string, TagNameNode, StringNode, StringTagFilter) StringTagConditionNode;
alias ComparisonNode!(long, IntegerFieldNode, IntegerNode, IntegerFieldFilter) IntegerFieldConditionNode;
alias ComparisonNode!(string, StringFieldNode, StringNode, StringFieldFilter) StringFieldConditionNode;

class OrConditionNode : ConditionNode {
    private ConditionNode _a, _b;
    this(Node a, Node b) {
        _a = cast(ConditionNode)a;
        _b = cast(ConditionNode)b;
        enforce(_a !is null, "'" ~ a.toString() ~ "' is not a condition");
        enforce(_a !is null, "'" ~ b.toString() ~ "' is not a condition");
        condition = new OrFilter(_a.condition, _b.condition);
    }

    override string toString() const {
        return "(" ~ _a.toString()  ~ ") or (" ~ _b.toString() ~ ")";
    }
}

class AndConditionNode : ConditionNode {
    private ConditionNode _a, _b;
    this(Node a, Node b) {
        _a = cast(ConditionNode)a;
        _b = cast(ConditionNode)b;
        enforce(_a !is null, "'" ~ a.toString() ~ "' is not a condition");
        enforce(_a !is null, "'" ~ b.toString() ~ "' is not a condition");
        condition = new AndFilter(_a.condition, _b.condition);
    }

    override string toString() const {
        return "(" ~ _a.toString()  ~ ") and (" ~ _b.toString() ~ ")";
    }
}

class NegateConditionNode : ConditionNode {
    private ConditionNode _a;
    this(Node a) {
        _a = cast(ConditionNode)a;
        enforce(_a !is null, "'" ~ a.toString() ~ "' is not a condition");
        condition = new NotFilter(_a.condition);
    }

    override string toString() const {
        return "not (" ~ _a.toString() ~ ")";
    }
}

class FlagConditionNode : ConditionNode {
    private string _flagname;
    this(in string flagname) {
        _flagname = flagname;
        switch (flagname) {
            case "paired":
                condition = new FlagFilter!"is_paired"(); break;
            case "proper_pair":
                condition = new FlagFilter!"proper_pair"(); break;
            case "unmapped":
                condition = new FlagFilter!"is_unmapped"(); break;
            case "mate_is_unmapped":
                condition = new FlagFilter!"mate_is_unmapped"(); break;
            case "reverse_strand":
                condition = new FlagFilter!"is_reverse_strand"(); break;
            case "mate_is_reverse_strand":
                condition = new FlagFilter!"mate_is_reverse_strand"(); break;
            case "first_of_pair":
                condition = new FlagFilter!"is_first_of_pair"(); break;
            case "second_of_pair":
                condition = new FlagFilter!"is_second_of_pair"(); break;
            case "secondary_alignment":
                condition = new FlagFilter!"is_secondary_alignment"(); break;
            case "failed_quality_control":
                condition = new FlagFilter!"failed_quality_control"(); break;
            case "duplicate":
                condition = new FlagFilter!"is_duplicate"(); break;
            default:
                throw new Exception("unknown flag '" ~ flagname ~ "'");
        }
    }

    override string toString() const {
        return _flagname; 
    }
}

final class QueryGrammar : Grammar!Node {
    this() {
        super("(end)");

        static auto makeScanner(string[] values) {
            return (in string str, size_t pos) {
                foreach (value; values) {
                    if (startsWith(str[pos .. $], value)) {
                        return pos + value.length;
                    }
                }
                return pos;
            };
        }

        auto flagnames = ["paired", "proper_pair", "unmapped", "mate_is_unmapped",
                          "reverse_strand", "mate_is_reverse_strand", "first_of_pair",
                          "second_of_pair", "secondary_alignment", "failed_quality_control",
                          "duplicate"];

        addSymbolToDict("(flag condition)", 0)
            .setScanner(makeScanner(flagnames))
            .setParser((in string str) { return cast(Node) new FlagConditionNode(str);;});

        auto integer_fields = ["ref_id", "position", "mapping_quality", 
                               "sequence_length", "mate_ref_id", "mate_position",
                               "template_length"];

        addSymbolToDict("(integer field)", 0)
            .setScanner(makeScanner(integer_fields))
            .setParser((in string str) { return cast(Node) new IntegerFieldNode(str);});
               
        auto string_fields = ["read_name"];

        addSymbolToDict("(string field)", 0)
            .setScanner(makeScanner(string_fields))
            .setParser((in string str) { return cast(Node) new StringFieldNode(str); });

        addSymbolToDict("(tag name)", 0)
            .setScanner(
                (in string str, size_t pos) {
                    auto s = str[pos .. $];
                    if (s.length == 0) return pos;
                    if (s[0] == '[' && s.length >= 4 && s[3] == ']') {
                        return pos + 4;
                    } else {
                        return pos;
                    }
                }
            ).setParser(
                (in string str) { return cast(Node) new TagNameNode(str[1..3]); }
            );

        addSymbolToDict("(integer)", 0)
            .setScanner(
                (in string str, size_t pos) {
                    if (pos >= str.length) return pos;
                    size_t i = pos;
                    if (str[i] == '-' || str[i] == '+') {
                        ++i;
                        if (str.length == 1) return pos;
                    }
                    while (i < str.length && isDigit(str[i])) ++i;
                    return i;
                }
            ).setParser(
                (in string str) { return cast(Node) new IntegerNode(to!long(str)); }
            );

        addSymbolToDict("(string)", 0)
            .setScanner(
                (in string str, size_t pos) {
                    if (pos >= str.length || str[pos] != '\'') return pos;
                    size_t i = pos + 1;
                    while (i < str.length) {
                        if (str[i] == '\\' && i + 1 < str.length && str[i+1] == '\'') {
                            i += 2;
                        } else if (str[i] == '\'') {
                            i += 1;
                            break;
                        } else {
                            i += 1;
                        }
                    }
                    if (i == str.length && str[$-1] != '\'') return pos;
                    return i;
                }
            ).setParser(
                (in string str) { 
                    assert(str[0] == '\'');
                    assert(str[$ - 1] == '\'');
                    auto strbuilder = appender!(char[])();
                    strbuilder.reserve(str.length - 2);
                    auto s = str[1 .. $ - 1];
                    for (size_t i = 0; i < s.length; i++) {
                        if (s[i] == '\\' && i + 1 < s.length && s[i+1] == '\'') {
                            strbuilder.put('\'');
                            i += 1;
                        } else {
                            strbuilder.put(s[i]);
                        }
                    }
                    return cast(Node) new StringNode(cast(string)(strbuilder.data));
                }
            );

        static auto makeComparisonOperator(string op) {
            return (Node a, Node b) { 
                auto integer = cast(IntegerNode) b;
                auto str = cast(StringNode) b;

                auto tag = cast(TagNameNode) a;
                if (integer !is null) {
                    if (tag !is null) {
                        auto node = new IntegerTagConditionNode(op, tag, integer);
                        return cast(Node) node;
                    } else {
                        auto field = cast(IntegerFieldNode) a;
                        if (field is null) {
                            throw new Exception("expected tag or integer field name instead of '"~
                                                a.toString() ~ "'"); 
                        }
                        return cast(Node) new IntegerFieldConditionNode(op, field, integer);
                    }
                } else if (str !is null) {
                    if (tag !is null) {
                        auto node = new StringTagConditionNode(op, tag, str);
                        return cast(Node) node;
                    } else {
                        auto field = cast(StringFieldNode) a;
                        if (field is null) {
                            throw new Exception("expected tag or string field name instead of '"~
                                                a.toString() ~ "'");
                        }
                        return cast(Node) new StringFieldConditionNode(op, field, str);
                    }
                } else {
                    assert(0); // no way to get here
                }
            };
        }

        foreach (op; [">", "<", ">=", "<=", "==", "!="])
            infix(op, 110, makeComparisonOperator(op));

        infix("and", 80, (Node a, Node b) { return cast(Node) new AndConditionNode(a, b); });
        infix("or", 60, (Node a, Node b) { return cast(Node) new OrConditionNode(a, b);  });
        prefix("not", 100, (Node a) { return cast(Node) new NegateConditionNode(a); });
        brackets("(", ")", int.max);
    }
}

unittest {
    import std.stdio;
    auto grammar = new QueryGrammar();
    assert(grammar.parse("unmapped and mate_is_unmapped").toString() == "(unmapped) and (mate_is_unmapped)");
    writeln(grammar.parse("[BQ] >= 'ab\\'asdga' and ([NM] < 12 or unmapped and [PS] == '321') or mapping_quality < 150 and read_name >= 'abcd'").toString());
}
