/*
    This file is part of Sambamba.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
module sambamba.utils.common.pratt_parser;

import std.string;
import std.exception;
import std.conv;
import std.algorithm;
import std.ascii;

struct SourcePosition {
    size_t position;
    size_t line;
    size_t column;
}

class PrattParser(T, alias WSSkipper=DefaultWSSkipper) {

    private {
        alias Token!(T, WSSkipper) Tok;
        alias SymbolDict!(T, WSSkipper) Dict;
        alias PrattParser!(T, WSSkipper) Parser;
        alias TokenRange!(T, WSSkipper) Range;
    }

    this(string str, Dict symbols) {
        this.str = str;
        token_range = Range(str, symbols);
        token = next();
    }

    T parse(int rbp=0) {
        prev_token = token;
        if (!prev_token_initialized) {
            prev_token_initialized = true;
        }
        token = next();

        T left = prev_token.nud(this);
        while (rbp < token.lbp) {
            prev_token = token;
            token = next();
            left = prev_token.led(this, left);
        }
        return left;
    }

    auto next_token() {
        return token;
    }

    string next_token_as_string() @property const {
        return str[token.start_position .. token.end_position];
    }

    Parser advance() { 
        token = next();
        return this;
    }

    Parser advance(string s) {
        if (next_token_as_string() != s) {
            throw new Exception("unexpected character at position " ~ to!string(current_position.line) ~ ":" ~
                    to!string(current_position.column));
        }
        token = next();
        return this;
    }

    @property SourcePosition current_position() const {
        SourcePosition sp;
        auto tok = token_initialized ? token : prev_token;
        sp.position = prev_token_initialized ? tok.start_position : 0;
        sp.line = token_range.current_line;
        sp.column = sp.position - token_range.last_new_line + 1;
        return sp;
    }

    @property string code() const { return str; }

private:
    string str;

    Range token_range;
    Tok token = void;
    Tok prev_token = void;
    bool token_initialized = false;
    bool prev_token_initialized = false;

    auto next() {
        auto token = token_range.front();
        if (!token_initialized) {
            token_initialized = true;
        }
        token_range.popFront();
        return token;
    }
}

class Symbol(T, alias WSSkipper=DefaultWSSkipper) {

    alias size_t delegate(in string str, size_t pos) TokenScanner;
    alias T delegate(in string str) TokenParser;
    alias PrattParser!(T, WSSkipper) Parser;

    private TokenScanner scanner;
    private TokenParser parser;

    string id;
    int lbp;
    T delegate(Parser) nud;
    T delegate(Parser, T) led;

    this(string id, int lbp=0) {
        this.id = id;
        this.lbp = lbp;
    }

    size_t scan(in string str, size_t pos) const {
        if (scanner) {
            return scanner(str, pos);
        } else {
            size_t len = id.length;
            for (size_t i = 0; i < len; ++i) {
                if (i + pos >= str.length || id[i] != str[i + pos])
                    return pos;
            }
            return pos + len;
        }
    }

    T parse(in string str) const {
        if (parser) {
            return parser(str);
        } else {
            T obj;
            return obj;
        }
    }

    bool has_scanner() @property const { return scanner !is null; }
    bool has_parser() @property const { return parser !is null; }

    Symbol setScanner(TokenScanner s) { 
        scanner = s; 
        return this;
    }

    Symbol setParser(TokenParser p) {
        parser = p; 
        return this;
    }
}

class SymbolDict(T, alias WSSkipper) {
    alias Symbol!(T, WSSkipper) Sym;
    private Sym[string] dict;
    private string end_id_;

    this(string end_id) { 
        this.end_id_ = end_id;
        dict[end_id] = new Sym(end_id, int.min);
    }

    Sym opIndex(string id) {
        return dict[id];
    }

    Sym opIndexAssign(Sym sym, string id) {
        dict[id] = sym;
        return sym;
    }

    bool opBinaryRight(string op)(string id) if (op == "in") { 
        return ((id in dict) !is null); 
    }

    @property auto values() { return dict.byValue(); }
    @property Sym end_symbol() { return dict[end_id_]; }
    @property string end_id() { return end_id_; }
}

struct Token(T, alias WSSkipper=DefaultWSSkipper) {
    private alias Symbol!(T, WSSkipper) Sym;
    private alias PrattParser!(T, WSSkipper) Parser;

    private Sym symbol_;
    size_t start_position;
    size_t length;
    T value = void;
    bool has_value = false;

    this(Sym sym, size_t start=0, size_t end=0) {
        symbol_ = sym;
        start_position = start;
        length = end - start;
    }
   
    this(Sym sym, T value, size_t start=0, size_t end=0) {
        this(sym, start, end);
        this.value = value;
        has_value = true;
    }

    @property string id() const { return symbol_.id; }
    @property Sym symbol() { return symbol_; }
    @property int lbp() const { return symbol_.lbp; }
    @property size_t end_position() const { return start_position + length; }

    T nud(Parser parser) {
        if (has_value) return value;
        if (symbol_.nud is null) {
            throw new Exception("parsing error near line " ~ to!string(parser.current_position.line)
                                ~ ": expected prefix operator");
        }
        return symbol_.nud(parser);
    }

    T led(Parser parser, T left) {
        if (symbol_.led is null) {
            throw new Exception("parsing error near line " ~ to!string(parser.current_position.line)
                                ~ ": expected infix/postfix operator");
        }
        return symbol_.led(parser, left);
    }
} 

struct TokenRange(T, alias WSSkipper) {

    private alias Symbol!(T, WSSkipper) Sym;
    private alias PrattParser!(T, WSSkipper) Parser;
    private alias SymbolDict!(T, WSSkipper) Dict;
    private alias Token!(T, WSSkipper) Tok;

private:
    string str;
    Dict symbols;
    size_t start;
    size_t end;
    Sym match;
    bool match_found;
    size_t last_new_line_;
    size_t current_line_;
public:

    @property size_t current_line() const { return current_line_; }
    @property size_t last_new_line() const { return last_new_line_; }

    this(string str_, Dict symbols_) {
        str = str_;
        symbols = symbols_;
        start = end = last_new_line_ = 0;
        current_line_ = 1;
        popFront();
    }

    bool empty() {
        return false;
    }

    void popFront() {
        WSSkipper(str, start, last_new_line_, current_line_);

        if (start < str.length) {
            end = start;
            match_found = false;
            foreach (sym; symbols.values) {
                size_t p = sym.scan(str, start);

                if (p > end || (match_found && sym.lbp > match.lbp && p == end)) {
                    match = sym;
                    match_found = true;
                    end = p;
                }
            }
            if (!match_found) {
                throw new Exception("invalid symbol in input stream at position " ~ to!string(start));
            }
        }
    }

    auto front() {
        if (start >= str.length) {
            return Tok(symbols.end_symbol);
        }
        size_t old_start = start;
        start = end;

        if (match.has_parser) {
            return Tok(match, match.parse(str[old_start .. end]), old_start, end);
        } else {
            return Tok(match, old_start, end);
        }
    }
} 

private struct keep_symbol_lbp_t {}
static keep_symbol_lbp_t keep_symbol_lbp;

void DefaultWSSkipper(const char[] str, ref size_t start,
                                        ref size_t last_new_line,
                                        ref size_t current_line)
{
    while (start < str.length && std.ascii.isWhite(str[start])) {
        if (str[start] == '\n') {
            last_new_line = start;
            ++current_line;
        }
        ++start;
    }
};

class Grammar(T, alias WSSkipper=DefaultWSSkipper) {
    protected alias Symbol!(T, WSSkipper) Sym;
    protected alias PrattParser!(T, WSSkipper) Parser;
    protected alias SymbolDict!(T, WSSkipper) Dict;

    protected Dict symbols_;

    this(string end_id) {
        symbols_ = new Dict(end_id);
        size_t delegate(in string, size_t) scanner = (in string s, size_t pos) { return pos; };
        symbols_[end_id].setScanner(scanner);
    }

    auto addSymbolToDict(string sym, int lbp=0) { 
        if (sym in symbols_) {
            auto s = symbols_[sym];
            s.lbp = max(s.lbp, lbp);
            return s;
        } else {
            auto s = new Sym(sym, lbp);
            return symbols_[sym] = s;
        }
    }

    struct Prefix {
        alias T delegate(T) handler_type;
        alias T delegate(Parser) func_type;
    }

    struct Postfix {
        alias T delegate(T) handler_type;
        alias T delegate(Parser, T) func_type;
    }

    struct LeftAssociative {
        alias T delegate(T, T) handler_type;
        alias T delegate(Parser, T) func_type;
    }

    struct RightAssociative {
        alias T delegate(T, T) handler_type;
        alias T delegate(Parser, T) func_type;
    }

    private static void setBehaviour(Semantics)(Sym sym, int rbp, Semantics.handler_type func) {
        static if (is (Semantics == Prefix)) 
            sym.nud = (Parser p) { return func(p.parse(rbp)); };
        else static if (is (Semantics == Postfix))
            sym.led = (Parser p, T left) { return func(left); };
        else static if (is (Semantics == LeftAssociative))
            sym.led = (Parser p, T left) { return func(left, p.parse(rbp)); };
        else static if (is (Semantics == RightAssociative))
            sym.led = (Parser p, T left) { return func(left, p.parse(rbp - 1)); };
        else 
            static assert(0, "expected Prefix | Postfix | LeftAssociative | RightAssociative");
    }

    private static void setBehaviour(Semantics)(Sym sym, Semantics.handler_type func) {
        setBehaviour!(Semantics)(sym, sym.lbp, func);
    }

    auto prefix(string op, int binding_power, T delegate(T) handler) {
        auto sym = addSymbolToDict(op, binding_power);
        setBehaviour!(Prefix)(sym, binding_power, handler);
        return sym;
    }

    auto prefix(string op, int binding_power, T delegate(T) handler, keep_symbol_lbp_t) {
        auto sym = addSymbolToDict(op);
        setBehaviour!(Prefix)(sym, binding_power, handler);
        return sym;
    }

    private template make_fix(string name, string semantics) {
        mixin("auto " ~ name ~ "(string op, int binding_power, " ~ semantics ~ ".handler_type handler) {
            auto sym = addSymbolToDict(op, binding_power);
            setBehaviour!(" ~ semantics ~ ")(sym, binding_power, handler);
            return sym;
        }");
    }

    mixin make_fix!("prefix", "Prefix");
    mixin make_fix!("postfix", "Postfix");
    mixin make_fix!("infix", "LeftAssociative");
    mixin make_fix!("infix_r", "RightAssociative");

    auto brackets(string ob, string cb, int binding_power, T delegate(T) handler=null) {
        auto open_sym = addSymbolToDict(ob, binding_power);
        addSymbolToDict(cb, 0);
        if (handler is null) {
            handler = (T val) { return val; };
        }
        open_sym.nud = (Parser p) { 
            T val = p.parse(0); 
            p.advance(cb); 
            return handler(val); 
        };
        return open_sym;
    }

    T parse(string text) {
        auto p = new Parser(text, symbols);
        return p.parse();
    }

    @property Dict symbols() { return symbols_; }
}

unittest {

    class Calculator : Grammar!int {
        this() {
            super("(end)");
            addSymbolToDict("(number)", 0).setScanner(
                (in string str, size_t pos) { 
                    size_t i = pos;
                    while (i < str.length && std.ascii.isDigit(str[i]))
                        ++i;
                    return i;
                }).setParser(
                (in string str) { 
                    return to!int(str);
                });
            infix("+", 10, (int a, int b){ return a + b;});
            infix("-", 10, (int a, int b){ return a - b;});
            infix("*", 20, (int a, int b){ return a * b;});
            infix("/", 20, (int a, int b){ return a / b;});
            prefix("+", 100, (int a){ return a; }, keep_symbol_lbp);
            prefix("-", 100, (int a){ return -a;}, keep_symbol_lbp);
            brackets("(", ")", int.max);
        }
    }

    auto calc = new Calculator();
    assert(calc.parse("2 + 3 * (4 - 5)") == -1);
    assert(calc.parse("2 * (3 + 2) / 5") == 2);
    assert(calc.parse("36 / (-2 + 2 * (2 + 3) * (1 + 1))") == 2);
}
