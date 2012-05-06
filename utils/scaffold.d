import std.traits;
import std.algorithm;
import std.conv;

template convertArithmeticType(T) 
    if (__traits(isArithmetic, T))
{
    alias Unqual!T U;

    /* https://github.com/ffi/ffi/wiki/Types
     * http://dlang.org/abi.html
     */
       
    static if (is(U == ubyte)) {
        enum convertArithmeticType = "uchar";
    } else static if (is(U == long)) {
        enum convertArithmeticType = "int64";
    } else static if (is(U == ulong)) {
        enum convertArithmeticType = "uint64";
    } else static if (!is(U == real)) {
        /* char
         * bool
         * short
         * ushort
         * int
         * uint
         * float
         * double
         */
        enum convertArithmeticType = U.stringof;
    }
}

string toRuby(S)() 
    if (is(S == struct))
{
    string[] layout;
    string[] methods;

    foreach (member; __traits(allMembers, S)) {
        enum fullname = S.stringof ~ "." ~ member;
        alias Unqual!(typeof(mixin(fullname))) FieldType;

        static if (isCallable!(mixin(S.stringof ~ "." ~ member))) {
            /* skip methods and properties */
        } else static if (__traits(isArithmetic, FieldType)) {
            layout ~= ":" ~ member ~ ", :" ~ convertArithmeticType!FieldType;
            methods ~= "  def " ~ member ~ "\n" ~
                       "    self[:" ~ member ~ "]\n" ~
                       "  end";
        } else static if (is(FieldType == string)) {
            layout ~= ":" ~ member ~ "_len, :size_t";
            layout ~= ":" ~ member ~ "_data, :pointer";
            methods ~= "  def " ~ member ~ "\n" ~
                       "    return nil if self[:" ~ member ~ "_data].address == 0\n" ~
                       "    self[:" ~ member ~ "_data].read_string(self[:" ~ 
                                member ~ "_len])\n" ~
                       "  end";
        } else static if (isDynamicArray!FieldType) {

            layout ~= ":" ~ member ~ "_len, :size_t";
            layout ~= ":" ~ member ~ "_data, :pointer";
            alias Unqual!(typeof(FieldType.init[0])) ElementType;

            static if (__traits(isArithmetic, typeof(FieldType.init[0]))) {
                enum type = convertArithmeticType!ElementType;
                methods ~= "  def " ~ member ~ "\n" ~
                           "    return nil if self[:" ~ member ~ "_data].address == 0\n" ~
                           "    self[:" ~ member ~ "_data].read_array_of_" ~ type ~ 
                                    "(self[:" ~ member ~ "_len])\n" ~
                           "  end";
            } else {
                /* TODO: add some checks */
                /* the user is responsible for generating bindings for ElementType */
                /* it's assumed that there's a Ruby type with name ElementType.stringof */
                methods ~= "  def " ~ member ~ "\n" ~
                           "    return nil if self[:" ~ member ~ "_data].address == 0\n" ~
                           "    (0...self[:" ~ member ~ "_len]).map {|i|\n" ~
                           "      " ~ ElementType.stringof ~ ".new(" ~ 
                                            "self[:" ~ member ~ "_data] + i * " ~ 
                                            ElementType.stringof ~ ".size)\n" ~
                           "    }\n" ~
                           "  end";
            }
        } else {
            mixin("\"No support for type '" ~ FieldType.stringof ~ "' yet.\";");
        }
    }

    return "class " ~ S.stringof ~ " < FFI::Struct\n" ~
           "  layout " ~ to!string(joiner(layout, ",\n         ")) ~ "\n\n" ~
           to!string(joiner(methods, "\n\n")) ~ "\nend";
}
