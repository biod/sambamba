// Written in the D programming language.

/**
   Copyright: Copyright Digital Mars 2000-2013.

   License: $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).

   Authors: $(HTTP walterbright.com, Walter Bright), $(HTTP erdani.com,
   Andrei Alexandrescu), and Kenji Hara

   Source: $(PHOBOSSRC std/_format.d)
 */
module contrib.undead.doformat;

//debug=format;                // uncomment to turn on debugging printf's

import core.vararg;
import std.exception;
import std.meta;
import std.range.primitives;
import std.traits;
import std.format;

version(CRuntime_DigitalMars)
{
    version = DigitalMarsC;
}

version (DigitalMarsC)
{
    // This is DMC's internal floating point formatting function
    extern (C)
    {
        extern shared char* function(int c, int flags, int precision,
                in real* pdval,
                char* buf, size_t* psl, int width) __pfloatfmt;
    }
}

/**********************************************************************
 * Signals a mismatch between a format and its corresponding argument.
 */
class FormatException : Exception
{
    @safe pure nothrow
    this()
    {
        super("format error");
    }

    @safe pure nothrow
    this(string msg, string fn = __FILE__, size_t ln = __LINE__, Throwable next = null)
    {
        super(msg, fn, ln, next);
    }
}


// Legacy implementation

enum Mangle : char
{
    Tvoid     = 'v',
    Tbool     = 'b',
    Tbyte     = 'g',
    Tubyte    = 'h',
    Tshort    = 's',
    Tushort   = 't',
    Tint      = 'i',
    Tuint     = 'k',
    Tlong     = 'l',
    Tulong    = 'm',
    Tfloat    = 'f',
    Tdouble   = 'd',
    Treal     = 'e',

    Tifloat   = 'o',
    Tidouble  = 'p',
    Tireal    = 'j',
    Tcfloat   = 'q',
    Tcdouble  = 'r',
    Tcreal    = 'c',

    Tchar     = 'a',
    Twchar    = 'u',
    Tdchar    = 'w',

    Tarray    = 'A',
    Tsarray   = 'G',
    Taarray   = 'H',
    Tpointer  = 'P',
    Tfunction = 'F',
    Tident    = 'I',
    Tclass    = 'C',
    Tstruct   = 'S',
    Tenum     = 'E',
    Ttypedef  = 'T',
    Tdelegate = 'D',

    Tconst    = 'x',
    Timmutable = 'y',
}

// return the TypeInfo for a primitive type and null otherwise.  This
// is required since for arrays of ints we only have the mangled char
// to work from. If arrays always subclassed TypeInfo_Array this
// routine could go away.
private TypeInfo primitiveTypeInfo(Mangle m)
{
    // BUG: should fix this in static this() to avoid double checked locking bug
    __gshared TypeInfo[Mangle] dic;
    if (!dic.length)
    {
        dic = [
            Mangle.Tvoid : typeid(void),
            Mangle.Tbool : typeid(bool),
            Mangle.Tbyte : typeid(byte),
            Mangle.Tubyte : typeid(ubyte),
            Mangle.Tshort : typeid(short),
            Mangle.Tushort : typeid(ushort),
            Mangle.Tint : typeid(int),
            Mangle.Tuint : typeid(uint),
            Mangle.Tlong : typeid(long),
            Mangle.Tulong : typeid(ulong),
            Mangle.Tfloat : typeid(float),
            Mangle.Tdouble : typeid(double),
            Mangle.Treal : typeid(real),
            Mangle.Tifloat : typeid(ifloat),
            Mangle.Tidouble : typeid(idouble),
            Mangle.Tireal : typeid(ireal),
            Mangle.Tcfloat : typeid(cfloat),
            Mangle.Tcdouble : typeid(cdouble),
            Mangle.Tcreal : typeid(creal),
            Mangle.Tchar : typeid(char),
            Mangle.Twchar : typeid(wchar),
            Mangle.Tdchar : typeid(dchar)
            ];
    }
    auto p = m in dic;
    return p ? *p : null;
}

// This stuff has been removed from the docs and is planned for deprecation.
/*
 * Interprets variadic argument list pointed to by argptr whose types
 * are given by arguments[], formats them according to embedded format
 * strings in the variadic argument list, and sends the resulting
 * characters to putc.
 *
 * The variadic arguments are consumed in order.  Each is formatted
 * into a sequence of chars, using the default format specification
 * for its type, and the characters are sequentially passed to putc.
 * If a $(D char[]), $(D wchar[]), or $(D dchar[]) argument is
 * encountered, it is interpreted as a format string. As many
 * arguments as specified in the format string are consumed and
 * formatted according to the format specifications in that string and
 * passed to putc. If there are too few remaining arguments, a
 * $(D FormatException) is thrown. If there are more remaining arguments than
 * needed by the format specification, the default processing of
 * arguments resumes until they are all consumed.
 *
 * Params:
 *        putc =        Output is sent do this delegate, character by character.
 *        arguments = Array of $(D TypeInfo)s, one for each argument to be formatted.
 *        argptr = Points to variadic argument list.
 *
 * Throws:
 *        Mismatched arguments and formats result in a $(D FormatException) being thrown.
 *
 * Format_String:
 *        <a name="format-string">$(I Format strings)</a>
 *        consist of characters interspersed with
 *        $(I format specifications). Characters are simply copied
 *        to the output (such as putc) after any necessary conversion
 *        to the corresponding UTF-8 sequence.
 *
 *        A $(I format specification) starts with a '%' character,
 *        and has the following grammar:

$(CONSOLE
$(I FormatSpecification):
    $(B '%%')
    $(B '%') $(I Flags) $(I Width) $(I Precision) $(I FormatChar)

$(I Flags):
    $(I empty)
    $(B '-') $(I Flags)
    $(B '+') $(I Flags)
    $(B '#') $(I Flags)
    $(B '0') $(I Flags)
    $(B ' ') $(I Flags)

$(I Width):
    $(I empty)
    $(I Integer)
    $(B '*')

$(I Precision):
    $(I empty)
    $(B '.')
    $(B '.') $(I Integer)
    $(B '.*')

$(I Integer):
    $(I Digit)
    $(I Digit) $(I Integer)

$(I Digit):
    $(B '0')
    $(B '1')
    $(B '2')
    $(B '3')
    $(B '4')
    $(B '5')
    $(B '6')
    $(B '7')
    $(B '8')
    $(B '9')

$(I FormatChar):
    $(B 's')
    $(B 'b')
    $(B 'd')
    $(B 'o')
    $(B 'x')
    $(B 'X')
    $(B 'e')
    $(B 'E')
    $(B 'f')
    $(B 'F')
    $(B 'g')
    $(B 'G')
    $(B 'a')
    $(B 'A')
)
    $(DL
    $(DT $(I Flags))
    $(DL
        $(DT $(B '-'))
        $(DD
        Left justify the result in the field.
        It overrides any $(B 0) flag.)

        $(DT $(B '+'))
        $(DD Prefix positive numbers in a signed conversion with a $(B +).
        It overrides any $(I space) flag.)

        $(DT $(B '#'))
        $(DD Use alternative formatting:
        $(DL
            $(DT For $(B 'o'):)
            $(DD Add to precision as necessary so that the first digit
            of the octal formatting is a '0', even if both the argument
            and the $(I Precision) are zero.)
            $(DT For $(B 'x') ($(B 'X')):)
            $(DD If non-zero, prefix result with $(B 0x) ($(B 0X)).)
            $(DT For floating point formatting:)
            $(DD Always insert the decimal point.)
            $(DT For $(B 'g') ($(B 'G')):)
            $(DD Do not elide trailing zeros.)
        ))

        $(DT $(B '0'))
        $(DD For integer and floating point formatting when not nan or
        infinity, use leading zeros
        to pad rather than spaces.
        Ignore if there's a $(I Precision).)

        $(DT $(B ' '))
        $(DD Prefix positive numbers in a signed conversion with a space.)
    )

    $(DT $(I Width))
    $(DD
    Specifies the minimum field width.
    If the width is a $(B *), the next argument, which must be
    of type $(B int), is taken as the width.
    If the width is negative, it is as if the $(B -) was given
    as a $(I Flags) character.)

    $(DT $(I Precision))
    $(DD Gives the precision for numeric conversions.
    If the precision is a $(B *), the next argument, which must be
    of type $(B int), is taken as the precision. If it is negative,
    it is as if there was no $(I Precision).)

    $(DT $(I FormatChar))
    $(DD
    $(DL
        $(DT $(B 's'))
        $(DD The corresponding argument is formatted in a manner consistent
        with its type:
        $(DL
            $(DT $(B bool))
            $(DD The result is <tt>'true'</tt> or <tt>'false'</tt>.)
            $(DT integral types)
            $(DD The $(B %d) format is used.)
            $(DT floating point types)
            $(DD The $(B %g) format is used.)
            $(DT string types)
            $(DD The result is the string converted to UTF-8.)
            A $(I Precision) specifies the maximum number of characters
            to use in the result.
            $(DT classes derived from $(B Object))
            $(DD The result is the string returned from the class instance's
            $(B .toString()) method.
            A $(I Precision) specifies the maximum number of characters
            to use in the result.)
            $(DT non-string static and dynamic arrays)
            $(DD The result is [s<sub>0</sub>, s<sub>1</sub>, ...]
            where s<sub>k</sub> is the kth element
            formatted with the default format.)
        ))

        $(DT $(B 'b','d','o','x','X'))
        $(DD The corresponding argument must be an integral type
        and is formatted as an integer. If the argument is a signed type
        and the $(I FormatChar) is $(B d) it is converted to
        a signed string of characters, otherwise it is treated as
        unsigned. An argument of type $(B bool) is formatted as '1'
        or '0'. The base used is binary for $(B b), octal for $(B o),
        decimal
        for $(B d), and hexadecimal for $(B x) or $(B X).
        $(B x) formats using lower case letters, $(B X) uppercase.
        If there are fewer resulting digits than the $(I Precision),
        leading zeros are used as necessary.
        If the $(I Precision) is 0 and the number is 0, no digits
        result.)

        $(DT $(B 'e','E'))
        $(DD A floating point number is formatted as one digit before
        the decimal point, $(I Precision) digits after, the $(I FormatChar),
        &plusmn;, followed by at least a two digit exponent: $(I d.dddddd)e$(I &plusmn;dd).
        If there is no $(I Precision), six
        digits are generated after the decimal point.
        If the $(I Precision) is 0, no decimal point is generated.)

        $(DT $(B 'f','F'))
        $(DD A floating point number is formatted in decimal notation.
        The $(I Precision) specifies the number of digits generated
        after the decimal point. It defaults to six. At least one digit
        is generated before the decimal point. If the $(I Precision)
        is zero, no decimal point is generated.)

        $(DT $(B 'g','G'))
        $(DD A floating point number is formatted in either $(B e) or
        $(B f) format for $(B g); $(B E) or $(B F) format for
        $(B G).
        The $(B f) format is used if the exponent for an $(B e) format
        is greater than -5 and less than the $(I Precision).
        The $(I Precision) specifies the number of significant
        digits, and defaults to six.
        Trailing zeros are elided after the decimal point, if the fractional
        part is zero then no decimal point is generated.)

        $(DT $(B 'a','A'))
        $(DD A floating point number is formatted in hexadecimal
        exponential notation 0x$(I h.hhhhhh)p$(I &plusmn;d).
        There is one hexadecimal digit before the decimal point, and as
        many after as specified by the $(I Precision).
        If the $(I Precision) is zero, no decimal point is generated.
        If there is no $(I Precision), as many hexadecimal digits as
        necessary to exactly represent the mantissa are generated.
        The exponent is written in as few digits as possible,
        but at least one, is in decimal, and represents a power of 2 as in
        $(I h.hhhhhh)*2<sup>$(I &plusmn;d)</sup>.
        The exponent for zero is zero.
        The hexadecimal digits, x and p are in upper case if the
        $(I FormatChar) is upper case.)
    )

    Floating point NaN's are formatted as $(B nan) if the
    $(I FormatChar) is lower case, or $(B NAN) if upper.
    Floating point infinities are formatted as $(B inf) or
    $(B infinity) if the
    $(I FormatChar) is lower case, or $(B INF) or $(B INFINITY) if upper.
    ))

Example:

-------------------------
import core.stdc.stdio;
import std.format;

void myPrint(...)
{
    void putc(dchar c)
    {
        fputc(c, stdout);
    }

    std.format.doFormat(&putc, _arguments, _argptr);
}

void main()
{
    int x = 27;

    // prints 'The answer is 27:6'
    myPrint("The answer is %s:", x, 6);
}
------------------------
 */
void doFormat()(scope void delegate(dchar) putc, TypeInfo[] arguments, va_list ap)
{
    import std.utf : encode, toUCSindex, isValidDchar, UTFException, toUTF8;
    import core.stdc.string : strlen;
    import core.stdc.stdlib : alloca, malloc, realloc, free;
    import core.stdc.stdio : snprintf;

    size_t bufLength = 1024;
    void* argBuffer = malloc(bufLength);
    scope(exit) free(argBuffer);

    size_t bufUsed = 0;
    foreach (ti; arguments)
    {
        // Ensure the required alignment
        bufUsed += ti.talign - 1;
        bufUsed -= (cast(size_t)argBuffer + bufUsed) & (ti.talign - 1);
        auto pos = bufUsed;
        // Align to next word boundary
        bufUsed += ti.tsize + size_t.sizeof - 1;
        bufUsed -= (cast(size_t)argBuffer + bufUsed) & (size_t.sizeof - 1);
        // Resize buffer if necessary
        while (bufUsed > bufLength)
        {
            bufLength *= 2;
            argBuffer = realloc(argBuffer, bufLength);
        }
        // Copy argument into buffer
        va_arg(ap, ti, argBuffer + pos);
    }

    auto argptr = argBuffer;
    void* skipArg(TypeInfo ti)
    {
        // Ensure the required alignment
        argptr += ti.talign - 1;
        argptr -= cast(size_t)argptr & (ti.talign - 1);
        auto p = argptr;
        // Align to next word boundary
        argptr += ti.tsize + size_t.sizeof - 1;
        argptr -= cast(size_t)argptr & (size_t.sizeof - 1);
        return p;
    }
    auto getArg(T)()
    {
        return *cast(T*)skipArg(typeid(T));
    }

    TypeInfo ti;
    Mangle m;
    uint flags;
    int field_width;
    int precision;

    enum : uint
    {
        FLdash = 1,
        FLplus = 2,
        FLspace = 4,
        FLhash = 8,
        FLlngdbl = 0x20,
        FL0pad = 0x40,
        FLprecision = 0x80,
    }

    static TypeInfo skipCI(TypeInfo valti)
    {
        for (;;)
        {
            if (typeid(valti).name.length == 18 &&
                    typeid(valti).name[9..18] == "Invariant")
                valti = (cast(TypeInfo_Invariant)valti).base;
            else if (typeid(valti).name.length == 14 &&
                    typeid(valti).name[9..14] == "Const")
                valti = (cast(TypeInfo_Const)valti).base;
            else
                break;
        }

        return valti;
    }

    void formatArg(char fc)
    {
        bool vbit;
        ulong vnumber;
        char vchar;
        dchar vdchar;
        Object vobject;
        real vreal;
        creal vcreal;
        Mangle m2;
        int signed = 0;
        uint base = 10;
        int uc;
        char[ulong.sizeof * 8] tmpbuf; // long enough to print long in binary
        const(char)* prefix = "";
        string s;

        void putstr(const char[] s)
        {
            //printf("putstr: s = %.*s, flags = x%x\n", s.length, s.ptr, flags);
            ptrdiff_t padding = field_width -
                (strlen(prefix) + toUCSindex(s, s.length));
            ptrdiff_t prepad = 0;
            ptrdiff_t postpad = 0;
            if (padding > 0)
            {
                if (flags & FLdash)
                    postpad = padding;
                else
                    prepad = padding;
            }

            if (flags & FL0pad)
            {
                while (*prefix)
                    putc(*prefix++);
                while (prepad--)
                    putc('0');
            }
            else
            {
                while (prepad--)
                    putc(' ');
                while (*prefix)
                    putc(*prefix++);
            }

            foreach (dchar c; s)
                putc(c);

            while (postpad--)
                putc(' ');
        }

        void putreal(real v)
        {
            //printf("putreal %Lg\n", vreal);

            switch (fc)
            {
                case 's':
                    fc = 'g';
                    break;

                case 'f', 'F', 'e', 'E', 'g', 'G', 'a', 'A':
                    break;

                default:
                    //printf("fc = '%c'\n", fc);
                Lerror:
                    throw new FormatException("incompatible format character for floating point type");
            }
            version (DigitalMarsC)
            {
                uint sl;
                char[] fbuf = tmpbuf;
                if (!(flags & FLprecision))
                    precision = 6;
                while (1)
                {
                    sl = fbuf.length;
                    prefix = (*__pfloatfmt)(fc, flags | FLlngdbl,
                            precision, &v, cast(char*)fbuf, &sl, field_width);
                    if (sl != -1)
                        break;
                    sl = fbuf.length * 2;
                    fbuf = (cast(char*)alloca(sl * char.sizeof))[0 .. sl];
                }
                putstr(fbuf[0 .. sl]);
            }
            else
            {
                ptrdiff_t sl;
                char[] fbuf = tmpbuf;
                char[12] format;
                format[0] = '%';
                int i = 1;
                if (flags & FLdash)
                    format[i++] = '-';
                if (flags & FLplus)
                    format[i++] = '+';
                if (flags & FLspace)
                    format[i++] = ' ';
                if (flags & FLhash)
                    format[i++] = '#';
                if (flags & FL0pad)
                    format[i++] = '0';
                format[i + 0] = '*';
                format[i + 1] = '.';
                format[i + 2] = '*';
                format[i + 3] = 'L';
                format[i + 4] = fc;
                format[i + 5] = 0;
                if (!(flags & FLprecision))
                    precision = -1;
                while (1)
                {
                    sl = fbuf.length;
                    int n;
                    version (CRuntime_Microsoft)
                    {
                        import std.math : isNaN, isInfinity;
                        if (isNaN(v)) // snprintf writes 1.#QNAN
                            n = snprintf(fbuf.ptr, sl, "nan");
                        else if (isInfinity(v)) // snprintf writes 1.#INF
                            n = snprintf(fbuf.ptr, sl, v < 0 ? "-inf" : "inf");
                        else
                            n = snprintf(fbuf.ptr, sl, format.ptr, field_width,
                                         precision, cast(double)v);
                    }
                    else
                        n = snprintf(fbuf.ptr, sl, format.ptr, field_width,
                                precision, v);
                    //printf("format = '%s', n = %d\n", cast(char*)format, n);
                    if (n >= 0 && n < sl)
                    {        sl = n;
                        break;
                    }
                    if (n < 0)
                        sl = sl * 2;
                    else
                        sl = n + 1;
                    fbuf = (cast(char*)alloca(sl * char.sizeof))[0 .. sl];
                }
                putstr(fbuf[0 .. sl]);
            }
            return;
        }

        static Mangle getMan(TypeInfo ti)
        {
          auto m = cast(Mangle)typeid(ti).name[9];
          if (typeid(ti).name.length == 20 &&
              typeid(ti).name[9..20] == "StaticArray")
                m = cast(Mangle)'G';
          return m;
        }

        /* p = pointer to the first element in the array
         * len = number of elements in the array
         * valti = type of the elements
         */
        void putArray(void* p, size_t len, TypeInfo valti)
        {
          //printf("\nputArray(len = %u), tsize = %u\n", len, valti.tsize);
          putc('[');
          valti = skipCI(valti);
          size_t tsize = valti.tsize;
          auto argptrSave = argptr;
          auto tiSave = ti;
          auto mSave = m;
          ti = valti;
          //printf("\n%.*s\n", typeid(valti).name.length, typeid(valti).name.ptr);
          m = getMan(valti);
          while (len--)
          {
            //doFormat(putc, (&valti)[0 .. 1], p);
            argptr = p;
            formatArg('s');
            p += tsize;
            if (len > 0) putc(',');
          }
          m = mSave;
          ti = tiSave;
          argptr = argptrSave;
          putc(']');
        }

        void putAArray(ubyte[long] vaa, TypeInfo valti, TypeInfo keyti)
        {
            putc('[');
            bool comma=false;
            auto argptrSave = argptr;
            auto tiSave = ti;
            auto mSave = m;
            valti = skipCI(valti);
            keyti = skipCI(keyti);
            foreach (ref fakevalue; vaa)
            {
                if (comma) putc(',');
                comma = true;
                void *pkey = &fakevalue;
                version (D_LP64)
                    pkey -= (long.sizeof + 15) & ~(15);
                else
                    pkey -= (long.sizeof + size_t.sizeof - 1) & ~(size_t.sizeof - 1);

                // the key comes before the value
                auto keysize = keyti.tsize;
                version (D_LP64)
                    auto keysizet = (keysize + 15) & ~(15);
                else
                    auto keysizet = (keysize + size_t.sizeof - 1) & ~(size_t.sizeof - 1);

                void* pvalue = pkey + keysizet;

                //doFormat(putc, (&keyti)[0..1], pkey);
                m = getMan(keyti);
                argptr = pkey;

                ti = keyti;
                formatArg('s');

                putc(':');
                //doFormat(putc, (&valti)[0..1], pvalue);
                m = getMan(valti);
                argptr = pvalue;

                ti = valti;
                formatArg('s');
            }
            m = mSave;
            ti = tiSave;
            argptr = argptrSave;
            putc(']');
        }

        //printf("formatArg(fc = '%c', m = '%c')\n", fc, m);
        int mi;
        switch (m)
        {
            case Mangle.Tbool:
                vbit = getArg!(bool)();
                if (fc != 's')
                {   vnumber = vbit;
                    goto Lnumber;
                }
                putstr(vbit ? "true" : "false");
                return;

            case Mangle.Tchar:
                vchar = getArg!(char)();
                if (fc != 's')
                {   vnumber = vchar;
                    goto Lnumber;
                }
            L2:
                putstr((&vchar)[0 .. 1]);
                return;

            case Mangle.Twchar:
                vdchar = getArg!(wchar)();
                goto L1;

            case Mangle.Tdchar:
                vdchar = getArg!(dchar)();
            L1:
                if (fc != 's')
                {   vnumber = vdchar;
                    goto Lnumber;
                }
                if (vdchar <= 0x7F)
                {   vchar = cast(char)vdchar;
                    goto L2;
                }
                else
                {   if (!isValidDchar(vdchar))
                        throw new UTFException("invalid dchar in format");
                    char[4] vbuf;
                    putstr(vbuf[0 .. encode(vbuf, vdchar)]);
                }
                return;

            case Mangle.Tbyte:
                signed = 1;
                vnumber = getArg!(byte)();
                goto Lnumber;

            case Mangle.Tubyte:
                vnumber = getArg!(ubyte)();
                goto Lnumber;

            case Mangle.Tshort:
                signed = 1;
                vnumber = getArg!(short)();
                goto Lnumber;

            case Mangle.Tushort:
                vnumber = getArg!(ushort)();
                goto Lnumber;

            case Mangle.Tint:
                signed = 1;
                vnumber = getArg!(int)();
                goto Lnumber;

            case Mangle.Tuint:
            Luint:
                vnumber = getArg!(uint)();
                goto Lnumber;

            case Mangle.Tlong:
                signed = 1;
                vnumber = cast(ulong)getArg!(long)();
                goto Lnumber;

            case Mangle.Tulong:
            Lulong:
                vnumber = getArg!(ulong)();
                goto Lnumber;

            case Mangle.Tclass:
                vobject = getArg!(Object)();
                if (vobject is null)
                    s = "null";
                else
                    s = vobject.toString();
                goto Lputstr;

            case Mangle.Tpointer:
                vnumber = cast(ulong)getArg!(void*)();
                if (fc != 'x')  uc = 1;
                flags |= FL0pad;
                if (!(flags & FLprecision))
                {   flags |= FLprecision;
                    precision = (void*).sizeof;
                }
                base = 16;
                goto Lnumber;

            case Mangle.Tfloat:
            case Mangle.Tifloat:
                if (fc == 'x' || fc == 'X')
                    goto Luint;
                vreal = getArg!(float)();
                goto Lreal;

            case Mangle.Tdouble:
            case Mangle.Tidouble:
                if (fc == 'x' || fc == 'X')
                    goto Lulong;
                vreal = getArg!(double)();
                goto Lreal;

            case Mangle.Treal:
            case Mangle.Tireal:
                vreal = getArg!(real)();
                goto Lreal;

            case Mangle.Tcfloat:
                vcreal = getArg!(cfloat)();
                goto Lcomplex;

            case Mangle.Tcdouble:
                vcreal = getArg!(cdouble)();
                goto Lcomplex;

            case Mangle.Tcreal:
                vcreal = getArg!(creal)();
                goto Lcomplex;

            case Mangle.Tsarray:
                putArray(argptr, (cast(TypeInfo_StaticArray)ti).len, (cast(TypeInfo_StaticArray)ti).next);
                return;

            case Mangle.Tarray:
                mi = 10;
                if (typeid(ti).name.length == 14 &&
                    typeid(ti).name[9..14] == "Array")
                { // array of non-primitive types
                  TypeInfo tn = (cast(TypeInfo_Array)ti).next;
                  tn = skipCI(tn);
                  switch (cast(Mangle)typeid(tn).name[9])
                  {
                    case Mangle.Tchar:  goto LarrayChar;
                    case Mangle.Twchar: goto LarrayWchar;
                    case Mangle.Tdchar: goto LarrayDchar;
                    default:
                        break;
                  }
                  void[] va = getArg!(void[])();
                  putArray(va.ptr, va.length, tn);
                  return;
                }
                if (typeid(ti).name.length == 25 &&
                    typeid(ti).name[9..25] == "AssociativeArray")
                { // associative array
                  ubyte[long] vaa = getArg!(ubyte[long])();
                  putAArray(vaa,
                        (cast(TypeInfo_AssociativeArray)ti).next,
                        (cast(TypeInfo_AssociativeArray)ti).key);
                  return;
                }

                while (1)
                {
                    m2 = cast(Mangle)typeid(ti).name[mi];
                    switch (m2)
                    {
                        case Mangle.Tchar:
                        LarrayChar:
                            s = getArg!(string)();
                            goto Lputstr;

                        case Mangle.Twchar:
                        LarrayWchar:
                            wchar[] sw = getArg!(wchar[])();
                            s = toUTF8(sw);
                            goto Lputstr;

                        case Mangle.Tdchar:
                        LarrayDchar:
                            s = toUTF8(getArg!(dstring)());
                        Lputstr:
                            if (fc != 's')
                                throw new FormatException("string");
                            if (flags & FLprecision && precision < s.length)
                                s = s[0 .. precision];
                            putstr(s);
                            break;

                        case Mangle.Tconst:
                        case Mangle.Timmutable:
                            mi++;
                            continue;

                        default:
                            TypeInfo ti2 = primitiveTypeInfo(m2);
                            if (!ti2)
                              goto Lerror;
                            void[] va = getArg!(void[])();
                            putArray(va.ptr, va.length, ti2);
                    }
                    return;
                }
                assert(0);

            case Mangle.Tenum:
                ti = (cast(TypeInfo_Enum)ti).base;
                m = cast(Mangle)typeid(ti).name[9];
                formatArg(fc);
                return;

            case Mangle.Tstruct:
            {        TypeInfo_Struct tis = cast(TypeInfo_Struct)ti;
                if (tis.xtoString is null)
                    throw new FormatException("Can't convert " ~ tis.toString()
                            ~ " to string: \"string toString()\" not defined");
                s = tis.xtoString(skipArg(tis));
                goto Lputstr;
            }

            default:
                goto Lerror;
        }

    Lnumber:
        switch (fc)
        {
            case 's':
            case 'd':
                if (signed)
                {   if (cast(long)vnumber < 0)
                    {        prefix = "-";
                        vnumber = -vnumber;
                    }
                    else if (flags & FLplus)
                        prefix = "+";
                    else if (flags & FLspace)
                        prefix = " ";
                }
                break;

            case 'b':
                signed = 0;
                base = 2;
                break;

            case 'o':
                signed = 0;
                base = 8;
                break;

            case 'X':
                uc = 1;
                if (flags & FLhash && vnumber)
                    prefix = "0X";
                signed = 0;
                base = 16;
                break;

            case 'x':
                if (flags & FLhash && vnumber)
                    prefix = "0x";
                signed = 0;
                base = 16;
                break;

            default:
                goto Lerror;
        }

        if (!signed)
        {
            switch (m)
            {
                case Mangle.Tbyte:
                    vnumber &= 0xFF;
                    break;

                case Mangle.Tshort:
                    vnumber &= 0xFFFF;
                    break;

                case Mangle.Tint:
                    vnumber &= 0xFFFFFFFF;
                    break;

                default:
                    break;
            }
        }

        if (flags & FLprecision && fc != 'p')
            flags &= ~FL0pad;

        if (vnumber < base)
        {
            if (vnumber == 0 && precision == 0 && flags & FLprecision &&
                !(fc == 'o' && flags & FLhash))
            {
                putstr(null);
                return;
            }
            if (precision == 0 || !(flags & FLprecision))
            {        vchar = cast(char)('0' + vnumber);
                if (vnumber < 10)
                    vchar = cast(char)('0' + vnumber);
                else
                    vchar = cast(char)((uc ? 'A' - 10 : 'a' - 10) + vnumber);
                goto L2;
            }
        }

        {
            ptrdiff_t n = tmpbuf.length;
            char c;
            int hexoffset = uc ? ('A' - ('9' + 1)) : ('a' - ('9' + 1));

            while (vnumber)
            {
                c = cast(char)((vnumber % base) + '0');
                if (c > '9')
                    c += hexoffset;
                vnumber /= base;
                tmpbuf[--n] = c;
            }
            if (tmpbuf.length - n < precision && precision < tmpbuf.length)
            {
                ptrdiff_t m = tmpbuf.length - precision;
                tmpbuf[m .. n] = '0';
                n = m;
            }
            else if (flags & FLhash && fc == 'o')
                prefix = "0";
            putstr(tmpbuf[n .. tmpbuf.length]);
            return;
        }

    Lreal:
        putreal(vreal);
        return;

    Lcomplex:
        putreal(vcreal.re);
        if (vcreal.im >= 0)
        {
            putc('+');
        }
        putreal(vcreal.im);
        putc('i');
        return;

    Lerror:
        throw new FormatException("formatArg");
    }

    for (int j = 0; j < arguments.length; )
    {
        ti = arguments[j++];
        //printf("arg[%d]: '%.*s' %d\n", j, typeid(ti).name.length, typeid(ti).name.ptr, typeid(ti).name.length);
        //ti.print();

        flags = 0;
        precision = 0;
        field_width = 0;

        ti = skipCI(ti);
        int mi = 9;
        do
        {
            if (typeid(ti).name.length <= mi)
                goto Lerror;
            m = cast(Mangle)typeid(ti).name[mi++];
        } while (m == Mangle.Tconst || m == Mangle.Timmutable);

        if (m == Mangle.Tarray)
        {
            if (typeid(ti).name.length == 14 &&
                    typeid(ti).name[9..14] == "Array")
            {
                TypeInfo tn = (cast(TypeInfo_Array)ti).next;
                tn = skipCI(tn);
                switch (cast(Mangle)typeid(tn).name[9])
                {
                case Mangle.Tchar:
                case Mangle.Twchar:
                case Mangle.Tdchar:
                    ti = tn;
                    mi = 9;
                    break;
                default:
                    break;
                }
            }
          L1:
            Mangle m2 = cast(Mangle)typeid(ti).name[mi];
            string  fmt;                        // format string
            wstring wfmt;
            dstring dfmt;

            /* For performance reasons, this code takes advantage of the
             * fact that most format strings will be ASCII, and that the
             * format specifiers are always ASCII. This means we only need
             * to deal with UTF in a couple of isolated spots.
             */

            switch (m2)
            {
            case Mangle.Tchar:
                fmt = getArg!(string)();
                break;

            case Mangle.Twchar:
                wfmt = getArg!(wstring)();
                fmt = toUTF8(wfmt);
                break;

            case Mangle.Tdchar:
                dfmt = getArg!(dstring)();
                fmt = toUTF8(dfmt);
                break;

            case Mangle.Tconst:
            case Mangle.Timmutable:
                mi++;
                goto L1;

            default:
                formatArg('s');
                continue;
            }

            for (size_t i = 0; i < fmt.length; )
            {        dchar c = fmt[i++];

                dchar getFmtChar()
                {   // Valid format specifier characters will never be UTF
                    if (i == fmt.length)
                        throw new FormatException("invalid specifier");
                    return fmt[i++];
                }

                int getFmtInt()
                {   int n;

                    while (1)
                    {
                        n = n * 10 + (c - '0');
                        if (n < 0)        // overflow
                            throw new FormatException("int overflow");
                        c = getFmtChar();
                        if (c < '0' || c > '9')
                            break;
                    }
                    return n;
                }

                int getFmtStar()
                {   Mangle m;
                    TypeInfo ti;

                    if (j == arguments.length)
                        throw new FormatException("too few arguments");
                    ti = arguments[j++];
                    m = cast(Mangle)typeid(ti).name[9];
                    if (m != Mangle.Tint)
                        throw new FormatException("int argument expected");
                    return getArg!(int)();
                }

                if (c != '%')
                {
                    if (c > 0x7F)        // if UTF sequence
                    {
                        i--;                // back up and decode UTF sequence
                        import std.utf : decode;
                        c = decode(fmt, i);
                    }
                  Lputc:
                    putc(c);
                    continue;
                }

                // Get flags {-+ #}
                flags = 0;
                while (1)
                {
                    c = getFmtChar();
                    switch (c)
                    {
                    case '-':        flags |= FLdash;        continue;
                    case '+':        flags |= FLplus;        continue;
                    case ' ':        flags |= FLspace;        continue;
                    case '#':        flags |= FLhash;        continue;
                    case '0':        flags |= FL0pad;        continue;

                    case '%':        if (flags == 0)
                                          goto Lputc;
                                     break;

                    default:         break;
                    }
                    break;
                }

                // Get field width
                field_width = 0;
                if (c == '*')
                {
                    field_width = getFmtStar();
                    if (field_width < 0)
                    {   flags |= FLdash;
                        field_width = -field_width;
                    }

                    c = getFmtChar();
                }
                else if (c >= '0' && c <= '9')
                    field_width = getFmtInt();

                if (flags & FLplus)
                    flags &= ~FLspace;
                if (flags & FLdash)
                    flags &= ~FL0pad;

                // Get precision
                precision = 0;
                if (c == '.')
                {   flags |= FLprecision;
                    //flags &= ~FL0pad;

                    c = getFmtChar();
                    if (c == '*')
                    {
                        precision = getFmtStar();
                        if (precision < 0)
                        {   precision = 0;
                            flags &= ~FLprecision;
                        }

                        c = getFmtChar();
                    }
                    else if (c >= '0' && c <= '9')
                        precision = getFmtInt();
                }

                if (j == arguments.length)
                    goto Lerror;
                ti = arguments[j++];
                ti = skipCI(ti);
                mi = 9;
                do
                {
                    m = cast(Mangle)typeid(ti).name[mi++];
                } while (m == Mangle.Tconst || m == Mangle.Timmutable);

                if (c > 0x7F)                // if UTF sequence
                    goto Lerror;        // format specifiers can't be UTF
                formatArg(cast(char)c);
            }
        }
        else
        {
            formatArg('s');
        }
    }
    return;

  Lerror:
    throw new FormatException();
}


private bool needToSwapEndianess(Char)(ref FormatSpec!Char f)
{
    import std.system : endian, Endian;

    return endian == Endian.littleEndian && f.flPlus
        || endian == Endian.bigEndian && f.flDash;
}

unittest
{   
    string res;
    void putc(dchar c)
    {
        res ~= c;
    }

    void myPrint(...)
    {
        doFormat(&putc, _arguments, _argptr);
    }

    myPrint("The answer is %s:", 27, 6);
    assert(res == "The answer is 27:6");
}
