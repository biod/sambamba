/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.std.hts.bam.validation.alignment;

public import bio.std.hts.bam.read;
public import bio.std.hts.bam.tagvalue;
import bio.core.utils.algo;

import std.algorithm;
import std.ascii;
import std.conv;
import std.typetuple;

/**
    Alignment validation error types.

    InvalidCigar error is accompanied by some CigarError,
    InvalidTags is accompanied by some TagError.
*/
enum AlignmentError {
    EmptyReadName, ///
    TooLongReadName, /// 
    ReadNameContainsInvalidCharacters, ///
    PositionIsOutOfRange, /// 
    QualityDataContainsInvalidElements, ///
    InvalidCigar, ///
    InvalidTags, ///
    DuplicateTagKeys ///
}

/// CIGAR string validation error types.
enum CigarError {
    InternalHardClipping, /// 
    InternalSoftClipping, ///
    InconsistentLength ///
}

/// Auxiliary data validation error types.
///
/// Refers to an individual tag.
enum TagError {
    EmptyString, ///
    EmptyHexadecimalString, ///
    NonPrintableString, ///
    NonPrintableCharacter, ///
    InvalidHexadecimalString, ///
    ExpectedIntegerValue, ///
    ExpectedStringValue, ///
    InvalidValueType, ///
    InvalidQualityString, ///
    ExpectedStringWithSameLengthAsSequence ///
}

/// Designates pair of predefined key from SAM/BAM specification
/// and expected type of tags with that key.
struct TagType(string key, T) {
    enum Key = key;
    alias T Type;
}

/// Compile-time available information about predefined tags
alias TypeTuple!(TagType!("AM", int),
                 TagType!("AS", int),
                 TagType!("BC", string),
                 TagType!("BQ", string),
                 TagType!("CC", string),
                 TagType!("CM", int),
                 TagType!("CP", int),
                 TagType!("CQ", string),
                 TagType!("CS", string),
                 TagType!("E2", string),
                 TagType!("FI", int),
                 TagType!("FS", string),
                 TagType!("FZ", ushort[]),
                 TagType!("LB", string),
                 TagType!("H0", int),
                 TagType!("H1", int),
                 TagType!("H2", int),
                 TagType!("HI", int),
                 TagType!("IH", int),
                 TagType!("MD", string),
                 TagType!("MQ", int),
                 TagType!("NH", int),
                 TagType!("NM", int),
                 TagType!("OQ", string),
                 TagType!("OP", int),
                 TagType!("OC", string),
                 TagType!("PG", string),
                 TagType!("PQ", int),
                 TagType!("PU", string),
                 TagType!("Q2", string),
                 TagType!("R2", string),
                 TagType!("RG", string),
                 TagType!("SM", int),
                 TagType!("TC", int),
                 TagType!("U2", string),
                 TagType!("UQ", int))
    PredefinedTags;


private template GetKey(U) {
    enum GetKey = U.Key;
}

private template PredefinedTagTypeHelper(string s) {
    alias PredefinedTags[staticIndexOf!(s, staticMap!(GetKey, PredefinedTags))] PredefinedTagTypeHelper;
}

/// Get predefined tag type by its key, in compile-time.
template PredefinedTagType(string s) {
    alias PredefinedTagTypeHelper!(s).Type PredefinedTagType;
}

/**
  Abstract class encapsulating visitation of SAM header elements.
*/
abstract class AbstractAlignmentValidator {
    /// Start validation process.
    ///
    /// Passing by reference is not only for doing less copying, 
    /// one might want to attempt to fix invalid fields 
    /// in onError() methods.
    void validate(ref BamRead alignment) {
        _visitAlignment(alignment);
    }

    /** Implement those methods to define your own behaviour.

        During validation process, in case of an error the corresponding
        method gets called, and is provided the object where the error occurred,
        and type of the error. Objects are passed by reference so that they
        can be changed (fixed / cleaned up / etc.)

        If onError() returns true, that means to continue validation process
        for this particular entity. Otherwise, all other validation checks are
        skipped and next entity is processed.
    */
    abstract bool onError(ref BamRead al, AlignmentError error); 
    abstract bool onError(ref BamRead al, CigarError error); /// ditto
    abstract bool onError(string key, const ref Value value, TagError error); /// ditto

private:

    // Method names are a bit misleading,
    // their return value is NOT whether a field is invalid or not
    // but rather whether onError() handlers decide to stop validation
    // when the field is invalid.

    bool invalidReadName(ref BamRead al) {
        // Read name (a.k.a. QNAME) must =~ /^[!-?A-~]{1,255}$/ 
        // according to specification.
        if (al.name.length == 0) {
            if (!onError(al, AlignmentError.EmptyReadName)) return true;
        } else if (al.name.length > 255) {
            if (!onError(al, AlignmentError.TooLongReadName)) return true;
        } else {
            foreach (char c; al.name) 
            {
                if ((c < '!') || (c > '~') || (c == '@')) {
                    if (!onError(al, AlignmentError.ReadNameContainsInvalidCharacters)) {
                        return true;
                    } else {
                        break;
                    }
                }
            }
        }
        return false;
    }

    bool invalidPosition(ref BamRead al) {
        /// Check that position is in range [-1 .. 2^29 - 2]
        if (al.position < -1 || al.position > ((1<<29) - 2)) {
            if (!onError(al, AlignmentError.PositionIsOutOfRange)) {
                return true;
            }
        }
        return false;
    }

    bool invalidQualityData(ref BamRead al) {
        /// Check quality data
        if (!all!"a == 0xFF"(al.base_qualities) &&
            !all!"0 <= a && a <= 93"(al.base_qualities)) 
        {
            if (!onError(al, AlignmentError.QualityDataContainsInvalidElements)) {
                return true;
            }
        }
        return false;
    }

    static bool internalHardClipping(ref BamRead al) {
        return (al.cigar.length > 2 && 
                any!"a.type == 'H'"(al.cigar[1..$-1]));
    }

    static bool internalSoftClipping(ref BamRead al) {
        if (al.cigar.length <= 2) return false;

        auto cigar = al.cigar;

        /// strip H operations from ends
        if (cigar[0].type == 'H') {
            cigar = cigar[1..$];
        }
        if (cigar[$-1].type == 'H') {
            cigar = cigar[0..$-1];
        }

        /// check that S operations are at the ends only
        return (cigar.length > 2 &&
                any!"a.type == 'S'"(cigar[1..$-1]));
    } 

    //  Sum of M/I/S/=/X operations must be equal to the sequence length
    //  if both sequence and CIGAR string are presented.
    static bool inconsistentLength(ref BamRead al) {
        return (al.sequence_length > 0 &&
                al.sequence_length != reduce!`a + b`(0, 
                                        map!`a.length`(
                                          filter!`canFind("MIS=X", a.type)`(
                                            al.cigar))));
    }
 
    bool invalidCigar(ref BamRead al) {
        
        if (al.cigar.length == 0) return false;

        static string check(string s) {
            import std.ascii : toUpper;
            return (`if (`~s.dup~`(al)`~
                   `    && !onError(al, CigarError.`~(cast(char)(s[0]-32))~s[1..$]~`)`~
                   `    && (called_on_error || onError(al, AlignmentError.InvalidCigar)))`~
                   `{`~
                   `    return true;`~
                   `}`).idup;
        }

        bool called_on_error = false;

        mixin(check("internalHardClipping"));
        mixin(check("internalSoftClipping"));
        mixin(check("inconsistentLength"));

        return false;
    }
   
    // Check tags, a lot of them are predefined in the specification
    // and have to satisfy certain requirements.
    bool invalidTags(ref BamRead al) {

        bool all_tags_are_good = true;
        
        void someTagIsBad() {
            if (all_tags_are_good) {
                if (!onError(al, AlignmentError.InvalidTags)) return;
            }
            all_tags_are_good = false;
        }

        /// Check that all tag keys are distinct.

        bool all_distinct = true;

        // Optimize for small number of tags
        ushort[256] keys = void;
        size_t i = 0;

        // Check each tag in turn.
        foreach (k, v; al) {
            if (!isValid(k, v, al)) {
                someTagIsBad();
            }
           
            if (i < keys.length) {
                keys[i] = *cast(ushort*)(k.ptr);

                if (all_distinct) {
                    for (size_t j = 0; j < i; ++j) {
                        if (keys[i] == keys[j]) {
                            all_distinct = false;
                            break;
                        }
                    }
                }

                i += 1;
            } else {
                if (all_distinct) {
                    // must be exactly one
                    int found = 0;
                    foreach (k2, v2; al) {
                        if (*cast(ushort*)(k2.ptr) == *cast(ushort*)(k.ptr)) {
                            if (found == 1) {
                                all_distinct = false;
                                break;
                            } else {
                                ++found;
                            }
                        }
                    }
                }
            }
        }
       
        if (!all_distinct) {
            if (!onError(al, AlignmentError.DuplicateTagKeys)) return true;
        }

        return false;
    }

    void _visitAlignment(ref BamRead al) {
        if (invalidReadName(al)) return;
        if (invalidPosition(al)) return;
        if (invalidQualityData(al)) return;
        if (invalidCigar(al)) return;
        if (invalidTags(al)) return;
    }

    bool isValid(string key, const ref Value value, const ref BamRead al) {

        bool result = true;

        if (value.is_hexadecimal_string()) {
            auto str = cast(string)value;
            if (str.length == 0) {
                if (!onError(key, value, TagError.EmptyHexadecimalString)) {
                    return false;
                }
                result = false;
            }
            /// check that it contains only 0..9a-fA-F characters
            if (!all!(isHexDigit)(str)) {
                if (!onError(key, value, TagError.InvalidHexadecimalString)) {
                    return false;
                }
                result = false;
            }
        } else if (value.is_character()) {
            /// character must be [!-~]
            auto c = cast(char)value;
            if (!(c >= '!' && c <= '~')) {
                if (!onError(key, value, TagError.NonPrintableCharacter)) {
                    return false;
                }
                result = false;
            }
        } else if (value.is_string()) {
            auto str = cast(string)value; 
            if (str.length == 0) {
                if (!onError(key, value, TagError.EmptyString)) {
                    return false;
                }
                result = false;
            }
            /// string must be [ !-~]+
            if (!all!"a >= ' ' && a <= '~'"(str)) {
                if (!onError(key, value, TagError.NonPrintableString)) {
                    return false;
                }
                result = false;
            }
        }

        /// check various tags from SAM/BAM specification
        if (!additionalChecksIfTheTagIsPredefined(key, value, al)) {
            result = false;
        }

        return result;
    }

    // There're some requirements for predefined tags to be checked
    // such as type, length in some cases, or even some regular expression.
    // See page 6 of SAM/BAM specification.
    bool additionalChecksIfTheTagIsPredefined(string key, const ref Value value,
                                              const ref BamRead al) 
    {
        bool result = true;

        // Creates a switch for all predefined tag keys.
        string switchTagKey() {
            char[] cs;
            foreach (t; PredefinedTags) {
                cs ~= `case "`~t.Key~`":`~
                      `  if (!checkTagValue!"`~t.Key~`"(value, al)) {`~
                      `    result = false;`~
                      `  }`~
                      `  break;`.dup;
            }
            return "switch (key) { " ~ cs.idup ~ " default : break; }";
        }

        mixin(switchTagKey());

        return result;
    }

    // Supposed to be inlined in the above switch
    bool checkTagValue(string s)(const ref Value value, const ref BamRead al) {

        bool result = true;
        
        /// 1. Check type.

        static if (is(PredefinedTagType!s == int)) {
            if (!value.is_integer) {
                if (!onError(s, value, TagError.ExpectedIntegerValue)) {
                    return false;
                }
                result = false;
            }
        } else if (is(PredefinedTagType!s == string)) {
            // Notice that there are no 'H'-typed predefined tags,
            // and they are almost unused and therefore are not likely
            // to appear in the future. 
            if (!value.is_string || value.bam_typeid == 'H') {
                if (!onError(s, value, TagError.ExpectedStringValue)) {
                    return false;
                }
                result = false;
            }
        } else {
            if (value.tag != GetTypeId!(PredefinedTagType!s)) {
                if (!onError(s, value, TagError.InvalidValueType)) {
                    return false;
                }
                result = false;
            }
        }
        
        /// 2. For tags which contain quality as a string,
        ///    check that all characters are valid 
        
        static if (staticIndexOf!(s, "CQ", "E2", "OQ", "Q2", "U2") != -1) {
            auto str = cast(string)value;
            if (str != "*" && !all!"a >= '!' && a <= '~'"(str)) {
                if (!onError(s, value, TagError.InvalidQualityString)) {
                    return false;
                }
                result = false;
            }
        }

        /// 3. In a couple of cases values are required to be 
        ///    of the same length as the read sequence.

        static if (staticIndexOf!(s, "BQ", "E2") != -1) {
            if ((cast(string)value).length != al.sequence_length) {
                if (!onError(s, value, TagError.ExpectedStringWithSameLengthAsSequence)) {
                    return false;
                }
            }
        }


        /// 4. MD tag ought to: a) match /^[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*$/
        ///                     b) match CIGAR string (TODO?)
        
        static if (s == "MD") {

            /// a) check regular expression

            auto s = cast(string)value;
            bool valid = true;
            if (s.length == 0) valid = false;
            if (!isDigit(s[0])) valid = false;
            size_t i = 1;
            while (i < s.length && isDigit(s[i])) 
                ++i;
            while (i < s.length) {
                if (isUpper(s[i])) {
                    ++i; // [A-Z]
                } else if (s[i] == '^') { // ^[A-Z]+
                    ++i;
                    if (i == s.length || !isUpper(s[i])) {
                        valid = false;
                        break;
                    }
                    while (i < s.length && isUpper(s[i]))
                        ++i;
                } else {
                    valid = false;
                    break;
                }
                // now [0-9]+
                if (i == s.length || !isDigit(s[i])) {
                    valid = false;
                    break;
                }
                while (i < s.length && isDigit(s[i]))
                    ++i;
            }

            if (i < s.length) {
                valid = false;
            }

            if (!valid) result = false;
        }

        return result;
    }
}

final private class BooleanValidator : AbstractAlignmentValidator {

    bool result;

    override void validate(ref BamRead al) {
        result = true;
        super.validate(al);
    }

    override bool onError(ref BamRead al, AlignmentError e) {
        return (result = false);
    }

    override bool onError(ref BamRead al, CigarError e) {
        return (result = false);
    }

    override bool onError(string key, const ref Value val, TagError e) {
        return (result = false);
    }

}

private static BooleanValidator booleanValidator;

static this() {
    booleanValidator = new BooleanValidator();
}

/// Check if alignment is valid
bool isValid(BamRead alignment) { 
    booleanValidator.validate(alignment);
    return booleanValidator.result;
}
