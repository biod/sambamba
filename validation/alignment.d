module validation.alignment;

public import alignment;
public import tagvalue;
import utils.algo;

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
    void validate(ref Alignment alignment) {
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
    abstract bool onError(ref Alignment al, AlignmentError error); 
    abstract bool onError(ref Alignment al, CigarError error); /// ditto
    abstract bool onError(string key, ref Value value, TagError error); /// ditto

private:
    void _visitAlignment(ref Alignment al) {

        /// Read name (a.k.a. QNAME) must =~ /^[!-?A-~]{1,255}$/ 
        /// according to specification.
        if (al.read_name.length == 0) {
            if (!onError(al, AlignmentError.EmptyReadName)) return;
        } else if (al.read_name.length > 255) {
            if (!onError(al, AlignmentError.TooLongReadName)) return;
        } else {
            if (!all!"(a >= '!' && a <= '?') || (a >= 'A' && a <= '~')"(al.read_name)) 
            {
                if (!onError(al, AlignmentError.ReadNameContainsInvalidCharacters)) return;
            }
        }

        /// Check that position is in range [-1 .. 2^29 - 2]
        if (al.position < -1 || al.position > ((1<<29) - 2)) {
            if (!onError(al, AlignmentError.PositionIsOutOfRange)) return;
        }

        /// Check quality data
        if (!all!"a == 0xFF"(al.phred_base_quality) &&
            !all!"0 <= a && a <= 93"(al.phred_base_quality)) 
        {
            if (!onError(al, AlignmentError.QualityDataContainsInvalidElements)) return;
        }
        
        /// Check CIGAR string
        
        if (al.cigar.length > 0) {

        bool cigar_is_good = true;
        
        void cigarIsBad() {
            if (cigar_is_good) {
                if (!onError(al, AlignmentError.InvalidCigar)) return;
            }
            cigar_is_good = false;
        }

            /// 1. H may only be present as first/last operation.
            if (al.cigar.length > 2 && 
                any!"a.operation == 'H'"(al.cigar[1..$-1]))
            {
                if (!onError(al, CigarError.InternalHardClipping)) return; 
                cigarIsBad(); 
            }

            /// 2. The same holds for S operations except that
            ///    H may be before or after them.
            if (al.cigar.length > 2) {
                auto cigar = al.cigar;

                /// strip H operations from ends
                if (cigar[0].operation == 'H') {
                    cigar = cigar[1..$];
                }
                if (cigar[$-1].operation == 'H') {
                    cigar = cigar[0..$-1];
                }

                /// check that S operations are at the ends only
                if (cigar.length > 2 &&
                    any!"a.operation == 'S'"(cigar[1..$-1]))
                {
                    if (!onError(al, CigarError.InternalSoftClipping)) return;    
                    cigarIsBad();
                }
            }
            
            /// 3. Sum of M/I/S/=/X operations shall equal sequence length
            ///    if both sequence and CIGAR string are presented.

            if (al.sequence_length > 0 &&
                al.sequence_length != reduce!`a + b`(0, 
                                        map!`a.length`(
                                          filter!`canFind("MIS=X", a.operation)`(
                                            al.cigar))))
            {
                if (!onError(al, CigarError.InconsistentLength)) return;
                cigarIsBad();
            }

        } 
        /// end of CIGAR checking

        //-----------------------------------------------------------------

        /// Check tags, a lot of them are predefined in the specification
        /// and have to satisfy certain requirements.
        
        bool all_tags_are_good = true;
        
        void someTagIsBad() {
            if (all_tags_are_good) {
                if (!onError(al, AlignmentError.InvalidTags)) return;
            }
            all_tags_are_good = false;
        }

        string[] keys;

        /// Check each tag in turn.
        foreach (k, v; al.tags) {
            if (!isValid(k, v, al)) {
                someTagIsBad();
            }

            keys ~= k;
        }
       
        /// Check that all tag keys are distinct.
        if (!allDistinct(keys)) {
            if (!onError(al, AlignmentError.DuplicateTagKeys)) return;
        }

    }

    bool isValid(string key, Value value, ref Alignment al) {

        bool result = true;

        if (value.is_hexadecimal_string()) {
            auto str = to!string(value);
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
            auto c = to!char(value);
            if (!(c >= '!' && c <= '~')) {
                if (!onError(key, value, TagError.NonPrintableCharacter)) {
                    return false;
                }
                result = false;
            }
        } else if (value.is_string()) {
            auto str = to!string(value); 
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
    bool additionalChecksIfTheTagIsPredefined(string key, Value value,
                                              ref Alignment al) 
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
    bool checkTagValue(string s)(Value value, ref Alignment al) {

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
            auto str = to!string(value);
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
            if (to!string(value).length != al.sequence_length) {
                if (!onError(s, value, TagError.ExpectedStringWithSameLengthAsSequence)) {
                    return false;
                }
            }
        }


        /// 4. MD tag ought to: a) match /^[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*$/
        ///                     b) match CIGAR string (TODO?)
        
        static if (s == "MD") {

            /// a) check regular expression

            auto s = to!string(value);
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

    final override void validate(ref Alignment al) {
        result = true;
        super.validate(al);
    }

    bool onError(ref Alignment al, AlignmentError e) {
        return (result = false);
    }

    bool onError(ref Alignment al, CigarError e) {
        return (result = false);
    }

    bool onError(string key, ref Value val, TagError e) {
        return (result = false);
    }

}

/// Check if alignment is valid
bool isValid(Alignment alignment) {
    scope validator = new BooleanValidator();
    validator.validate(alignment);
    return validator.result;
}
