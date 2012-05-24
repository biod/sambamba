module validation.alignment;

public import alignment;
public import tagvalue;
import utils.algo;

import std.algorithm;
import std.ascii;
import std.conv;

/**
    Alignment validation error types.

    InvalidCigar error is accompanied by some CigarError,
    InvalidTag is accompanied by some TagError.
*/
enum AlignmentError {
    EmptyReadName,
    TooLongReadName,
    ReadNameContainsInvalidCharacters,
    PositionIsOutOfRange,
    QualityDataContainsInvalidElements,
    InvalidCigar,
    InvalidTag,
    DuplicateTagKeys
}

/// CIGAR string validation error types.
enum CigarError {
    InternalHardClipping,
    InternalSoftClipping,
    InconsistentLength
}

/// Auxiliary data validation error types.
///
/// Refers to an individual tag.
enum TagError {
    EmptyString,
    EmptyHexadecimalString,
    NonPrintableString,
    NonPrintableCharacter,
    InvalidHexadecimalString
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
    */
    abstract void onError(ref Alignment al, AlignmentError error); // ditto
    abstract void onError(ref Alignment al, CigarError error); // ditto
    abstract void onError(string key, ref Value value, TagError error); // ditto

private:
    void _visitAlignment(ref Alignment al) {

        /// Read name (a.k.a. QNAME) must =~ /^[!-?A-~]{1,255}$/ 
        /// according to specification.
        if (al.read_name.length == 0) {
            onError(al, AlignmentError.EmptyReadName);
        } else if (al.read_name.length > 255) {
            onError(al, AlignmentError.TooLongReadName);
        } else {
            if (!all!"(a >= '!' && a <= '?') || (a >= 'A' && a <= '~')"(al.read_name)) 
            {
                onError(al, AlignmentError.ReadNameContainsInvalidCharacters);
            }
        }

        /// Check that position is in range [-1 .. 2^29 - 2]
        if (al.position < -1 || al.position > ((1<<29) - 2)) {
            onError(al, AlignmentError.PositionIsOutOfRange);
        }

        /// Check quality data
        if (!all!"a == 0xFF"(al.phred_base_quality) &&
            !all!"0 <= a && a <= 93"(al.phred_base_quality)) 
        {
            onError(al, AlignmentError.QualityDataContainsInvalidElements);
        }
        
        /// Check CIGAR string
        
        if (al.cigar.length > 0) {

        bool cigar_is_good = true;
        
        void cigarIsBad() {
            if (cigar_is_good) {
                onError(al, AlignmentError.InvalidCigar);
            }
            cigar_is_good = false;
        }

            /// 1. H may only be present as first/last operation.
            if (al.cigar.length > 2 && 
                any!"a.operation == 'H'"(al.cigar[1..$-1]))
            {
                onError(al, CigarError.InternalHardClipping); 
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
                    onError(al, CigarError.InternalSoftClipping);    
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
                onError(al, CigarError.InconsistentLength);
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
                onError(al, AlignmentError.InvalidTag);
            }
            all_tags_are_good = false;
        }

        string[] keys;

        /// Check each tag in turn.
        foreach (k, v; al.tags) {
            if (!isValid(k, v)) {
                someTagIsBad();
            }

            keys ~= k;
        }
       
        /// Check that all tag keys are distinct.
        if (!allDistinct(keys)) {
            onError(al, AlignmentError.DuplicateTagKeys);
        }

    }

    bool isValid(string key, Value value) {

        bool result = true;

        if (value.is_hexadecimal_string()) {
            auto str = to!string(value);
            if (str.length == 0) {
                onError(key, value, TagError.EmptyHexadecimalString);
            }
            /// check that it contains only 0..9a-fA-F characters
            if (!all!(isHexDigit)(str)) {
                onError(key, value, TagError.InvalidHexadecimalString);
                result = false;
            }
        } else if (value.is_character()) {
            /// character must be [!-~]
            auto c = to!char(value);
            if (!(c >= '!' && c <= '~')) {
                onError(key, value, TagError.NonPrintableCharacter);
                result = false;
            }
        } else if (value.is_string()) {
            auto str = to!string(value); 
            if (str.length == 0) {
                onError(key, value, TagError.EmptyString);
                result = false;
            }
            /// string must be [ !-~]+
            if (!all!"a >= ' ' && a <= '~'"(str)) {
                onError(key, value, TagError.NonPrintableString);
                result = false;
            }
        }
        
        /// TODO: add checks for various tags from SAM/BAM specification

        return result;
    }
}

private class BooleanValidator : AbstractAlignmentValidator {

    class BooleanValidationException : Exception {
        this() { super(""); }
    }

    bool result;

    override void validate(ref Alignment al) {
        result = true;
        try {
            super.validate(al);
        } catch (BooleanValidationException e) {
            result = false;
        }
    }

    void onError(ref Alignment al, AlignmentError e) {
        throw new BooleanValidationException();
    }

    void onError(ref Alignment al, CigarError e) {
        throw new BooleanValidationException();
    }

    void onError(string key, ref Value val, TagError e) {
        throw new BooleanValidationException();
    }

}

/// Check if alignment is valid
bool isValid(Alignment alignment) {
    auto validator = new BooleanValidator();
    validator.validate(alignment);
    return validator.result;
}
