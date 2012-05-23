/**
  Module for SAM header validation.
  
  In order to implement your own validation behaviour,
  subclass SamHeaderValidator and define your own 
  onError() methods.
*/
module validation.samheader;

import samheader;
import utils.algo : allDistinct;

import std.algorithm;
import std.functional;
import std.ascii;

/// High-level SAM header errors
enum SamHeaderError {
    InvalidSqLine,
    InvalidPgLine,
    InvalidRgLine,
    InvalidFormatVersion,
    InvalidSortingOrder,
    NotUniqueSequenceNames,
    NotUniqueReadGroupIdentifiers,
    NotUniqueProgramIdentifiers
}

/// Errors in @SQ lines
enum SqLineError {
    EmptySequenceName,
    InvalidSequenceName,
    SequenceLengthOutOfRange
}

/// Errors in @RG lines
enum RgLineError {
    UnknownPlatform
}

/// Errors in @PG lines
enum PgLineError {
    NoMatchForPreviousProgram
}

/**
  Abstract class encapsulating visitation of SAM header elements.
*/
abstract class SamHeaderValidator {

    /// Passing by reference is not only for doing less copying, 
    /// one might want to attempt to fix invalid fields 
    /// in onError() methods (see below).
    this(ref SamHeader header) {
        _header = header;
    }

    /// Start validation process. Override if you need to.
    void validate() {
        _visitHeader();
    }

    /** Implement those methods to define your own behaviour.

        During validation process, in case of error the corresponding
        method gets called, and is provided the object where error occurred,
        and type of error. Objects are passed by reference so that they
        can be changed (fixed / cleaned / etc.)

        Example:
        ----------------------

        class BooleanValidator : SamHeaderValidator {
            bool result = true;
            this(ref SamHeader header) {
                super(header);
            }
            void onError(T, U)(T t, U u) {
                result = false;
            }
        }

        bool isValid(ref SamHeader header) {
            auto validator = new BooleanValidator(header);
            validator.validate();
            return validator.result;
        }
    */
    abstract void onError(ref SamHeader header, SamHeaderError error);
    abstract void onError(ref SqLine line, SqLineError error); /// ditto
    abstract void onError(ref PgLine line, PgLineError error); /// ditto
    abstract void onError(ref RgLine line, RgLineError error); /// ditto

private:
    SamHeader _header;

    bool isValid(ref SqLine sq) {

        bool result = true;

        if (sq.sequence_name.length == 0) {
            onError(sq, SqLineError.EmptySequenceName);
            result = false;
        }

        // @SQ/LN must be in range 1 .. (1<<29)-1
        // (sequence_length is uint)
        if (sq.sequence_length == 0 || sq.sequence_length >= (1<<29)) 
        {
            onError(sq, SqLineError.SequenceLengthOutOfRange);
            result = false;
        }

        // check that sequence_name is /^[!-)+-<>-~][!-~]*$/
        auto first = sq.sequence_name[0];
        if (!((first >= '!' && first <= ')') ||
              (first >= '+' && first <= '<') ||
              (first >= '>' && first <= '~'))) 
        {
            onError(sq, SqLineError.InvalidSequenceName);
            result = false;
        }
        
        if (!all!"a >= '!' && a <= '~'"(sq.sequence_name[1..$])) {
            onError(sq, SqLineError.InvalidSequenceName);
            result = false;
        }

        return result;
    }

    bool isValid(ref RgLine rg) {
        bool res = canFind(["ILLUMINA",
                            "SOLID",
                            "LS454",
                            "HELICOS",
                            "PACBIO"],
                           rg.platform);
        if (!res) {
            onError(rg, RgLineError.UnknownPlatform);
        }

        return res;
    }

    bool isValid(ref PgLine pg) {

        // checking PP tag occurs in visit() 
        // because it involves other @PG lines

        return true;
    }

    void _visitHeader() {

        foreach (sq; _header.sq_lines) {
            if (!isValid(sq)) onError(_header, SamHeaderError.InvalidSqLine);
        }

        foreach (rg; _header.rg_lines) {
            if (!isValid(rg)) onError(_header, SamHeaderError.InvalidRgLine);
        }

        foreach (pg; _header.pg_lines) {
            if (!isValid(pg)) onError(_header, SamHeaderError.InvalidPgLine);
        }

        if (_header.hasHeaderLine()) {
            if (!checkFormatVersion(_header.format_version)) {
                onError(_header, SamHeaderError.InvalidFormatVersion);
            }
            if (!checkSortingOrder(_header.sorting_order)) {
                onError(_header, SamHeaderError.InvalidSortingOrder);
            }
        }

        // check uniqueness of @SQ/SN
        if (!allDistinct(map!"a.sequence_name"(_header.sq_lines))) {
            onError(_header, SamHeaderError.NotUniqueSequenceNames);
        }

        // check uniqueness of @RG/ID
        if (!allDistinct(map!"a.identifier"(_header.rg_lines))) {
            onError(_header, SamHeaderError.NotUniqueReadGroupIdentifiers);
        }

        // check uniqueness of @PG/ID
        if (!allDistinct(map!"a.identifier"(_header.pg_lines))) {
            onError(_header, SamHeaderError.NotUniqueProgramIdentifiers);
        }

        // check that each @PG/PP matches some @PG/ID
        foreach (pg; _header.pg_lines) {
            if (pg.previous_program.length != 0) {
                if (!canFind(map!"a.identifier"(_header.pg_lines),
                             pg.previous_program)) 
                {
                    onError(pg, PgLineError.NoMatchForPreviousProgram);
                }
            }
        }
    } // visitHeader

} // SamHeaderValidator

private {

/// check that @HD/VN is /^[0-9]+\.[0-9]+$/ 
bool checkFormatVersion(string ver) nothrow {

    if (ver.length == 0) {
        return false; // must be non-empty
    }

    if (!isDigit(ver[0])) {
        return false; // and it must start with digit
    }

    ver = ver[1..$];

    bool passed_dot = false;

    while (ver.length > 0) {
        if (isDigit(ver[0])) {
            ver = ver[1..$]; //\\// skip digits
        } else if (ver[0] == '.') {
            if (passed_dot) {
                return false; //\\// must contain only one dot
            }
            passed_dot = true;
            ver = ver[1..$];
            if (ver.length == 0 || !isDigit(ver[0])) {
                return false; //\\// there must be a digit after dot
            }
        }
    }

    return true;
}

unittest {
    assert(checkFormatVersion("1.53") == true);
    assert(checkFormatVersion("a.23") == false);
    assert(checkFormatVersion("1.2.3") == false);
    assert(checkFormatVersion("5.") == false);
    assert(checkFormatVersion("3.141592653589793") == true);
    assert(checkFormatVersion("100500.42") == true);
    assert(checkFormatVersion("2.71828.3.5") == false);
}

/// check @HD/SO
bool checkSortingOrder(string sorting_order) {
    return canFind(["unknown",
                    "unsorted",
                    "queryname",
                    "coordinate"],
                   sorting_order);
}

/// The most simple validator possible
class BooleanValidator : SamHeaderValidator {

    class BooleanValidationException : Exception {
        this() { super(""); }
    }

    bool result = true;

    this(ref SamHeader header) {
        super(header);
    }

    override void validate() {
        try {
            super.validate();
        } catch (BooleanValidationException e) {
            result = false;
        }
    }

    void onError(ref SamHeader header, SamHeaderError e) {
        throw new BooleanValidationException();
    }

    void onError(ref SqLine line, SqLineError e) {
        throw new BooleanValidationException();
    }

    void onError(ref RgLine header, RgLineError e) {
        throw new BooleanValidationException();
    }

    void onError(ref PgLine header, PgLineError e) {
        throw new BooleanValidationException();
    }
}

} // private

/// Check if header is valid
bool isValid(ref SamHeader header) {
    auto validator = new BooleanValidator(header);
    validator.validate();
    return validator.result;
}

unittest {
    auto header = SamHeader("@HD\tVN:1.3\tSO:coordinate\n@SQ\tSN:chr1\tLN:1575\n@SQ\tSN:chr2\tLN:1584");
    assert(isValid(header));

    auto empty_seq_name = SamHeader("@HD\tVN:1.3\tSO:coordinate\n@SQ\tSN:chr1\tLN:1575\n@SQ\tSN:\tLN:1584");
    assert(!isValid(empty_seq_name));

}
