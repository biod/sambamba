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
/**
  Module for SAM header validation.
  
  In order to implement your own validation behaviour,
  subclass AbstractSamHeaderValidator and define your own 
  onError() methods.
*/
module bio.std.hts.bam.validation.samheader;

public import bio.std.hts.sam.header;
import bio.core.utils.algo;

import std.algorithm;
import std.functional;
import std.ascii;

/// SAM header validation error types.
///
/// Each Invalid??Line error is accompanied by 
/// corresponding ??LineError.
enum SamHeaderError {
    InvalidSqLine,
    InvalidPgLine,
    InvalidRgLine,
    InvalidFormatVersion
}

/// @SQ line validation error types.
enum SqLineError {
    MissingSequenceName,
    InvalidSequenceName,
    SequenceLengthOutOfRange
}

/// @RG line validation error types.
enum RgLineError {
    UnknownPlatform,
    MissingIdentifier
}

/// @PG line validation error types.
enum PgLineError {
    NoMatchForPreviousProgram,
    MissingIdentifier
}

/**
  Abstract class encapsulating visitation of SAM header elements.
*/
abstract class AbstractSamHeaderValidator {

    /// Start validation process.
    ///
    /// Passing by reference is not only for doing less copying, 
    /// one might want to attempt to fix invalid fields 
    /// in onError() methods.
    void validate(ref SamHeader header) {
        _visitHeader(header);
    }

    /** Implement those methods to define your own behaviour.

        During validation process, in case of an error the corresponding
        method gets called, and is provided the object where the error occurred,
        and type of the error. Objects are passed by reference so that they
        can be changed (fixed / cleaned up / etc.)

        'False' return value means to stop further validation checks for the 
        current entity and skip to the next one.
    */
    abstract bool onError(ref SamHeader header, SamHeaderError error);
    abstract bool onError(ref SqLine line, SqLineError error); /// ditto
    abstract bool onError(ref PgLine line, PgLineError error); /// ditto
    abstract bool onError(ref RgLine line, RgLineError error); /// ditto

private:

    bool isValid(ref SqLine sq) {

        /// All members of SqLine get initialized.
        /// Initial value for name is an empty string,
        /// and for sequence_length is 0

        bool result = true;

        if (sq.name.length == 0) {
            onError(sq, SqLineError.MissingSequenceName);
            result = false;
        } else {
            // check that sequence_name is /^[!-)+-<>-~][!-~]*$/
            auto first = sq.name[0];
            if (!((first >= '!' && first <= ')') ||
                  (first >= '+' && first <= '<') ||
                  (first >= '>' && first <= '~'))) 
            {
                onError(sq, SqLineError.InvalidSequenceName);
                result = false;
            }
            
            if (!all!"a >= '!' && a <= '~'"(sq.name[1..$])) {
                onError(sq, SqLineError.InvalidSequenceName);
                result = false;
            }
        }

        // @SQ/LN must be in range 1 .. (1<<29)-1
        // (sequence_length is uint)
        if (sq.length == 0 || sq.length >= (1<<29)) 
        {
            onError(sq, SqLineError.SequenceLengthOutOfRange);
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
        
        if (rg.identifier.length == 0) {
            onError(rg, RgLineError.MissingIdentifier);
            res = false;
        }

        return res;
    }

    bool isValid(ref PgLine pg) {

        // checking PP tag occurs in _visitHeader()
        // because it involves other @PG lines
        
        if (pg.identifier.length == 0) {
            onError(pg, PgLineError.MissingIdentifier);
            return false;
        }
        
        return true;
    }

    void _visitHeader(ref SamHeader header) {

        foreach (sq; header.sequences) {
            if (!isValid(sq)) if (!onError(header, SamHeaderError.InvalidSqLine)) return;
        }

        foreach (rg; header.read_groups) {
            if (!isValid(rg)) if (!onError(header, SamHeaderError.InvalidRgLine)) return;
        }

        foreach (pg; header.programs) {
            if (!isValid(pg)) if (!onError(header, SamHeaderError.InvalidPgLine)) return;
        }

        if (!checkFormatVersion(header.format_version)) {
            if (!onError(header, SamHeaderError.InvalidFormatVersion)) return;
        }

        // uniqueness of @SQ/SN, @RG/ID, and @PG/ID
        // is guaranteed by design of HeaderLineDictionary template class

        // check that each @PG/PP matches some @PG/ID
        foreach (pg; header.programs) {
            if (pg.previous_program.length != 0) {
                if (!canFind(map!"a.identifier"(header.programs.values),
                             pg.previous_program)) 
                {
                    if (!onError(pg, PgLineError.NoMatchForPreviousProgram)) return;
                }
            }
        }
    } // visitHeader

} // AbstractSamHeaderValidator

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
            ver = ver[1..$]; // skip digits
        } else if (ver[0] == '.') {
            if (passed_dot) {
                return false; // must contain only one dot
            }
            passed_dot = true;
            ver = ver[1..$];
            if (ver.length == 0 || !isDigit(ver[0])) {
                return false; // there must be a digit after dot
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

final private class BooleanValidator : AbstractSamHeaderValidator {

    bool result;

    override void validate(ref SamHeader header) {
        result = true;
        super.validate(header);
    }

    override bool onError(ref SamHeader header, SamHeaderError e) {
        return (result = false);
    }

    override bool onError(ref SqLine line, SqLineError e) {
        return (result = false);
    }

    override bool onError(ref RgLine header, RgLineError e) {
        return (result = false);
    }

    override bool onError(ref PgLine header, PgLineError e) {
        return (result = false);
    }
}

static BooleanValidator booleanValidator;

} // private

static this() {
    booleanValidator = new BooleanValidator();
}

/// Check if header is valid
bool isValid(SamHeader header) {
    booleanValidator.validate(header);
    return booleanValidator.result;
}

unittest {
    auto valid_header = new SamHeader("@HD\tVN:1.3\tSO:coordinate\n@SQ\tSN:chr1\tLN:1575");
    assert(isValid(valid_header));

    auto empty_seq_name = new SamHeader("@HD\tVN:1.3\tSO:coordinate\n@SQ\tSN:\tLN:1575");
    assert(!isValid(empty_seq_name));

    auto missing_seq_name = new SamHeader("@HD\tVN:1.3\tSO:coordinate\n@SQ\tLN:1575");
    assert(!isValid(missing_seq_name));

    auto missing_seq_length = new SamHeader("@HD\tVN:1.3\tSO:coordinate\n@SQ\tSN:chr1");
    assert(!isValid(missing_seq_length));

    auto seq_length_out_of_range = new SamHeader("@HD\tVN:1.3\tSO:coordinate\n@SQ\tSN:chr1\tLN:876543210");
    assert(!isValid(seq_length_out_of_range));

    auto invalid_seq_name = new SamHeader("@HD\tVN:1.3\tSO:coordinate\n@SQ\tSN:chr \tLN:1575");
    assert(!isValid(invalid_seq_name));

    auto missing_version = new SamHeader("@HD\tSO:coordinate");
    assert(!isValid(missing_version));

    auto invalid_version_format = new SamHeader("@HD\tVN:6.7.8");
    assert(!isValid(invalid_version_format));

    auto unknown_platform = new SamHeader("@RG\tID:678\tPL:TROLOLO");
    assert(!isValid(unknown_platform));

    auto missing_rg_id = new SamHeader("@RG\tPL:ILLUMINA");
    assert(!isValid(missing_rg_id));

    auto missing_pg_id = new SamHeader("@PG\tPN:bwa\tVN:0.5.9-r16");
    assert(!isValid(missing_pg_id));

    auto unknown_previous_program = new SamHeader("@PG\tID:bwa_aln_fastq\tPN:bwa\tPP:bwa_index");
    assert(!isValid(unknown_previous_program));

    auto another_valid_header = new SamHeader(q"[@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:249250621	M5:1b22b98cdeb4a9304cb5d48026a85128	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:2	LN:243199373	M5:a0d9851da00400dec1098a9255ac712e	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:3	LN:198022430	M5:fdfd811849cc2fadebc929bb925902e5	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:4	LN:191154276	M5:23dccd106897542ad87d2765d28a19a1	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:5	LN:180915260	M5:0740173db9ffd264d728f32784845cd7	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:6	LN:171115067	M5:1d3a93a248d92a729ee764823acbbc6b	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:7	LN:159138663	M5:618366e953d6aaad97dbe4777c29375e	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:8	LN:146364022	M5:96f514a9929e410c6651697bded59aec	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:9	LN:141213431	M5:3e273117f15e0a400f01055d9f393768	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:10	LN:135534747	M5:988c28e000e84c26d552359af1ea2e1d	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:11	LN:135006516	M5:98c59049a2df285c76ffb1c6db8f8b96	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:12	LN:133851895	M5:51851ac0e1a115847ad36449b0015864	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:13	LN:115169878	M5:283f8d7892baa81b510a015719ca7b0b	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:14	LN:107349540	M5:98f3cae32b2a2e9524bc19813927542e	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:15	LN:102531392	M5:e5645a794a8238215b2cd77acb95a078	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:16	LN:90354753	M5:fc9b1a7b42b97a864f56b348b06095e6	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:17	LN:81195210	M5:351f64d4f4f9ddd45b35336ad97aa6de	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:18	LN:78077248	M5:b15d4b2d29dde9d3e4f93d1d0f2cbc9c	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:19	LN:59128983	M5:1aacd71f30db8e561810913e0b72636d	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:20	LN:63025520	M5:0dec9660ec1efaaf33281c0d5ea2560f	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:21	LN:48129895	M5:2979a6085bfe28e3ad6f552f361ed74d	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:22	LN:51304566	M5:a718acaa6135fdca8357d5bfe94211dd	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:X	LN:155270560	M5:7e0e2e580297b7764e31dbc80c2540dd	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:Y	LN:59373566	M5:1fa3474750af0948bdf97d5a0ee52e51	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:MT	LN:16569	M5:c68f52674c9fb33aef52dcf399755519	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000207.1	LN:4262	M5:f3814841f1939d3ca19072d9e89f3fd7	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000226.1	LN:15008	M5:1c1b2cd1fccbc0a99b6a447fa24d1504	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000229.1	LN:19913	M5:d0f40ec87de311d8e715b52e4c7062e1	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000231.1	LN:27386	M5:ba8882ce3a1efa2080e5d29b956568a4	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000210.1	LN:27682	M5:851106a74238044126131ce2a8e5847c	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000239.1	LN:33824	M5:99795f15702caec4fa1c4e15f8a29c07	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000235.1	LN:34474	M5:118a25ca210cfbcdfb6c2ebb249f9680	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000201.1	LN:36148	M5:dfb7e7ec60ffdcb85cb359ea28454ee9	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000247.1	LN:36422	M5:7de00226bb7df1c57276ca6baabafd15	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000245.1	LN:36651	M5:89bc61960f37d94abf0df2d481ada0ec	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000197.1	LN:37175	M5:6f5efdd36643a9b8c8ccad6f2f1edc7b	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000203.1	LN:37498	M5:96358c325fe0e70bee73436e8bb14dbd	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000246.1	LN:38154	M5:e4afcd31912af9d9c2546acf1cb23af2	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000249.1	LN:38502	M5:1d78abec37c15fe29a275eb08d5af236	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000196.1	LN:38914	M5:d92206d1bb4c3b4019c43c0875c06dc0	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000248.1	LN:39786	M5:5a8e43bec9be36c7b49c84d585107776	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000244.1	LN:39929	M5:0996b4475f353ca98bacb756ac479140	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000238.1	LN:39939	M5:131b1efc3270cc838686b54e7c34b17b	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000202.1	LN:40103	M5:06cbf126247d89664a4faebad130fe9c	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000234.1	LN:40531	M5:93f998536b61a56fd0ff47322a911d4b	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000232.1	LN:40652	M5:3e06b6741061ad93a8587531307057d8	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000206.1	LN:41001	M5:43f69e423533e948bfae5ce1d45bd3f1	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000240.1	LN:41933	M5:445a86173da9f237d7bcf41c6cb8cc62	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000236.1	LN:41934	M5:fdcd739913efa1fdc64b6c0cd7016779	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000241.1	LN:42152	M5:ef4258cdc5a45c206cea8fc3e1d858cf	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000243.1	LN:43341	M5:cc34279a7e353136741c9fce79bc4396	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000242.1	LN:43523	M5:2f8694fc47576bc81b5fe9e7de0ba49e	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000230.1	LN:43691	M5:b4eb71ee878d3706246b7c1dbef69299	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000237.1	LN:45867	M5:e0c82e7751df73f4f6d0ed30cdc853c0	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000233.1	LN:45941	M5:7fed60298a8d62ff808b74b6ce820001	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000204.1	LN:81310	M5:efc49c871536fa8d79cb0a06fa739722	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000198.1	LN:90085	M5:868e7784040da90d900d2d1b667a1383	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000208.1	LN:92689	M5:aa81be49bf3fe63a79bdc6a6f279abf6	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000191.1	LN:106433	M5:d75b436f50a8214ee9c2a51d30b2c2cc	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000227.1	LN:128374	M5:a4aead23f8053f2655e468bcc6ecdceb	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000228.1	LN:129120	M5:c5a17c97e2c1a0b6a9cc5a6b064b714f	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000214.1	LN:137718	M5:46c2032c37f2ed899eb41c0473319a69	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000221.1	LN:155397	M5:3238fb74ea87ae857f9c7508d315babb	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000209.1	LN:159169	M5:f40598e2a5a6b26e84a3775e0d1e2c81	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000218.1	LN:161147	M5:1d708b54644c26c7e01c2dad5426d38c	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000220.1	LN:161802	M5:fc35de963c57bf7648429e6454f1c9db	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000213.1	LN:164239	M5:9d424fdcc98866650b58f004080a992a	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000211.1	LN:166566	M5:7daaa45c66b288847b9b32b964e623d3	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000199.1	LN:169874	M5:569af3b73522fab4b40995ae4944e78e	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000217.1	LN:172149	M5:6d243e18dea1945fb7f2517615b8f52e	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000216.1	LN:172294	M5:642a232d91c486ac339263820aef7fe0	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000215.1	LN:172545	M5:5eb3b418480ae67a997957c909375a73	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000205.1	LN:174588	M5:d22441398d99caf673e9afb9a1908ec5	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000219.1	LN:179198	M5:f977edd13bac459cb2ed4a5457dba1b3	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000224.1	LN:179693	M5:d5b2fc04f6b41b212a4198a07f450e20	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000223.1	LN:180455	M5:399dfa03bf32022ab52a846f7ca35b30	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000195.1	LN:182896	M5:5d9ec007868d517e73543b005ba48535	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000212.1	LN:186858	M5:563531689f3dbd691331fd6c5730a88b	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000222.1	LN:186861	M5:6fe9abac455169f50470f5a6b01d0f59	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000200.1	LN:187035	M5:75e4c8d17cd4addf3917d1703cacaf25	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000193.1	LN:189789	M5:dbb6e8ece0b5de29da56601613007c2a	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000194.1	LN:191469	M5:6ac8f815bf8e845bb3031b73f812c012	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000225.1	LN:211173	M5:63945c3e6962f28ffd469719a747e73c	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:GL000192.1	LN:547496	M5:325ba9e808f669dfeee210fdd7b470ac	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:NC_007605	LN:171823	M5:6743bd63b3ff2b5b8985d8933c53290a	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@SQ	SN:hs37d5	LN:35477943	M5:5b6a4b3a81a2d3c134b7d14bf6ad39f1	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz        AS:NCBI37       SP:Human
@RG	ID:ERR016155	LB:HUMgdtRAGDIAAPE	SM:HG00125	PI:488	CN:BGI	PL:ILLUMINA	DS:SRP001294
@RG	ID:ERR016156	LB:HUMgdtRAGDIAAPE	SM:HG00125	PI:489	CN:BGI	PL:ILLUMINA	DS:SRP001294
@RG	ID:ERR016157	LB:HUMgdtRAGDIAAPE	SM:HG00125	PI:488	CN:BGI	PL:ILLUMINA	DS:SRP001294
@PG	ID:bwa_index	PN:bwa	VN:0.5.9-r16	CL:bwa index -a bwtsw $reference_fasta
@PG	ID:bwa_aln_fastq	PN:bwa	PP:bwa_index	VN:0.5.9-r16	CL:bwa aln -q 15 -f $sai_file $reference_fasta $fastq_file
@PG	ID:bwa_sam	PN:bwa	PP:bwa_aln_fastq	VN:0.5.9-r16	CL:bwa sampe -a 1464 -r $rg_line -f $sam_file $reference_fasta $sai_file(s) $fastq_file(s)
@PG	ID:bwa_sam.1	PN:bwa	PP:bwa_aln_fastq	VN:0.5.9-r16	CL:bwa sampe -a 1467 -r $rg_line -f $sam_file $reference_fasta $sai_file(s) $fastq_file(s)
@PG	ID:sam_to_fixed_bam	PN:samtools	PP:bwa_sam	VN:0.1.17 (r973:277)	CL:samtools view -bSu $sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - $reference_fasta > $fixed_bam_file
@PG	ID:sam_to_fixed_bam.1	PN:samtools	PP:bwa_sam.1	VN:0.1.17 (r973:277)	CL:samtools view -bSu $sam_file | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd -u - $reference_fasta > $fixed_bam_file
@PG	ID:gatk_target_interval_creator	PN:GenomeAnalysisTK	PP:sam_to_fixed_bam	VN:1.2-29-g0acaf2d	CL:java $jvm_args -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference_fasta -o $intervals_file -known $known_indels_file(s) 
@PG	ID:gatk_target_interval_creator.1	PN:GenomeAnalysisTK	PP:sam_to_fixed_bam.1	VN:1.2-29-g0acaf2d	CL:java $jvm_args -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference_fasta -o $intervals_file -known $known_indels_file(s) 
@PG	ID:bam_realignment_around_known_indels	PN:GenomeAnalysisTK	PP:gatk_target_interval_creator	VN:1.2-29-g0acaf2d	CL:java $jvm_args -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference_fasta -I $bam_file -o $realigned_bam_file -targetIntervals $intervals_file -known $known_indels_file(s) -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing
@PG	ID:bam_realignment_around_known_indels.1	PN:GenomeAnalysisTK	PP:gatk_target_interval_creator.1	VN:1.2-29-g0acaf2d	CL:java $jvm_args -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference_fasta -I $bam_file -o $realigned_bam_file -targetIntervals $intervals_file -known $known_indels_file(s) -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing
@PG	ID:bam_count_covariates	PN:GenomeAnalysisTK	PP:bam_realignment_around_known_indels	VN:1.2-29-g0acaf2d	CL:java $jvm_args -jar GenomeAnalysisTK.jar -T CountCovariates -R $reference_fasta -I $bam_file -recalFile $bam_file.recal_data.csv -knownSites $known_sites_file(s) -l INFO -L '1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;X;Y;MT' -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate
@PG	ID:bam_count_covariates.1	PN:GenomeAnalysisTK	PP:bam_realignment_around_known_indels.1	VN:1.2-29-g0acaf2d	CL:java $jvm_args -jar GenomeAnalysisTK.jar -T CountCovariates -R $reference_fasta -I $bam_file -recalFile $bam_file.recal_data.csv -knownSites $known_sites_file(s) -l INFO -L '1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;X;Y;MT' -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate
@PG	ID:bam_recalibrate_quality_scores	PN:GenomeAnalysisTK	PP:bam_count_covariates	VN:1.2-29-g0acaf2d	CL:java $jvm_args -jar GenomeAnalysisTK.jar -T TableRecalibration -R $reference_fasta -recalFile $bam_file.recal_data.csv -I $bam_file -o $recalibrated_bam_file -l INFO -compress 0 --disable_bam_indexing
@PG	ID:bam_recalibrate_quality_scores.1	PN:GenomeAnalysisTK	PP:bam_count_covariates.1	VN:1.2-29-g0acaf2d	CL:java $jvm_args -jar GenomeAnalysisTK.jar -T TableRecalibration -R $reference_fasta -recalFile $bam_file.recal_data.csv -I $bam_file -o $recalibrated_bam_file -l INFO -compress 0 --disable_bam_indexing
@PG	ID:bam_calculate_bq	PN:samtools	PP:bam_recalibrate_quality_scores	VN:0.1.17 (r973:277)	CL:samtools calmd -Erb $bam_file $reference_fasta > $bq_bam_file
@PG	ID:bam_calculate_bq.1	PN:samtools	PP:bam_recalibrate_quality_scores.1	VN:0.1.17 (r973:277)	CL:samtools calmd -Erb $bam_file $reference_fasta > $bq_bam_file
@PG	ID:bam_merge	PN:picard	PP:bam_calculate_bq	VN:1.53	CL:java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam VALIDATION_STRINGENCY=SILENT
@PG	ID:bam_merge.1	PN:picard	PP:bam_calculate_bq.1	VN:1.53	CL:java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam VALIDATION_STRINGENCY=SILENT
@PG	ID:bam_mark_duplicates	PN:picard	PP:bam_merge	VN:1.53	CL:java $jvm_args -jar MarkDuplicates.jar INPUT=$bam_file OUTPUT=$markdup_bam_file ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT
@PG	ID:bam_mark_duplicates.1	PN:picard	PP:bam_merge.1	VN:1.53	CL:java $jvm_args -jar MarkDuplicates.jar INPUT=$bam_file OUTPUT=$markdup_bam_file ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT
@PG	ID:bam_merge.2	PN:picard	PP:bam_mark_duplicates	VN:1.53	CL:java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam VALIDATION_STRINGENCY=SILENT
@PG	ID:bam_merge.1.2	PN:picard	PP:bam_mark_duplicates.1	VN:1.53	CL:java $jvm_args -jar MergeSamFiles.jar INPUT=$bam_file(s) OUTPUT=$merged_bam VALIDATION_STRINGENCY=SILENT
@CO	$known_indels_file(s) = ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz
@CO	$known_indels_file(s) .= ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf.gz
@CO	$known_sites_file(s) = ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz
]");
    assert(isValid(another_valid_header));
}
