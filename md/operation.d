module md.operation;

import std.conv;

/// MD tag operation types
enum MdOperationType : ubyte {
    Match,
    Mismatch,
    Deletion
}

/// Single MD operation.
struct MdOperation {
    MdOperationType type; /// Operation type
    union {
        uint match; /// If this is a match, contains the number of matched bases.
        string deletion; /// If this is deletion, contains deleted sequence.
        char mismatch; /// If this is a mismatch, contains the mismatched reference base.
    }

    static MdOperation createMatch(uint match) {
        MdOperation m;
        m.type = MdOperationType.Match;
        m.match = match;
        return m;
    }

    static MdOperation createDeletion(string deletion) {
        MdOperation m;
        m.type = MdOperationType.Deletion;
        m.deletion = deletion;
        return m;
    }

    static MdOperation createMismatch(char mismatch) {
        MdOperation m;
        m.type = MdOperationType.Mismatch;
        m.mismatch = mismatch;
        return m;
    }

    bool opEquals(ref const MdOperation other) const {

        if (type != other.type) {
            return false;
        }

        final switch (type) {
            case MdOperationType.Match:
                return match == other.match;
            case MdOperationType.Mismatch:
                return mismatch == other.mismatch;
            case MdOperationType.Deletion:
                return deletion == other.deletion;
        }
    }

    string toString() const {
        final switch (type) {
            case MdOperationType.Match:
                return "Match(" ~ to!string(match) ~ ")";
            case MdOperationType.Mismatch:
                return "Mismatch(" ~ to!string(mismatch) ~ ")";
            case MdOperationType.Deletion:
                return "Deletion(" ~ to!string(deletion) ~ ")";
        }
    }
}
