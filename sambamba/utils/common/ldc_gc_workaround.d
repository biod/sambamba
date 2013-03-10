/// Remove address ranges related to system libraries, from GC scanning.
/// That makes a huge speed difference in case of frequent allocations.
module sambamba.utils.common.ldc_gc_workaround;

import std.stdio;
import std.format;
import std.algorithm;

import core.memory;

/// Shared static constructor does the job, no need to do any calls manually.
shared static this() {
    version(LDC) {
        removeExtraGcRanges();
    }
}

/// A path is excluded from GC scanning if it contains any of these substrings.
immutable defaultExclusionList = ["/libc-", "/ld-", "/libpthread-", "/librt-", "/libm-", "/libgcc_s", "/libdl-"];

/// On start, LDC parses /proc/self/maps file 
/// and adds every range with writable permissions to GC.
/// However, standard libraries don't contain any pointers to D objects,
/// and garbage collector spends a lot of time dealing with these ranges,
/// while this is completely unnecessary.
///
/// This function scans /proc/self/maps again and excludes ranges
/// that come from libc, ld, libpthread, librt, and libm.
/// You can pass your own exclusion list instead of the default one.
void removeExtraGcRanges(const char[][] exclusionList=defaultExclusionList) {

    debug {
        writeln("[DEBUG] removing extra GC ranges...");
    }

    version(Posix) {
        foreach (line; File("/proc/self/maps").byLine())
        {
            ulong start;
            ulong end;
            char[] permissions;
            ulong offset;
            char[] dev;
            ulong inode;
            char[] path;
            auto l = line.dup;

            formattedRead(l, "%x-%x %s %x %s %d %s", &start, &end, &permissions, 
                                                     &offset, &dev, &inode, &path);

            if (permissions[1] == 'w') {

                // exclude empty entries
                if (path.length == 0) {
                    GC.removeRange(cast(void*)start);
                }

                bool is_system_lib = false;
                foreach (substr; exclusionList) {
                    if (path.canFind(substr)) {

                        debug {
                            writeln("removing path ", path);
                        }

                        GC.removeRange(cast(void*)start);
                        is_system_lib = true;
                        break;
                    }
                } 
                
                debug {
                    if (!is_system_lib)
                        writeln("skipping path ", path);
                }
            }
        }
    }
}
