#!/usr/bin/env python3

import msgpack

import subprocess
import os.path
import random
import tempfile
import csv

sambamba = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../build/sambamba')

def referenceSequences(bam_fn):
    header = subprocess.check_output([sambamba, "view", "-H", bam_fn,
                                      "-f", "msgpack"])
    return {item['SN']: item['LN'] for item in msgpack.unpackb(header)[2]}

def randomIntervals(ref_seqs, n, max_length=1000):
    keys = ref_seqs.keys()
    result = []
    for i in range(n):
        overlap = i > 0 and random.randint(1, 100) <= 20
        if overlap:
            ref, start, e = random.choice(result)
            start += random.randint(0, (e - start - 1))
        else:
            ref = random.choice(keys)
            length = ref_seqs[ref]
            start = random.randint(1, length)
        end = start + random.randint(1, max_length)
        result.append((ref, start, end))
    return result

def writeBed(intervals, bedfile):
    for ref, start, end in intervals:
        bedfile.write("{}\t{}\t{}\n".format(ref, start, end))
    bedfile.flush()

def depthRegionReport(bam_fn, intervals):
    with tempfile.NamedTemporaryFile("w+") as bed:
        writeBed(intervals, bed)
        devnull = open("/dev/null", "w")
        res = subprocess.check_output([sambamba, "depth", "region", "-t", "0",
                                       "--combined",
                                       bam_fn, "-L", bed.name, "-F", ""],
                                      stderr=devnull)
        lst = list(csv.DictReader(res[2:].split("\n"), dialect='excel-tab'))
        return [(i['chrom'], int(i['chromStart']), int(i['chromEnd']),
                 i['readCount']) for i in lst]

def expectedRegionReport(bam_fn, intervals):
    report = []
    for ref, start, end in intervals:
        # coordinate translation: half-open zero-based -> closed 1-based
        cov = subprocess.check_output([sambamba, "view", bam_fn, "-c",
                                       "-t", "0",
                                       "{}:{}-{}".format(ref, start+1, end)])
        report.append((ref, start, end, cov.strip()))
    return report

def saveRegionReport(report, fn):
    with open(fn, "w+") as f:
        for fields in report:
            f.write("{}\t{}\t{}\t{}\n".format(*fields))

def depthBaseReport(bam_fn, intervals):
    with tempfile.NamedTemporaryFile("w+") as bed:
        writeBed(intervals, bed)
        devnull = open("/dev/null", "w")
        res = subprocess.check_output([sambamba, "depth", "base", "-t", "0",
                                       "--combined",
                                       bam_fn, "-L", bed.name, "-F", "", "-z"],
                                       stderr=devnull)
        lst = list(csv.DictReader(res.split("\n"), dialect='excel-tab'))
        return [(i['REF'], int(i['POS']), int(i['COV'])) for i in lst]

bam_fn = "dm3PacBio_valid.bam"

# minor differences are still possible e.g. when region is wholly inside a deletion,
# but they are quite rare and can be checked manually
def test_depth_region():
    refs = referenceSequences(bam_fn)

    intervals = list(randomIntervals(refs, 10000, 1000))
    report = depthRegionReport(bam_fn, intervals)
    expected = expectedRegionReport(bam_fn, intervals)

    return compareResults(report, expected, "region_")

def compareResults(report, expected, prefix):
    n_correct = 0
    n_wrong = 0
    for expected_entry, entry in zip(expected, report):
        if entry != expected_entry:
            if n_wrong == 0:
                fn_failed = prefix + "failed_report.txt"
                fn_expected = prefix + "expected_report.txt"
                saveRegionReport(report, fn_failed)
                saveRegionReport(expected, fn_expected)
                print 'TEST FAILURE'
                print "results saved to " + fn_expected + " and " + fn_failed
            print "different results:"
            print "expected: ", expected_entry
            print "got: ", entry
            n_wrong += 1
            if n_wrong >= 10:
                print "10 or more errors detected, exiting"
                break
        else:
            n_correct += 1

    return n_correct == len(expected)

def test_depth_base():
    refs = referenceSequences(bam_fn)
    intervals = list(randomIntervals(refs, 1000, 100))
    report = depthBaseReport(bam_fn, intervals)
    expected = []
    for ref, pos, count in report:
        cov = subprocess.check_output([sambamba, "view", bam_fn, "-c",
                                       "-t", "0",
                                       "{}:{}-{}".format(ref, pos+1, pos+1)])
        cov = int(cov.strip())
        expected.append((ref, pos, cov))

    return compareResults(report, expected, "base_")

test_depth_base()
