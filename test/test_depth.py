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

def randomIntervals(ref_seqs, n):
    keys = ref_seqs.keys()
    for i in range(n):
        ref = random.choice(keys)
        length = ref_seqs[ref]
        start = random.randint(1, length)
        end = start + random.randint(1, 1000)
        yield (ref, start, end)

def writeBed(intervals, bedfile):
    for ref, start, end in intervals:
        bedfile.write("{}\t{}\t{}\n".format(ref, start, end))
    bedfile.flush()

def depthRegionReport(bam_fn, intervals):
    with tempfile.NamedTemporaryFile("w+") as bed:
        writeBed(intervals, bed)
        devnull = open("/dev/null", "w")
        res = subprocess.check_output([sambamba, "depth", "region", "-t", "0",
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

bam_fn = "dm3PacBio_valid.bam"
refs = referenceSequences(bam_fn)
intervals = list(randomIntervals(refs, 100000))
report = depthRegionReport(bam_fn, intervals)
expected = expectedRegionReport(bam_fn, intervals)

n_correct = 0
n_wrong = 0
for expected_entry, entry in zip(expected, report):
    if entry != expected_entry:
        if n_wrong == 0:
            saveRegionReport(report, "failed_report.txt")
            saveRegionReport(expected, "expected_report.txt")
            print 'TEST FAILURE'
            print 'results saved to expected_report.txt and failed_report.txt'
        print "different results:"
        print "expected: ", expected_entry
        print "got: ", entry
        n_wrong += 1
        if n_wrong >= 10:
            print "10 or more errors detected, exiting"
            break
    else:
        n_correct += 1

if n_correct == len(expected):
    print "tests passed!"
