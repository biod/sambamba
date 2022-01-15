# Sambamba can use meson with:
#   source guix-build
#   rm -rf build/ ; env D_LD=gold CC=gcc meson build
#   cd build/
#   env CC=gcc ninja
#
# See also INSTALL.md

project('Sambamba', ['c','d'],
    license : 'GPL-2.0',
    version : '0.8.2',
    # default_options : ['buildtype=debugoptimized']
)

# lang = 'd'
# dc = meson.get_compiler(lang)

source_root = meson.source_root()
build_root = meson.build_root()

# we need a C compiler for the static library search (for BioD) to work,
# as the D compiler abstraction in Meson doesn't have find_library() capabilities yet

cc = meson.get_compiler('c')

# if meson.get_compiler('d').get_id() != 'llvm'
#     error('We only support the LLVM D compiler at this time. Please compile with LDC.')
# endif

extra_dflags = [ '-I=../BioD/' ]
if get_option('optimize_strong')
    extra_dflags += ['-O3', '-release', '-enable-inlining', '-boundscheck=off']
endif

#
# Sources
#
undead_src = [
    'BioD/contrib/undead/cstream.d',
    'BioD/contrib/undead/doformat.d',
    'BioD/contrib/undead/internal/file.d',
    'BioD/contrib/undead/stream.d',
    'BioD/contrib/undead/utf.d'
]

biod_src = [
    'BioD/bio/core/base.d',
    'BioD/bio/core/bgzf/block.d',
    'BioD/bio/core/bgzf/chunk.d',
    'BioD/bio/core/bgzf/compress.d',
    'BioD/bio/core/bgzf/constants.d',
    'BioD/bio/core/bgzf/inputstream.d',
    'BioD/bio/core/bgzf/outputstream.d',
    'BioD/bio/core/bgzf/virtualoffset.d',
    'BioD/bio/core/call.d',
    'BioD/bio/core/decompress.d',
    'BioD/bio/core/genotype.d',
    'BioD/bio/core/kmer.d',
    'BioD/bio/core/region.d',
    'BioD/bio/core/sequence.d',
    'BioD/bio/core/tinymap.d',
    'BioD/bio/core/utils/algo.d',
    'BioD/bio/core/utils/bylinefast.d',
    'BioD/bio/core/utils/exception.d',
    'BioD/bio/core/utils/format.d',
    'BioD/bio/core/utils/memoize.d',
    'BioD/bio/core/utils/outbuffer.d',
    'BioD/bio/core/utils/range.d',
    'BioD/bio/core/utils/roundbuf.d',
    'BioD/bio/core/utils/stream.d',
    'BioD/bio/core/utils/switchendianness.d',
    'BioD/bio/core/utils/tmpfile.d',
    'BioD/bio/core/utils/zlib.d',
    'BioD/bio/std/experimental/hts/bam/header.d',
    'BioD/bio/std/experimental/hts/bam/reader.d',
    'BioD/bio/std/experimental/hts/bam/writer.d',
    'BioD/bio/std/experimental/hts/bgzf.d',
    'BioD/bio/std/experimental/hts/bgzf_writer.d',
    'BioD/bio/std/experimental/hts/constants.d',
    'BioD/bio/std/experimental/hts/hashing.d',
    'BioD/bio/std/experimental/hts/logger.d',
    'BioD/bio/std/experimental/hts/pileup.d',
    'BioD/bio/std/experimental/hts/reads.d',
    'BioD/bio/std/experimental/hts/unpack.d',
    'BioD/bio/std/file/fai.d',
    'BioD/bio/std/file/fasta.d',
    'BioD/bio/std/file/fastq.d',
    'BioD/bio/std/genotype/maf.d',
    'BioD/bio/std/genotype/snp.d',
    'BioD/bio/std/hts/bam/abstractreader.d',
    'BioD/bio/std/hts/bam/bai/bin.d',
    'BioD/bio/std/hts/bam/baifile.d',
    'BioD/bio/std/hts/bam/bai/indexing.d',
    'BioD/bio/std/hts/bam/baseinfo.d',
    'BioD/bio/std/hts/bam/cigar.d',
    'BioD/bio/std/hts/bam/constants.d',
    'BioD/bio/std/hts/bam/md/core.d',
    'BioD/bio/std/hts/bam/md/operation.d',
    'BioD/bio/std/hts/bam/md/parse.d',
    'BioD/bio/std/hts/bam/md/reconstruct.d',
    'BioD/bio/std/hts/bam/multireader.d',
    'BioD/bio/std/hts/bam/pileup.d',
    'BioD/bio/std/hts/bam/randomaccessmanager.d',
    'BioD/bio/std/hts/bam/read.d',
    'BioD/bio/std/hts/bam/reader.d',
    'BioD/bio/std/hts/bam/readrange.d',
    'BioD/bio/std/hts/bam/reference.d',
    'BioD/bio/std/hts/bam/referenceinfo.d',
    'BioD/bio/std/hts/bam/region.d',
    'BioD/bio/std/hts/bam/splitter.d',
    'BioD/bio/std/hts/bam/tagvalue.d',
    'BioD/bio/std/hts/bam/validation/alignment.d',
    'BioD/bio/std/hts/bam/validation/samheader.d',
    'BioD/bio/std/hts/bam/writer.d',
    'BioD/bio/std/hts/iontorrent/flowcall.d',
    'BioD/bio/std/hts/iontorrent/flowindex.d',
    'BioD/bio/std/hts/sam/header.d',
    'BioD/bio/std/hts/sam/reader.d',
    'BioD/bio/std/hts/sam/utils/fastrecordparser.d',
    'BioD/bio/std/hts/sam/utils/recordparser.d',
    'BioD/bio/std/hts/snpcallers/maq.d',
    'BioD/bio/std/hts/snpcallers/simple.d',
    'BioD/bio/std/hts/thirdparty/msgpack.d',
    'BioD/bio/std/hts/utils/array.d',
    'BioD/bio/std/hts/utils/graph.d',
    'BioD/bio/std/hts/utils/samheadermerger.d',
    'BioD/bio/std/hts/utils/value.d',
    'BioD/bio/std/maf/block.d',
    'BioD/bio/std/maf/parser.d',
    'BioD/bio/std/maf/reader.d',
    'BioD/bio/std/range/splitter.d',
    'BioD/bio/std/sff/constants.d',
    'BioD/bio/std/sff/index.d',
    'BioD/bio/std/sff/read.d',
    'BioD/bio/std/sff/reader.d',
    'BioD/bio/std/sff/readrange.d',
    'BioD/bio/std/sff/utils/roundup.d',
    'BioD/bio/std/sff/writer.d'
]

sambamba_src = [
    'sambamba/depth.d',
    'sambamba/fixbins.d',
    'sambamba/flagstat.d',
    'sambamba/index.d',
    'sambamba/markdup2.d',
    'sambamba/markdup.d',
    'sambamba/merge.d',
    'sambamba/pileup.d',
    'sambamba/slice.d',
    'sambamba/sort.d',
    'sambamba/subsample.d',
    'sambamba/utils/common/bed.d',
    'sambamba/utils/common/file.d',
    'sambamba/utils/common/filtering.d',
    'sambamba/utils/common/intervaltree.d',
    'sambamba/utils/common/ldc_gc_workaround.d',
    'sambamba/utils/common/overwrite.d',
    'sambamba/utils/common/pratt_parser.d',
    'sambamba/utils/common/progressbar.d',
    'sambamba/utils/common/queryparser.d',
    'sambamba/utils/common/readstorage.d',
    'sambamba/utils/common/tmpdir.d',
    'sambamba/utils/view/alignmentrangeprocessor.d',
    'sambamba/utils/view/headerserializer.d',
    'sambamba/validate.d',
    'sambamba/view.d',
]

utils_src = [
    'utils/lz4.d',
    'utils/strip_bcf_header.d',
    'utils/version_.d'
]

cram_src = [
]

thirdparty_src = [
    'thirdparty/mergesort.d',
    'thirdparty/unstablesort.d'
]

manpages = [
    'man/sambamba.1',
    'man/sambamba-flagstat.1',
    'man/sambamba-index.1',
    'man/sambamba-markdup.1',
    'man/sambamba-merge.1',
    'man/sambamba-pileup.1',
    'man/sambamba-slice.1',
    'man/sambamba-sort.1',
    'man/sambamba-view.1'
]

#
# Dependencies
#
lz4     = dependency('liblz4')
zlib    = dependency('zlib')

#
# Configure
#

# Write LDC version to file

ldmd_prog = find_program('ldc2')
mkdir_prog = find_program('mkdir')

r = run_command(mkdir_prog.path(), '-p', build_root + '/utils/')
if r.returncode() != 0
  error('Unable to create "utils/" directory in build root: ' + r.stderr().strip())
endif

version_info_d_fname = build_root + '/utils/ldc_version_info_.d'

r = run_command('python3', source_root + '/gen_ldc_version_info.py', ldmd_prog.path(), version_info_d_fname)
if r.returncode() != 0
  error('Unable to write LDC version file: ' + r.stderr().strip())
endif

#
# Targets
#
sambamba_exe = executable('sambamba',
    ['sambamba/main.d',
     sambamba_src,
     biod_src,
     undead_src,
     utils_src,
     cram_src,
     thirdparty_src,
     version_info_d_fname],
    dependencies: [zlib,lz4],
    d_args: extra_dflags,
    d_import_dirs: [include_directories('.'),include_directories('BioD')],
    install: true
)

test_exe = find_program(join_paths(source_root, 'test', 'test_suite.sh'))
test('sambamba_test',test_exe,
     args: sambamba_exe.full_path(),
     workdir: source_root
)

#
# Install extra files
#
install_man(manpages)
install_data(['etc/bash_completion.d/sambamba'], install_dir: '/usr/share/bash-completion/completions')
