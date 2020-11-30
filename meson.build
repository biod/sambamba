project('Sambamba', 'd',
    meson_version : '>=0.48',
    license : 'GPL-2.0',
    version : '0.7.1',
    default_options : ['buildtype=debugoptimized']
)

source_root = meson.source_root()
build_root = meson.build_root()

# we need a C compiler for the static library search (for BioD) to work,
# as the D compiler abstraction in Meson doesn't have find_library() capabilities yet
add_languages('c')

if meson.get_compiler('d').get_id() != 'llvm'
    error('We only support the LLVM D compiler at time. Please compile with LDC.')
endif
extra_dflags = []
if get_option('optimize_strong')
    extra_dflags = ['-O3', '-release', '-enable-inlining', '-boundscheck=off']
endif

#
# Sources
#
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
biod_dep   = dependency('biod', version: '>=0.1.0', static: true)
lz4_dep    = dependency('liblz4')

#
# Configure
#

# Write LDC version to file
ldmd_prog = find_program('ldmd2')
mkdir_prog = find_program('mkdir')
r = run_command(mkdir_prog.path(), '-p', build_root + '/utils/')
if r.returncode() != 0
  error('Unable to create "utils/" directory in build root: ' + r.stderr().strip())
endif
version_info_d_fname = build_root + '/utils/ldc_version_info_.d'
r = run_command('sh', '-c', source_root + '/gen_ldc_version_info.py ' + ldmd_prog.path() + ' > ' + version_info_d_fname)
if r.returncode() != 0
  error('Unable to write LDC version file: ' + r.stderr().strip())
endif

#
# Targets
#
sambamba_exe = executable('sambamba',
    ['sambamba/main.d',
     sambamba_src,
     utils_src,
     cram_src,
     thirdparty_src,
     version_info_d_fname],
    dependencies: [biod_dep,
                   lz4_dep,
    d_args: extra_dflags,
    d_import_dirs: [include_directories('.')],
    install: true
)

test_exe = find_program(join_paths(source_root, 'test', 'test_suite.sh'))
test('sambamba_test',
     test_exe,
     env: ['sambamba=' + sambamba_exe.full_path()],
     workdir: source_root
)

#
# Install extra files
#
install_man(manpages)
install_data(['etc/bash_completion.d/sambamba'], install_dir: '/usr/share/bash-completion/completions')
