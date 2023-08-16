 project('msgpack-d', 'd',
    meson_version: '>=0.47',
    license: 'BSL-1.0',
    version: '1.0.4'
)

project_soversion = '1'

pkgc = import('pkgconfig')

#
# Sources
#
msgpack_src = [
    'src/msgpack/attribute.d',
    'src/msgpack/buffer.d',
    'src/msgpack/common.d',
    'src/msgpack/exception.d',
    'src/msgpack/package.d',
    'src/msgpack/packer.d',
    'src/msgpack/register.d',
    'src/msgpack/streaming_unpacker.d',
    'src/msgpack/unpacker.d',
    'src/msgpack/value.d',
]

src_dir = include_directories('src/')

#
# Targets
#
msgpack_lib = library('msgpack-d',
        [msgpack_src],
        include_directories: [src_dir],
        install: true,
        version: meson.project_version(),
        soversion: project_soversion,
)

pkgc.generate(name: 'msgpack-d',
              libraries: [msgpack_lib],
              subdirs: 'd/msgpack',
              version: meson.project_version(),
              description: 'Library for lexing and parsing D source code.'
)

# for use by others which embed this as subproject
msgpack_dep = declare_dependency(
    link_with: [msgpack_lib],
    include_directories: [src_dir]
)

#
# Tests
#
if meson.get_compiler('d').get_id() == 'llvm'
  extra_args = ['-main', '-link-defaultlib-shared']
else
  extra_args = ['-main']
endif

msgpack_test_exe = executable('test_msgpack',
    [msgpack_src],
    include_directories: [src_dir],
    d_unittest: true,
    link_args: extra_args,
)
test('test_msgpack', msgpack_test_exe)

#
# Install
#
install_subdir('src/msgpack/', install_dir: 'include/d/msgpack/')
