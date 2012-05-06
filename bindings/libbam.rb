require 'ffi'

require './bindings/scaffolds/SqLine.rb'
require './bindings/scaffolds/RgLine.rb'
require './bindings/scaffolds/PgLine.rb'
require './bindings/scaffolds/SamHeader.rb'
require './bindings/SamHeader.rb'

module LibBAM 
	extend FFI::Library

	ffi_lib './libbam.so'

	attach_function :libbam_init, [], :void

	attach_function :bamfile_new, [:string], :pointer
	attach_function :bamfile_destroy, [:pointer], :void
	attach_function :bamfile_get_header, [:pointer], SamHeader.by_value

    @@initialized ||= false
    LibBAM.libbam_init unless @@initialized
end

class BamFile
    def initialize(filename)
        @ptr = LibBAM.bamfile_new filename
        ObjectSpace.define_finalizer @ptr, BamFile.finalize(@ptr)
    end

    def header
        LibBAM.bamfile_get_header @ptr
    end

    def self.finalize ptr
        proc { LibBAM.bamfile_destroy ptr }
    end
end
