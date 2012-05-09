require 'ffi'

require './bindings/scaffolds/SqLine.rb'
require './bindings/scaffolds/RgLine.rb'
require './bindings/scaffolds/PgLine.rb'
require './bindings/scaffolds/SamHeader.rb'
require './bindings/scaffolds/ReferenceSequenceInfo.rb'
require './bindings/SamHeader.rb'

class DArray < FFI::Struct 
	layout :length, :size_t,
		   :ptr, :pointer
end

module LibBAM 
	extend FFI::Library

	ffi_lib './libbam.so'

	attach_function :libbam_init, [], :void

    attach_function :get_last_error, [], :string

	attach_function :bamfile_new, [:string], :pointer
	attach_function :bamfile_destroy, [:pointer], :void
	attach_function :bamfile_get_header, [:pointer], SamHeader.by_value
	attach_function :bamfile_get_reference_sequences, [:pointer], DArray.by_value

    @@initialized ||= false
    LibBAM.libbam_init unless @@initialized
end

class BamFile
    def initialize(filename)
        @ptr = LibBAM.bamfile_new filename
        if @ptr.address.zero?
            raise Exception.new(LibBAM.get_last_error)
        else 
            ObjectSpace.define_finalizer @ptr, BamFile.finalize(@ptr)
        end
    end

    def header
        LibBAM.bamfile_get_header @ptr
    end

	def reference_sequences
		return @ref_seqs unless @ref_seqs.nil?
		array = LibBAM.bamfile_get_reference_sequences @ptr
		ptr = array[:ptr]
		return nil if ptr.address.zero?

		sz = ReferenceSequenceInfo.size
		(0...array[:length]).map {|k|
			ReferenceSequenceInfo.new (ptr + k * sz)
		}
	end

    def self.finalize ptr
        proc { LibBAM.bamfile_destroy ptr }
    end
end
