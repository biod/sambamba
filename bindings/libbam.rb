require 'ffi'

require './bindings/scaffolds/SqLine.rb'
require './bindings/scaffolds/RgLine.rb'
require './bindings/scaffolds/PgLine.rb'
require './bindings/scaffolds/SamHeader.rb'
require './bindings/scaffolds/ReferenceSequenceInfo.rb'
require './bindings/SamHeader.rb'
require './bindings/Alignment.rb'

class DArray < FFI::Struct 
	layout :length, :size_t,
		   :ptr, :pointer
end

module LibBAM 
	extend FFI::Library

	ffi_lib './libbam.so'
    
    # TODO: add typedefs so that bindings become more clear and readable

	attach_function :libbam_init, [], :void

    attach_function :get_last_error, [], :string

	attach_function :bamfile_new, [:string], :pointer
	attach_function :bamfile_destroy, [:pointer], :void
	attach_function :bamfile_get_header, [:pointer], SamHeader.by_value
	attach_function :bamfile_get_reference_sequences, [:pointer], DArray.by_value

    attach_function :bamfile_get_alignments, [:pointer], :pointer
    attach_function :alignment_range_destroy, [:pointer], :void
    attach_function :alignment_range_front, [:pointer], :pointer
    attach_function :alignment_range_popfront, [:pointer], :void
    attach_function :alignment_range_empty, [:pointer], :bool

    attach_function :alignment_destroy, [:pointer], :void
    attach_function :alignment_ref_id, [:pointer], :int32
    attach_function :alignment_position, [:pointer], :int32
    attach_function :alignment_bin, [:pointer], :uint16
    attach_function :alignment_mapping_quality, [:pointer], :uint8
    attach_function :alignment_flag, [:pointer], :uint16
    attach_function :alignment_sequence_length, [:pointer], :int32
    attach_function :alignment_next_ref_id, [:pointer], :int32
    attach_function :alignment_next_pos, [:pointer], :int32
    attach_function :alignment_template_length, [:pointer], :int32
    attach_function :alignment_read_name, [:pointer], DArray.by_value
    attach_function :alignment_cigar_string, [:pointer], :string
    attach_function :alignment_sequence, [:pointer], :string
    attach_function :alignment_phred_base_quality, [:pointer], DArray.by_value

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

    def alignments
        ptr = LibBAM.bamfile_get_alignments @ptr
        puts ptr
        AlignmentIterator.new(ptr, @ptr)
    end

    def self.finalize ptr
        proc { LibBAM.bamfile_destroy ptr }
    end
end

class AlignmentIterator
    include Enumerable

    def initialize(ptr, bam_ptr)
        @bam_ptr = bam_ptr # so that parent object won't get destroyed
                           # during lifetime of the iterator
        @ptr = ptr
        ObjectSpace.define_finalizer @ptr, AlignmentIterator.finalize(@ptr)
    end

    def each
        while not LibBAM.alignment_range_empty @ptr do
            yield Alignment.new(LibBAM.alignment_range_front @ptr)
            LibBAM.alignment_range_popfront @ptr
        end
    end

    def self.finalize ptr
        proc { LibBAM.alignment_range_destroy ptr }
    end
end
