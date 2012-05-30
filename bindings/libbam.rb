require 'ffi'

require './bindings/scaffolds/SqLine.rb'
require './bindings/scaffolds/RgLine.rb'
require './bindings/scaffolds/PgLine.rb'
require './bindings/scaffolds/SamHeader.rb'
require './bindings/scaffolds/ReferenceSequenceInfo.rb'
require './bindings/SamHeader.rb'
require './bindings/Alignment.rb'
require './bindings/TagValue.rb'
require './bindings/DTypes.rb'

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
    attach_function :bamfile_get_valid_alignments, [:pointer], :pointer
    attach_function :bamfile_rewind, [:pointer], :void

    attach_function :alignment_range_destroy, [:pointer], :void
    attach_function :alignment_range_front, [:pointer], :pointer
    attach_function :alignment_range_popfront, [:pointer], :int
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
    attach_function :alignment_get_tag_value, [:pointer, :string], TagValue
    attach_function :alignment_get_all_tags, [:pointer], DHash
    attach_function :alignment_is_valid, [:pointer], :bool

    attach_function :dhash_destroy, [:pointer], :void
    attach_function :tag_value_destroy, [:pointer], :void

    @@initialized ||= false
    LibBAM.libbam_init unless @@initialized
end

class BamFile
    def initialize(filename)
        @ptr = LibBAM.bamfile_new filename
        if @ptr.address.zero?
            raise Exception.new(LibBAM.get_last_error)
        else 
            ObjectSpace.define_finalizer self, BamFile.finalize(@ptr)
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
			ReferenceSequenceInfo.new(ptr + k * sz)
		}
	end

    def alignments
        AlignmentIterator.new(@ptr)
    end

    def rewind!
        LibBAM.bamfile_rewind @ptr
    end

    def self.finalize ptr
        proc { LibBAM.bamfile_destroy ptr }
    end
end

class AlignmentIterator
    include Enumerable

    def initialize(bam_ptr)
        @bam_ptr = bam_ptr # so that parent object won't get destroyed
                           # during lifetime of the iterator
    end

    def each(options={})
        if options[:valid] == true then
            @ptr = LibBAM.bamfile_get_valid_alignments @bam_ptr
        else
            @ptr = LibBAM.bamfile_get_alignments @bam_ptr
        end
        ObjectSpace.define_finalizer self, AlignmentIterator.finalize(@ptr)

        is_empty = LibBAM.alignment_range_empty(@ptr)
        while not is_empty do
            yield Alignment.new(LibBAM.alignment_range_front @ptr)
            res = LibBAM.alignment_range_popfront @ptr
            if res == 1 then
                is_empty = true
            end
            if res == -1 then
                raise LibBAM.get_last_error
            end
        end
    end

    def self.finalize ptr
        proc { LibBAM.alignment_range_destroy ptr }
    end
end
