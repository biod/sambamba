require 'ffi'
require './bindings/TagValue.rb'

class Alignment

  def initialize(ptr)
    @ptr = ptr
    ObjectSpace.define_finalizer @ptr, Alignment.finalize(@ptr)
  end

  def self.finalize(ptr)
    proc { puts 'finalization...'; LibBAM.alignment_destroy ptr }
  end

  def [](tag)
    raise 'Tag length must be 2' if tag.length != 2 
    val_ptr = LibBAM.alignment_get_tag_value @ptr, tag
    val = TagValue.new val_ptr
    val.to_ruby_value
  end

  def tags
    dhash = DHash.new(LibBAM.alignment_get_all_tags @ptr)
    dhash.to_ruby_value
  end

  def ref_id
    LibBAM.alignment_ref_id @ptr
  end

  def position
    LibBAM.alignment_position @ptr
  end

  def bin
    LibBAM.alignment_bin @ptr
  end

  def mapping_quality
    LibBAM.alignment_mapping_quality @ptr
  end

  def flag
    @flag ||= LibBAM.alignment_flag @ptr
  end

  def sequence_length
    LibBAM.alignment_sequence_length @ptr
  end

  def next_ref_id
    LibBAM.alignment_next_ref_id @ptr
  end

  def next_pos
    LibBAM.alignment_next_pos @ptr
  end

  def template_length
    LibBAM.alignment_template_length @ptr
  end

  def read_name
    arr = LibBAM.alignment_read_name @ptr
    arr[:ptr].read_string(arr[:length])
  end

  def cigar_string
    LibBAM.alignment_cigar_string @ptr
  end

  def sequence
    LibBAM.alignment_sequence @ptr
  end

  def quality
    arr = LibBAM.alignment_phred_base_quality @ptr
    arr[:ptr].read_array_of_uint8(arr[:length])
  end

  def valid?
    LibBAM.alignment_is_valid @ptr
  end

  ### Template having multiple segments in sequencing
  def is_paired                
    (flag & 0x1) != 0
  end

  ### Each segment properly aligned according to the aligner
  def proper_pair              
    (flag & 0x2) != 0
  end

  ### Segment unmapped
  def is_unmapped              
    (flag & 0x4) != 0
  end

  ### Next segment in the template unmapped
  def mate_is_unmapped         
    (flag & 0x8) != 0
  end

  ### Sequence being reverse complemented
  def is_reverse_strand        
    (flag & 0x10) != 0
  end

  ### Sequence of the next segment in the template being reversed
  def mate_is_reverse_strand   
    (flag & 0x20) != 0
  end

  ### The first segment in the template
  def is_first_of_pair         
    (flag & 0x40) != 0
  end

  ### The last segment in the template
  def is_second_of_pair        
    (flag & 0x80) != 0
  end

  ### Secondary alignment
  def is_secondary_alignment   
    (flag & 0x100) != 0
  end

  ### Not passing quality controls
  def failed_quality_control   
    (flag & 0x200) != 0
  end

  ### PCR or optical duplicate
  def is_duplicate             
    (flag & 0x400) != 0
  end
end
