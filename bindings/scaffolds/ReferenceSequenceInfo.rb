class ReferenceSequenceInfo < FFI::Struct
  layout :name_len, :size_t,
         :name_data, :pointer,
         :length, :int

  def name
    return nil if self[:name_data].address == 0
    self[:name_data].read_string(self[:name_len])
  end

  def length
    self[:length]
  end
end