class SqLine < FFI::Struct
  layout :sequence_name_len, :size_t,
         :sequence_name_data, :pointer,
         :sequence_length, :uint,
         :assembly_len, :size_t,
         :assembly_data, :pointer,
         :md5_len, :size_t,
         :md5_data, :pointer,
         :species_len, :size_t,
         :species_data, :pointer,
         :uri_len, :size_t,
         :uri_data, :pointer

  def sequence_name
    return nil if self[:sequence_name_data].address == 0
    self[:sequence_name_data].read_string(self[:sequence_name_len])
  end

  def sequence_length
    self[:sequence_length]
  end

  def assembly
    return nil if self[:assembly_data].address == 0
    self[:assembly_data].read_string(self[:assembly_len])
  end

  def md5
    return nil if self[:md5_data].address == 0
    self[:md5_data].read_string(self[:md5_len])
  end

  def species
    return nil if self[:species_data].address == 0
    self[:species_data].read_string(self[:species_len])
  end

  def uri
    return nil if self[:uri_data].address == 0
    self[:uri_data].read_string(self[:uri_len])
  end
end