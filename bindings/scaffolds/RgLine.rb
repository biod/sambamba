class RgLine < FFI::Struct
  layout :identifier_len, :size_t,
         :identifier_data, :pointer,
         :sequencing_center_len, :size_t,
         :sequencing_center_data, :pointer,
         :description_len, :size_t,
         :description_data, :pointer,
         :date_len, :size_t,
         :date_data, :pointer,
         :flow_order_len, :size_t,
         :flow_order_data, :pointer,
         :key_sequence_len, :size_t,
         :key_sequence_data, :pointer,
         :library_len, :size_t,
         :library_data, :pointer,
         :programs_len, :size_t,
         :programs_data, :pointer,
         :predicted_insert_size, :uint,
         :platform_len, :size_t,
         :platform_data, :pointer,
         :platform_unit_len, :size_t,
         :platform_unit_data, :pointer,
         :sample_len, :size_t,
         :sample_data, :pointer

  def identifier
    return nil if self[:identifier_data].address == 0
    self[:identifier_data].read_string(self[:identifier_len])
  end

  def sequencing_center
    return nil if self[:sequencing_center_data].address == 0
    self[:sequencing_center_data].read_string(self[:sequencing_center_len])
  end

  def description
    return nil if self[:description_data].address == 0
    self[:description_data].read_string(self[:description_len])
  end

  def date
    return nil if self[:date_data].address == 0
    self[:date_data].read_string(self[:date_len])
  end

  def flow_order
    return nil if self[:flow_order_data].address == 0
    self[:flow_order_data].read_string(self[:flow_order_len])
  end

  def key_sequence
    return nil if self[:key_sequence_data].address == 0
    self[:key_sequence_data].read_string(self[:key_sequence_len])
  end

  def library
    return nil if self[:library_data].address == 0
    self[:library_data].read_string(self[:library_len])
  end

  def programs
    return nil if self[:programs_data].address == 0
    self[:programs_data].read_string(self[:programs_len])
  end

  def predicted_insert_size
    self[:predicted_insert_size]
  end

  def platform
    return nil if self[:platform_data].address == 0
    self[:platform_data].read_string(self[:platform_len])
  end

  def platform_unit
    return nil if self[:platform_unit_data].address == 0
    self[:platform_unit_data].read_string(self[:platform_unit_len])
  end

  def sample
    return nil if self[:sample_data].address == 0
    self[:sample_data].read_string(self[:sample_len])
  end
end