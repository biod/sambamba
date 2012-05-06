class PgLine < FFI::Struct
  layout :identifier_len, :size_t,
         :identifier_data, :pointer,
         :program_name_len, :size_t,
         :program_name_data, :pointer,
         :command_line_len, :size_t,
         :command_line_data, :pointer,
         :previous_program_len, :size_t,
         :previous_program_data, :pointer,
         :program_version_len, :size_t,
         :program_version_data, :pointer

  def identifier
    return nil if self[:identifier_data].address == 0
    self[:identifier_data].read_string(self[:identifier_len])
  end

  def program_name
    return nil if self[:program_name_data].address == 0
    self[:program_name_data].read_string(self[:program_name_len])
  end

  def command_line
    return nil if self[:command_line_data].address == 0
    self[:command_line_data].read_string(self[:command_line_len])
  end

  def previous_program
    return nil if self[:previous_program_data].address == 0
    self[:previous_program_data].read_string(self[:previous_program_len])
  end

  def program_version
    return nil if self[:program_version_data].address == 0
    self[:program_version_data].read_string(self[:program_version_len])
  end
end