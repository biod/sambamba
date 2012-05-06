class SamHeader < FFI::Struct
  layout :header_len, :size_t,
         :header_data, :pointer,
         :_sq_lines_len, :size_t,
         :_sq_lines_data, :pointer,
         :_rg_lines_len, :size_t,
         :_rg_lines_data, :pointer,
         :_pg_lines_len, :size_t,
         :_pg_lines_data, :pointer

  def header
    return nil if self[:header_data].address == 0
    self[:header_data].read_string(self[:header_len])
  end

  def _sq_lines
    return nil if self[:_sq_lines_data].address == 0
    (0...self[:_sq_lines_len]).map {|i|
      SqLine.new(self[:_sq_lines_data] + i * SqLine.size)
    }
  end

  def _rg_lines
    return nil if self[:_rg_lines_data].address == 0
    (0...self[:_rg_lines_len]).map {|i|
      RgLine.new(self[:_rg_lines_data] + i * RgLine.size)
    }
  end

  def _pg_lines
    return nil if self[:_pg_lines_data].address == 0
    (0...self[:_pg_lines_len]).map {|i|
      PgLine.new(self[:_pg_lines_data] + i * PgLine.size)
    }
  end
end