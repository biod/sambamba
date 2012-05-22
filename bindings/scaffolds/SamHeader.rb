class SamHeader < FFI::Struct
  layout :_header_len, :size_t,
         :_header_data, :pointer,
         :_sq_lines_len, :size_t,
         :_sq_lines_data, :pointer,
         :_rg_lines_len, :size_t,
         :_rg_lines_data, :pointer,
         :_pg_lines_len, :size_t,
         :_pg_lines_data, :pointer,
         :_fasta_urls_len, :size_t,
         :_fasta_urls_data, :pointer,
         :_format_version_len, :size_t,
         :_format_version_data, :pointer,
         :_sorting_order_len, :size_t,
         :_sorting_order_data, :pointer

  def _header
    return nil if self[:_header_data].address == 0
    self[:_header_data].read_string(self[:_header_len])
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

  def _fasta_urls
    return nil if self[:_fasta_urls_data].address == 0
    (0...self[:_fasta_urls_len]).map {|i|
      string.new(self[:_fasta_urls_data] + i * string.size)
    }
  end

  def _format_version
    return nil if self[:_format_version_data].address == 0
    self[:_format_version_data].read_string(self[:_format_version_len])
  end

  def _sorting_order
    return nil if self[:_sorting_order_data].address == 0
    self[:_sorting_order_data].read_string(self[:_sorting_order_len])
  end
end