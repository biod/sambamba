class SamHeader
    alias_method :sq_lines, :_sq_lines
    alias_method :rg_lines, :_rg_lines
    alias_method :pg_lines, :_pg_lines
    undef_method :_sq_lines
    undef_method :_rg_lines
    undef_method :_pg_lines
end
