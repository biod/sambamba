class SamHeader
    alias_method :sq_lines, :_sq_lines
    alias_method :rg_lines, :_rg_lines
    alias_method :pg_lines, :_pg_lines
    undef_method :_sq_lines
    undef_method :_rg_lines
    undef_method :_pg_lines

    alias_method :raw_contents, :_header
    undef_method :_header

    alias_method :version, :_format_version
    undef_method :_format_version

    alias_method :sorting_order, :_sorting_order
    undef_method :_sorting_order

    #
    # Returns representation of @PG lines in dot format
    # so that the user can draw the graph with graphviz
    #
    def program_graph
        lines = []
        ids = {}

        self.pg_lines.each_with_index do |pg, i|
            ids[pg.identifier] = ["pg#{i}", pg.command_line || pg.program_name]
        end

        ids.values.each do |dot_id, command_line|
            lines << "#{dot_id} [label = \"#{command_line}\"]"
        end

        self.pg_lines.each do |pg|
            if pg.previous_program != nil
                lines << "#{ids[pg.previous_program].first} -> #{ids[pg.identifier].first}"
            end
        end

        "digraph G {\n\t" + lines.join(";\n\t") + "\n}"
    end
end
