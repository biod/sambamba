#!/usr/bin/env ruby
FILENAME = 'test.bam'

`echo 'nthreads\tcompression_level\tinput_bufsize\toutput_bufsize\tmap_bufsize\tmap_workunitsize\ttime' > stats.dat`

[1, 0, 2, 3].each do |compression_level|
    [12, 14, 16, 18].reverse.map {|n| 1<<n }.each do |input_buffer_size|
        [12, 16, 20, 24].reverse.map {|l| 1<<l }.each do |output_buffer_size|
            [2, 4, 6, 8].reverse.each do |nthreads|
                [3, 4, 5, 6].reverse.map {|k| [1<<k, k] }.each do |map_buf_size, k|
                    (0..k).map {|m| 1<<m }.each do |workunitsize|
                        `./writebam #{FILENAME} #{nthreads} #{compression_level} #{input_buffer_size} #{output_buffer_size} #{map_buf_size} #{workunitsize} >> stats.dat`
                    end
                end
            end
        end
    end
end
