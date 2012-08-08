#!/usr/bin/env ruby
require 'haml'

class Configuration
    def initialize(params)
        @params = params
    end

    def method_missing(name, *args)
        @params[name]
    end

    def to_s
        s = cpu + ' @ ' + freq + " (#{cores}"
        s += ' with hyperthreading' if hyperthreading
        s += ')'

        s + ', ' + memory
    end
end

Result = Struct.new :configuration, 
                    :time_elapsed, 
                    :peak_memory_usage,
                    :average_cpu_load

class CommandLine
    attr_reader :command_line
    attr_reader :results

    def initialize(command_line)
        @command_line = command_line
    end

    def add_results(*results)
        @results = *results
        self
    end
end

configurations = {
    :atom => Configuration.new(cpu: 'Intel Atom N450',
                               memory: '1GB',
                               freq: '1.66GHz',
                               cores: 1,
                               hyperthreading: true,
                               storage: 'HDD'),

    :xeon => Configuration.new(cpu: '2x Intel Xeon E5310',
                               memory: '8GB',
                               freq: '1.60GHz',
                               cores: 8,
                               storage: 'HDD'),

    :i5 => Configuration.new(cpu: '4x Intel Core i5-2500K',
                             memory: '8GB', ########## is it correct??? ### 
                             freq: '3.30GHz',
                             cores: 16,
                             storage: 'SSD')
}

sambamba = {}
samtools = {}

sambamba[:index_not_cached] = CommandLine.new("sambamba index $FILENAME").add_results(
    Result.new(configurations[:atom], 12.29,   32,     147),
    Result.new(configurations[:xeon], 6.96,    32,     139)
)

samtools[:index_not_cached] = CommandLine.new("samtools index $FILENAME").add_results(
    Result.new(configurations[:atom], 13.60,   1.4,     92),
    Result.new(configurations[:xeon], 8.73,    1.4,     93)
)

sambamba[:simple_filter_bam_output_not_cached] = CommandLine.new(%!sambamba view -f bam $FILENAME 20:10,000,000-20,000,000 -F "mapping_quality >= 50" -o test.bam!).add_results(
    Result.new(configurations[:atom], 22.96,    90,     98),
    Result.new(configurations[:xeon],  5.24,    90,    250)
)

samtools[:simple_filter_bam_output_not_cached] = CommandLine.new(%!samtools view -b $FILENAME 20:10,000,000-20,000,000 -q50 -o test.bam!).add_results(
    Result.new(configurations[:atom], 23.16,    1.8,    96),
    Result.new(configurations[:xeon], 10.83,    1.8,    98)
)

sambamba[:index_cached] = CommandLine.new("sambamba index $FILENAME").add_results(
    Result.new(configurations[:atom], 9.43,     32,    188),
    Result.new(configurations[:xeon], 2.21,     32,    433)
)

samtools[:index_cached] = CommandLine.new("samtools index $FILENAME").add_results(
    Result.new(configurations[:atom], 12.08,    1.4,    99),
    Result.new(configurations[:xeon],  7.98,    1.4,   100)
)

sambamba[:complex_filter_count_cached] = CommandLine.new(%!sambamba view $FILENAME -c -F "[RG] == 'ERR016156' and proper_pair and first_of_pair and not duplicate" 20:1000000-3000000!).add_results(
    Result.new(configurations[:atom], 0.53,     50,    144),
    Result.new(configurations[:xeon], 0.20,     50,    208)
)

samtools[:complex_filter_count_cached] = CommandLine.new(%!samtools view $FILENAME -c -r 'ERR016156' -f66 -F1024 20:1000000-3000000!).add_results(
    Result.new(configurations[:atom], 0.42,     1.3,    99),
    Result.new(configurations[:xeon], 0.27,     1.5,    100)
)

Operation = Struct.new :description, :sambamba, :samtools

engine = Haml::Engine.new(File.read 'table_template.haml')
puts engine.render(Object.new, :operations => [
                     Operation.new(
                       'Indexing BAM file (empty file cache)',
                       sambamba[:index_not_cached],
                       samtools[:index_not_cached]
                     ),
                     Operation.new(
                       'Indexing BAM file (file fully cached into RAM)',
                       sambamba[:index_cached],
                       samtools[:index_cached]
                     ),
                     Operation.new(
                       'Filtering reads from a region, with BAM output (empty file cache)',
                       sambamba[:simple_filter_bam_output_not_cached],
                       samtools[:simple_filter_bam_output_not_cached]
                     ),
                     Operation.new(
                       'Counting reads from a region (file fully cached into RAM)',
                       sambamba[:complex_filter_count_cached],
                       samtools[:complex_filter_count_cached]
                     )
                   ])
