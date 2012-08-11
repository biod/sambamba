import bamfile;
import bamoutput;
import std.conv;
import std.parallelism;
import std.stream;
import std.datetime;
import std.stdio;

void main(string[] args) {
    if (args.length < 7) {
        writeln("usage: ./writebam <filename> <n_threads> <compression_level> <inbufsz> <outbufsz> <parbufsz> [<parwusz>]");
        return;
    }
    auto fn = args[1];
    auto n_threads = to!int(args[2]);
    auto compression_level = to!int(args[3]);
    auto input_buffer_size = to!int(args[4]);
    auto output_buffer_size = to!int(args[5]);
    auto parallel_map_buffer_size = to!int(args[6]);
    auto parallel_map_work_unit_size = args.length >= 8 ? to!int(args[7]) : size_t.max;

    auto task_pool = new TaskPool(n_threads);
    scope(exit) task_pool.finish();

    auto bam = BamFile(fn, task_pool);
    bam.setBufferSize(input_buffer_size);
    auto stream = new BufferedFile("compressed.bam", FileMode.OutNew, output_buffer_size);

    StopWatch sw;
    sw.start();

    writeBAM(stream, bam.header.text, bam.reference_sequences, bam.alignments,
                compression_level, task_pool, parallel_map_buffer_size,
                parallel_map_work_unit_size);

    stream.close();

    sw.stop();
    writeln(n_threads, "\t", compression_level, "\t", input_buffer_size, "\t", 
            output_buffer_size, "\t", parallel_map_buffer_size, "\t", 
            parallel_map_work_unit_size, "\t", sw.peek().nsecs);
}
