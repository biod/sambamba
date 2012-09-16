import utils.msgpack;
import bamfile;

import std.stdio;
import std.array;

// Write alignments from a given file to stdout in MessagePack format
//
// Example of reading them from Ruby:
//
//   require 'msgpack'
//   unpacker = MessagePack::Unpacker.new STDIN
//   x = 0
//   begin
//     unpacker.each do |read|
//       x += 1
//       if (x & 0xFFF) == 0 then
//         ObjectSpace.garbage_collect
//       end
//     end
//   rescue EOFError
//   end
//
void main(string[] args) {
    auto bam = BamFile(args[1]);
    stdout.setvbuf(1_024_576);

    auto packer = packer(Appender!(ubyte[])());
    foreach (read; bam.alignments) {
        packer.pack(read);
        stdout.rawWrite(packer.stream.data);
        packer.stream.clear();
    }

    stdout.flush();
}
