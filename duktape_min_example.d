/*
  1) clone duktape repo and `sudo make -f Makefile.sharedlibrary install` from it
  2) export LD_LIBRARY_PATH=/usr/local/lib
  2) rdmd -d -IBioD -L-lduktape duktape_min_example.d <bam file>
 */
import bio.bam.reader;
import utils.duktape; // extern(C) declarations from duktype.h

import std.stdio;

extern(C) duk_ret_t
read_sequence_length(duk_context ctx) {
  auto read = cast(BamRead*)duk_get_pointer(ctx, 0);
  duk_push_int(ctx, read.sequence_length);
  return 1;
}

void main(string[] args) {
  auto read = new BamReader(args[1]).reads.front;

  duk_context ctx = duk_create_heap_default();
  scope(exit) duk_destroy_heap(ctx);

  // register the function globally
  duk_push_c_function(ctx, &read_sequence_length, 1);
  duk_put_global_string(ctx, "sequenceLength");

  duk_eval_string(ctx, "(function(r) {
       return sequenceLength(r) > 40;
  })");                         // stack: [func]

  duk_push_pointer(ctx, &read); // stack: [func bamread]
  duk_call(ctx, 1);             // stack: [result]

  duk_bool_t result = duk_get_boolean(ctx, -1);
  duk_pop(ctx);

  writeln("filter result:", to!bool(result));
}
