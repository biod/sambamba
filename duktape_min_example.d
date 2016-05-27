/*
  1) clone duktape repo and `sudo make -f Makefile.sharedlibrary install` from it
  2) export LD_LIBRARY_PATH=/usr/local/lib
  2) rdmd -d -IBioD -L-lduktape duktape_min_example.d <bam file> <JS expression>
 */
import bio.bam.reader;
import utils.duktape; // extern(C) declarations from duktype.h

import std.stdio;
import std.regex;
import std.string;
import std.meta : AliasSeq;

extern(C) duk_ret_t
getField(alias field)(duk_context ctx) {
  auto read = cast(BamRead*)duk_get_pointer(ctx, 0);
  mixin(`ctx.push(read.` ~ field ~ `);`);
  return 1;
}

void registerField(alias field)(duk_context ctx, string js_name=field) {
  js_name = js_name.replaceAll!(m => m.hit[1 .. $].toUpper)(regex(r"_\w"));

  duk_push_c_function(ctx, &getField!field, 1);
  duk_put_global_string(ctx, toStringz(js_name));
}

void main(string[] args) {
  auto read = new BamReader(args[1]).reads.front;

  duk_context ctx = duk_create_heap_default();
  scope(exit) duk_destroy_heap(ctx);

  // make each of these fields available as JS functions;
  // e.g. mate_ref_id can be accessed as mateRefId(r)
  foreach (field; AliasSeq!(
      "name", "ref_id", "ref_name", "position",
      "mapping_quality", "flag", "sequence_length",
      "mate_ref_id", "mate_ref_name", "mate_position",
      "template_length", "strand",
      "is_paired", "proper_pair",
      "is_unmapped", "mate_is_unmapped",
      "is_reverse_strand", "mate_is_reverse_strand",
      "is_first_of_pair", "is_second_of_pair",
      "is_secondary_alignment", "failed_quality_control",
      "is_duplicate", "is_supplementary",
      "basesCovered", "cigarString"))
      ctx.registerField!field();

  auto filter = toStringz("(function(r) {return " ~ args[2] ~ "})");
  duk_eval_string(ctx, filter); // stack: [func]

  duk_push_pointer(ctx, &read); // stack: [func bamread]
  duk_call(ctx, 1);             // stack: [result]

  duk_bool_t result = duk_get_boolean(ctx, -1);
  duk_pop(ctx);

  writeln("filter result:", to!bool(result));
}
