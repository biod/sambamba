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

import bio.bam.tagvalue;
void push(T)(duk_context ctx, T value) if (is(T == bio.bam.tagvalue.Value))
{
  if (value.is_nothing) ctx.duk_push_null();
  else if (value.is_character) utils.duktape.push(ctx, value.to!char);
  else if (value.is_float) utils.duktape.push(ctx, value.to!float);
  else if (value.is_integer) utils.duktape.push(ctx, value.to!long);
  else if (value.is_string) utils.duktape.push(ctx, value.to!string);
  else if (value.is_numeric_array) {
    ctx.duk_push_external_buffer();
    void[] buf = *cast(void[]*)(&value);
    uint type_size = bio.bam.tagvalue.charToSizeof(value.bam_typeid);
    ctx.duk_config_buffer(-1, buf.ptr, buf.length * type_size);
    duk_uint_t type_tag;
    switch (value.bam_typeid) {
      case 'c': type_tag = DUK_BUFOBJ_INT8ARRAY; break;
      case 'C': type_tag = DUK_BUFOBJ_UINT8ARRAY; break;
      case 's': type_tag = DUK_BUFOBJ_INT16ARRAY; break;
      case 'S': type_tag = DUK_BUFOBJ_UINT16ARRAY; break;
      case 'i': type_tag = DUK_BUFOBJ_INT32ARRAY; break;
      case 'I': type_tag = DUK_BUFOBJ_UINT32ARRAY; break;
      case 'f': type_tag = DUK_BUFOBJ_FLOAT32ARRAY; break;
      default: throw new Exception("unsupported tag type");
    }
    ctx.duk_push_buffer_object(-1, 0, buf.length * type_size, type_tag);
  }
}

extern(C) duk_ret_t
getField(alias field)(duk_context ctx) {
  auto read = cast(BamRead*)duk_require_pointer(ctx, 0);
  mixin(`utils.duktape.push(ctx, read.` ~ field ~ `);`);
  return 1;
}

extern(C) duk_ret_t
getTag(duk_context ctx) {
  auto read = cast(BamRead*)duk_require_pointer(ctx, 0);

  duk_size_t s_len;
  auto s_ptr = ctx.duk_require_lstring(1, &s_len);
  if (s_len != 2) {
    push(ctx, Value(null));
    return 1;
  }

  auto tag = cast(string)(s_ptr[0 .. s_len]);
  push(ctx, (*read)[tag]);
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

  duk_push_c_function(ctx, &getTag, 2);
  duk_put_global_string(ctx, "getTag");

  auto filter = toStringz("(function(r) {return " ~ args[2] ~ "})");
  duk_eval_string(ctx, filter); // stack: [func]

  duk_push_pointer(ctx, &read); // stack: [func bamread]
  duk_call(ctx, 1);             // stack: [result]

  duk_bool_t result = duk_get_boolean(ctx, -1);
  duk_pop(ctx);

  writeln("filter result:", to!bool(result));
}
