/*
  most of the file scaffolded with dstep on ubuntu with llvm 3.9 installed:
  ./dstep -c -Iduktape-1.5.0/src -I/usr/lib/llvm-3.9/lib/clang/3.9.0/include -o duktape.d duktape-1.5.0/src/duktape.h
*/
module utils.duktape;

import core.stdc.config;
import core.stdc.stdio;

extern (C) {
immutable DUK_COMPILE_EVAL =                   (1 << 3)    /* compile eval code (instead of global code) */;
immutable DUK_COMPILE_FUNCTION =               (1 << 4)    /* compile function code (instead of global code) */;
immutable DUK_COMPILE_STRICT =                 (1 << 5)    /* use strict (outer) context for global, eval, or function code */;
immutable DUK_COMPILE_SAFE =                   (1 << 6)    /* (internal) catch compilation errors */;
immutable DUK_COMPILE_NORESULT =               (1 << 7)    /* (internal) omit eval result */;
immutable DUK_COMPILE_NOSOURCE =               (1 << 8)    /* (internal) no source string on stack */;
immutable DUK_COMPILE_STRLEN =                 (1 << 9)    /* (internal) take strlen() of src_buffer (avoids double evaluation in macro) */;
immutable DUK_COMPILE_NOFILENAME =             (1 << 10)    /* (internal) no filename on stack */;

/* Flags for duk_def_prop() and its variants */
immutable DUK_DEFPROP_WRITABLE =               (1 << 0)    /* set writable (effective if DUK_DEFPROP_HAVE_WRITABLE set) */;
immutable DUK_DEFPROP_ENUMERABLE =             (1 << 1)    /* set enumerable (effective if DUK_DEFPROP_HAVE_ENUMERABLE set) */;
immutable DUK_DEFPROP_CONFIGURABLE =           (1 << 2)    /* set configurable (effective if DUK_DEFPROP_HAVE_CONFIGURABLE set) */;
immutable DUK_DEFPROP_HAVE_WRITABLE =         (1 << 3)    /* set/clear writable */;
immutable DUK_DEFPROP_HAVE_ENUMERABLE =       (1 << 4)    /* set/clear enumerable */;
immutable DUK_DEFPROP_HAVE_CONFIGURABLE =     (1 << 5)    /* set/clear configurable */;
immutable DUK_DEFPROP_HAVE_VALUE =            (1 << 6)    /* set value (given on value stack) */;
immutable DUK_DEFPROP_HAVE_GETTER =           (1 << 7)    /* set getter (given on value stack) */;
immutable DUK_DEFPROP_HAVE_SETTER =           (1 << 8)    /* set setter (given on value stack) */;
immutable DUK_DEFPROP_FORCE =                  (1 << 9)    /* force change if possible, may still fail for e.g. virtual properties */;
immutable DUK_DEFPROP_SET_WRITABLE =          (DUK_DEFPROP_HAVE_WRITABLE | DUK_DEFPROP_WRITABLE);
immutable DUK_DEFPROP_CLEAR_WRITABLE =        DUK_DEFPROP_HAVE_WRITABLE;
immutable DUK_DEFPROP_SET_ENUMERABLE =        (DUK_DEFPROP_HAVE_ENUMERABLE | DUK_DEFPROP_ENUMERABLE);
immutable DUK_DEFPROP_CLEAR_ENUMERABLE =      DUK_DEFPROP_HAVE_ENUMERABLE;
immutable DUK_DEFPROP_SET_CONFIGURABLE =      (DUK_DEFPROP_HAVE_CONFIGURABLE | DUK_DEFPROP_CONFIGURABLE);
immutable DUK_DEFPROP_CLEAR_CONFIGURABLE =    DUK_DEFPROP_HAVE_CONFIGURABLE;

/* Flags for duk_push_thread_raw() */
immutable DUK_THREAD_NEW_GLOBAL_ENV =       (1 << 0)    /* create a new global environment */;

/* Flags for duk_push_string_file_raw() */
immutable DUK_STRING_PUSH_SAFE =              (1 << 0)    /* no error if file does not exist */;

/* Duktape specific error codes (must be 8 bits at most, see duk_error.h) */
immutable DUK_ERR_NONE =                       0    /* no error (e.g. from duk_get_error_code()) */;
immutable DUK_ERR_UNIMPLEMENTED_ERROR =       50   /* UnimplementedError */  /* XXX: replace with TypeError? */;
immutable DUK_ERR_UNSUPPORTED_ERROR =         51   /* UnsupportedError */    /* XXX: replace with TypeError? */;
immutable DUK_ERR_INTERNAL_ERROR =            52   /* InternalError */       /* XXX: replace with plain Error? */;
immutable DUK_ERR_ALLOC_ERROR =               53   /* AllocError */          /* XXX: replace with RangeError? */;
immutable DUK_ERR_ASSERTION_ERROR =           54   /* AssertionError */      /* XXX: to be removed? */;
immutable DUK_ERR_API_ERROR =                 55   /* APIError */            /* XXX: replace with TypeError? */;
immutable DUK_ERR_UNCAUGHT_ERROR =            56   /* UncaughtError */       /* XXX: to be removed? */;

/* Ecmascript E5 specification error codes */
immutable DUK_ERR_ERROR =                      100  /* Error */;
immutable DUK_ERR_EVAL_ERROR =                101  /* EvalError */;
immutable DUK_ERR_RANGE_ERROR =               102  /* RangeError */;
immutable DUK_ERR_REFERENCE_ERROR =           103  /* ReferenceError */;
immutable DUK_ERR_SYNTAX_ERROR =              104  /* SyntaxError */;
immutable DUK_ERR_TYPE_ERROR =                105  /* TypeError */;
immutable DUK_ERR_URI_ERROR =                 106  /* URIError */;

/* Return codes for C functions (shortcut for throwing an error) */
immutable DUK_RET_UNIMPLEMENTED_ERROR =       (-DUK_ERR_UNIMPLEMENTED_ERROR);
immutable DUK_RET_UNSUPPORTED_ERROR =         (-DUK_ERR_UNSUPPORTED_ERROR);
immutable DUK_RET_INTERNAL_ERROR =            (-DUK_ERR_INTERNAL_ERROR);
immutable DUK_RET_ALLOC_ERROR =               (-DUK_ERR_ALLOC_ERROR);
immutable DUK_RET_ASSERTION_ERROR =           (-DUK_ERR_ASSERTION_ERROR);
immutable DUK_RET_API_ERROR =                 (-DUK_ERR_API_ERROR);
immutable DUK_RET_UNCAUGHT_ERROR =            (-DUK_ERR_UNCAUGHT_ERROR);
immutable DUK_RET_ERROR =                      (-DUK_ERR_ERROR);
immutable DUK_RET_EVAL_ERROR =                (-DUK_ERR_EVAL_ERROR);
immutable DUK_RET_RANGE_ERROR =               (-DUK_ERR_RANGE_ERROR);
immutable DUK_RET_REFERENCE_ERROR =           (-DUK_ERR_REFERENCE_ERROR);
immutable DUK_RET_SYNTAX_ERROR =              (-DUK_ERR_SYNTAX_ERROR);
immutable DUK_RET_TYPE_ERROR =                (-DUK_ERR_TYPE_ERROR);
immutable DUK_RET_URI_ERROR =                 (-DUK_ERR_URI_ERROR);

/* Return codes for protected calls (duk_safe_call(), duk_pcall()) */
immutable DUK_EXEC_SUCCESS =                   0;
immutable DUK_EXEC_ERROR =                     1;

/* Log levels */
immutable DUK_LOG_TRACE =                      0;
immutable DUK_LOG_DEBUG =                      1;
immutable DUK_LOG_INFO =                       2;
immutable DUK_LOG_WARN =                       3;
immutable DUK_LOG_ERROR =                      4;
immutable DUK_LOG_FATAL =                      5;

alias void* duk_context;
alias int duk_idx_t;
alias int duk_int_t;
alias int duk_small_int_t;
alias int duk_small_uint_t;
alias uint duk_uint_t;
alias duk_uint_t duk_uarridx_t;
alias double duk_double_t;
alias ubyte duk_uint8_t;
alias ushort duk_uint16_t;
alias int duk_int32_t;
alias uint duk_uint32_t;
alias ulong duk_uint64_t;
alias duk_small_uint_t duk_bool_t;
alias size_t duk_size_t;
alias duk_int_t duk_codepoint_t;

alias duk_small_int_t duk_ret_t;
alias duk_int_t duk_errcode_t;

alias int function (duk_context) duk_c_function;
alias void* function (void*, c_ulong) duk_alloc_function;
alias void* function (void*, void*, c_ulong) duk_realloc_function;
alias void function (void*, void*) duk_free_function;
alias void function (duk_context, int, const(char)*) duk_fatal_function;
alias void function (void*, int) duk_decode_char_function;
alias int function (void*, int) duk_map_char_function;
alias int function (duk_context) duk_safe_call_function;
alias c_ulong function (void*, char*, c_ulong) duk_debug_read_function;
alias c_ulong function (void*, const(char)*, c_ulong) duk_debug_write_function;
alias c_ulong function (void*) duk_debug_peek_function;
alias void function (void*) duk_debug_read_flush_function;
alias void function (void*) duk_debug_write_flush_function;
alias int function (duk_context, void*, int) duk_debug_request_function;
alias void function (void*) duk_debug_detached_function;

struct duk_memory_functions
{
    duk_alloc_function alloc_func;
    duk_realloc_function realloc_func;
    duk_free_function free_func;
    void* udata;
}

struct duk_function_list_entry
{
    const(char)* key;
    duk_c_function value;
    duk_idx_t nargs;
}

struct duk_number_list_entry
{
    const(char)* key;
    duk_double_t value;
}

union duk_double_union
{
    double d;
    float[2] f;
    duk_uint64_t[1] ull;
    duk_uint32_t[2] ui;
    duk_uint16_t[4] us;
    duk_uint8_t[8] uc;
}

duk_context duk_create_heap (duk_alloc_function alloc_func, duk_realloc_function realloc_func, duk_free_function free_func, void* heap_udata, duk_fatal_function fatal_handler);
void duk_destroy_heap (duk_context ctx);
void* duk_alloc_raw (duk_context ctx, duk_size_t size);
void duk_free_raw (duk_context ctx, void* ptr);
void* duk_realloc_raw (duk_context ctx, void* ptr, duk_size_t size);
void* duk_alloc (duk_context ctx, duk_size_t size);
void duk_free (duk_context ctx, void* ptr);
void* duk_realloc (duk_context ctx, void* ptr, duk_size_t size);
void duk_get_memory_functions (duk_context ctx, duk_memory_functions* out_funcs);
void duk_gc (duk_context ctx, duk_uint_t flags);
void duk_throw (duk_context ctx);
void duk_fatal (duk_context ctx, duk_errcode_t err_code, const(char)* err_msg);
void duk_error_raw (duk_context ctx, duk_errcode_t err_code, const(char)* filename, duk_int_t line, const(char)* fmt, ...);
duk_bool_t duk_is_strict_call (duk_context ctx);
duk_bool_t duk_is_constructor_call (duk_context ctx);
duk_idx_t duk_normalize_index (duk_context ctx, duk_idx_t index);
duk_idx_t duk_require_normalize_index (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_valid_index (duk_context ctx, duk_idx_t index);
void duk_require_valid_index (duk_context ctx, duk_idx_t index);
duk_idx_t duk_get_top (duk_context ctx);
void duk_set_top (duk_context ctx, duk_idx_t index);
duk_idx_t duk_get_top_index (duk_context ctx);
duk_idx_t duk_require_top_index (duk_context ctx);
duk_bool_t duk_check_stack (duk_context ctx, duk_idx_t extra);
void duk_require_stack (duk_context ctx, duk_idx_t extra);
duk_bool_t duk_check_stack_top (duk_context ctx, duk_idx_t top);
void duk_require_stack_top (duk_context ctx, duk_idx_t top);
void duk_swap (duk_context ctx, duk_idx_t index1, duk_idx_t index2);
void duk_swap_top (duk_context ctx, duk_idx_t index);
void duk_dup (duk_context ctx, duk_idx_t from_index);
void duk_dup_top (duk_context ctx);
void duk_insert (duk_context ctx, duk_idx_t to_index);
void duk_replace (duk_context ctx, duk_idx_t to_index);
void duk_copy (duk_context ctx, duk_idx_t from_index, duk_idx_t to_index);
void duk_remove (duk_context ctx, duk_idx_t index);
void duk_xcopymove_raw (duk_context to_ctx, duk_context from_ctx, duk_idx_t count, duk_bool_t is_copy);
void duk_push_undefined (duk_context ctx);
void duk_push_null (duk_context ctx);
void duk_push_boolean (duk_context ctx, duk_bool_t val);
void duk_push_true (duk_context ctx);
void duk_push_false (duk_context ctx);
void duk_push_number (duk_context ctx, duk_double_t val);
void duk_push_nan (duk_context ctx);
void duk_push_int (duk_context ctx, duk_int_t val);
void duk_push_uint (duk_context ctx, duk_uint_t val);
const(char)* duk_push_string (duk_context ctx, const(char)* str);
const(char)* duk_push_lstring (duk_context ctx, const(char)* str, duk_size_t len);
void duk_push_pointer (duk_context ctx, void* p);
const(char)* duk_push_sprintf (duk_context ctx, const(char)* fmt, ...);
const(char)* duk_push_string_file_raw (duk_context ctx, const(char)* path, duk_uint_t flags);
void duk_push_this (duk_context ctx);
void duk_push_current_function (duk_context ctx);
void duk_push_current_thread (duk_context ctx);
void duk_push_global_object (duk_context ctx);
void duk_push_heap_stash (duk_context ctx);
void duk_push_global_stash (duk_context ctx);
void duk_push_thread_stash (duk_context ctx, duk_context target_ctx);
duk_idx_t duk_push_object (duk_context ctx);
duk_idx_t duk_push_array (duk_context ctx);
duk_idx_t duk_push_c_function (duk_context ctx, duk_c_function func, duk_idx_t nargs);
duk_idx_t duk_push_c_lightfunc (duk_context ctx, duk_c_function func, duk_idx_t nargs, duk_idx_t length, duk_int_t magic);
duk_idx_t duk_push_thread_raw (duk_context ctx, duk_uint_t flags);
duk_idx_t duk_push_error_object_raw (duk_context ctx, duk_errcode_t err_code, const(char)* filename, duk_int_t line, const(char)* fmt, ...);
void* duk_push_buffer_raw (duk_context ctx, duk_size_t size, duk_small_uint_t flags);
void duk_push_buffer_object (duk_context ctx, duk_idx_t idx_buffer, duk_size_t byte_offset, duk_size_t byte_length, duk_uint_t flags);
duk_idx_t duk_push_heapptr (duk_context ctx, void* ptr);
void duk_pop (duk_context ctx);
void duk_pop_n (duk_context ctx, duk_idx_t count);
void duk_pop_2 (duk_context ctx);
void duk_pop_3 (duk_context ctx);
duk_int_t duk_get_type (duk_context ctx, duk_idx_t index);
duk_bool_t duk_check_type (duk_context ctx, duk_idx_t index, duk_int_t type);
duk_uint_t duk_get_type_mask (duk_context ctx, duk_idx_t index);
duk_bool_t duk_check_type_mask (duk_context ctx, duk_idx_t index, duk_uint_t mask);
duk_bool_t duk_is_undefined (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_null (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_null_or_undefined (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_boolean (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_number (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_nan (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_string (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_object (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_buffer (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_pointer (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_lightfunc (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_array (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_function (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_c_function (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_ecmascript_function (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_bound_function (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_thread (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_dynamic_buffer (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_fixed_buffer (duk_context ctx, duk_idx_t index);
duk_bool_t duk_is_external_buffer (duk_context ctx, duk_idx_t index);
duk_errcode_t duk_get_error_code (duk_context ctx, duk_idx_t index);
duk_bool_t duk_get_boolean (duk_context ctx, duk_idx_t index);
duk_double_t duk_get_number (duk_context ctx, duk_idx_t index);
duk_int_t duk_get_int (duk_context ctx, duk_idx_t index);
duk_uint_t duk_get_uint (duk_context ctx, duk_idx_t index);
const(char)* duk_get_string (duk_context ctx, duk_idx_t index);
const(char)* duk_get_lstring (duk_context ctx, duk_idx_t index, duk_size_t* out_len);
void* duk_get_buffer (duk_context ctx, duk_idx_t index, duk_size_t* out_size);
void* duk_get_buffer_data (duk_context ctx, duk_idx_t index, duk_size_t* out_size);
void* duk_get_pointer (duk_context ctx, duk_idx_t index);
duk_c_function duk_get_c_function (duk_context ctx, duk_idx_t index);
duk_context duk_get_context (duk_context ctx, duk_idx_t index);
void* duk_get_heapptr (duk_context ctx, duk_idx_t index);
duk_size_t duk_get_length (duk_context ctx, duk_idx_t index);
void duk_require_undefined (duk_context ctx, duk_idx_t index);
void duk_require_null (duk_context ctx, duk_idx_t index);
duk_bool_t duk_require_boolean (duk_context ctx, duk_idx_t index);
duk_double_t duk_require_number (duk_context ctx, duk_idx_t index);
duk_int_t duk_require_int (duk_context ctx, duk_idx_t index);
duk_uint_t duk_require_uint (duk_context ctx, duk_idx_t index);
const(char)* duk_require_string (duk_context ctx, duk_idx_t index);
const(char)* duk_require_lstring (duk_context ctx, duk_idx_t index, duk_size_t* out_len);
void* duk_require_buffer (duk_context ctx, duk_idx_t index, duk_size_t* out_size);
void* duk_require_buffer_data (duk_context ctx, duk_idx_t index, duk_size_t* out_size);
void* duk_require_pointer (duk_context ctx, duk_idx_t index);
duk_c_function duk_require_c_function (duk_context ctx, duk_idx_t index);
duk_context duk_require_context (duk_context ctx, duk_idx_t index);
void duk_require_function (duk_context ctx, duk_idx_t index);
void* duk_require_heapptr (duk_context ctx, duk_idx_t index);
void duk_to_undefined (duk_context ctx, duk_idx_t index);
void duk_to_null (duk_context ctx, duk_idx_t index);
duk_bool_t duk_to_boolean (duk_context ctx, duk_idx_t index);
duk_double_t duk_to_number (duk_context ctx, duk_idx_t index);
duk_int_t duk_to_int (duk_context ctx, duk_idx_t index);
duk_uint_t duk_to_uint (duk_context ctx, duk_idx_t index);
duk_int32_t duk_to_int32 (duk_context ctx, duk_idx_t index);
duk_uint32_t duk_to_uint32 (duk_context ctx, duk_idx_t index);
duk_uint16_t duk_to_uint16 (duk_context ctx, duk_idx_t index);
const(char)* duk_to_string (duk_context ctx, duk_idx_t index);
const(char)* duk_to_lstring (duk_context ctx, duk_idx_t index, duk_size_t* out_len);
void* duk_to_buffer_raw (duk_context ctx, duk_idx_t index, duk_size_t* out_size, duk_uint_t flags);
void* duk_to_pointer (duk_context ctx, duk_idx_t index);
void duk_to_object (duk_context ctx, duk_idx_t index);
void duk_to_defaultvalue (duk_context ctx, duk_idx_t index, duk_int_t hint);
void duk_to_primitive (duk_context ctx, duk_idx_t index, duk_int_t hint);
const(char)* duk_safe_to_lstring (duk_context ctx, duk_idx_t index, duk_size_t* out_len);
const(char)* duk_base64_encode (duk_context ctx, duk_idx_t index);
void duk_base64_decode (duk_context ctx, duk_idx_t index);
const(char)* duk_hex_encode (duk_context ctx, duk_idx_t index);
void duk_hex_decode (duk_context ctx, duk_idx_t index);
const(char)* duk_json_encode (duk_context ctx, duk_idx_t index);
void duk_json_decode (duk_context ctx, duk_idx_t index);
void* duk_resize_buffer (duk_context ctx, duk_idx_t index, duk_size_t new_size);
void* duk_steal_buffer (duk_context ctx, duk_idx_t index, duk_size_t* out_size);
void duk_config_buffer (duk_context ctx, duk_idx_t index, void* ptr, duk_size_t len);
duk_bool_t duk_get_prop (duk_context ctx, duk_idx_t obj_index);
duk_bool_t duk_get_prop_string (duk_context ctx, duk_idx_t obj_index, const(char)* key);
duk_bool_t duk_get_prop_index (duk_context ctx, duk_idx_t obj_index, duk_uarridx_t arr_index);
duk_bool_t duk_put_prop (duk_context ctx, duk_idx_t obj_index);
duk_bool_t duk_put_prop_string (duk_context ctx, duk_idx_t obj_index, const(char)* key);
duk_bool_t duk_put_prop_index (duk_context ctx, duk_idx_t obj_index, duk_uarridx_t arr_index);
duk_bool_t duk_del_prop (duk_context ctx, duk_idx_t obj_index);
duk_bool_t duk_del_prop_string (duk_context ctx, duk_idx_t obj_index, const(char)* key);
duk_bool_t duk_del_prop_index (duk_context ctx, duk_idx_t obj_index, duk_uarridx_t arr_index);
duk_bool_t duk_has_prop (duk_context ctx, duk_idx_t obj_index);
duk_bool_t duk_has_prop_string (duk_context ctx, duk_idx_t obj_index, const(char)* key);
duk_bool_t duk_has_prop_index (duk_context ctx, duk_idx_t obj_index, duk_uarridx_t arr_index);
void duk_def_prop (duk_context ctx, duk_idx_t obj_index, duk_uint_t flags);
duk_bool_t duk_get_global_string (duk_context ctx, const(char)* key);
duk_bool_t duk_put_global_string (duk_context ctx, const(char)* key);
void duk_get_prototype (duk_context ctx, duk_idx_t index);
void duk_set_prototype (duk_context ctx, duk_idx_t index);
void duk_get_finalizer (duk_context ctx, duk_idx_t index);
void duk_set_finalizer (duk_context ctx, duk_idx_t index);
void duk_set_global_object (duk_context ctx);
duk_int_t duk_get_magic (duk_context ctx, duk_idx_t index);
void duk_set_magic (duk_context ctx, duk_idx_t index, duk_int_t magic);
duk_int_t duk_get_current_magic (duk_context ctx);
void duk_put_function_list (duk_context ctx, duk_idx_t obj_index, const(duk_function_list_entry)* funcs);
void duk_put_number_list (duk_context ctx, duk_idx_t obj_index, const(duk_number_list_entry)* numbers);
void duk_get_var (duk_context ctx);
void duk_put_var (duk_context ctx);
duk_bool_t duk_del_var (duk_context ctx);
duk_bool_t duk_has_var (duk_context ctx);
void duk_compact (duk_context ctx, duk_idx_t obj_index);
void duk_enum (duk_context ctx, duk_idx_t obj_index, duk_uint_t enum_flags);
duk_bool_t duk_next (duk_context ctx, duk_idx_t enum_index, duk_bool_t get_value);
void duk_concat (duk_context ctx, duk_idx_t count);
void duk_join (duk_context ctx, duk_idx_t count);
void duk_decode_string (duk_context ctx, duk_idx_t index, duk_decode_char_function callback, void* udata);
void duk_map_string (duk_context ctx, duk_idx_t index, duk_map_char_function callback, void* udata);
void duk_substring (duk_context ctx, duk_idx_t index, duk_size_t start_char_offset, duk_size_t end_char_offset);
void duk_trim (duk_context ctx, duk_idx_t index);
duk_codepoint_t duk_char_code_at (duk_context ctx, duk_idx_t index, duk_size_t char_offset);
duk_bool_t duk_equals (duk_context ctx, duk_idx_t index1, duk_idx_t index2);
duk_bool_t duk_strict_equals (duk_context ctx, duk_idx_t index1, duk_idx_t index2);
duk_bool_t duk_instanceof (duk_context ctx, duk_idx_t index1, duk_idx_t index2);
void duk_call (duk_context ctx, duk_idx_t nargs);
void duk_call_method (duk_context ctx, duk_idx_t nargs);
void duk_call_prop (duk_context ctx, duk_idx_t obj_index, duk_idx_t nargs);
duk_int_t duk_pcall (duk_context ctx, duk_idx_t nargs);
duk_int_t duk_pcall_method (duk_context ctx, duk_idx_t nargs);
duk_int_t duk_pcall_prop (duk_context ctx, duk_idx_t obj_index, duk_idx_t nargs);
void duk_new (duk_context ctx, duk_idx_t nargs);
duk_int_t duk_pnew (duk_context ctx, duk_idx_t nargs);
duk_int_t duk_safe_call (duk_context ctx, duk_safe_call_function func, duk_idx_t nargs, duk_idx_t nrets);
duk_int_t duk_eval_raw (duk_context ctx, const(char)* src_buffer, duk_size_t src_length, duk_uint_t flags);
duk_int_t duk_compile_raw (duk_context ctx, const(char)* src_buffer, duk_size_t src_length, duk_uint_t flags);
void duk_dump_function (duk_context ctx);
void duk_load_function (duk_context ctx);
void duk_log (duk_context ctx, duk_int_t level, const(char)* fmt, ...);
void duk_push_context_dump (duk_context ctx);
void duk_debugger_attach_custom (duk_context ctx, duk_debug_read_function read_cb, duk_debug_write_function write_cb, duk_debug_peek_function peek_cb, duk_debug_read_flush_function read_flush_cb, duk_debug_write_flush_function write_flush_cb, duk_debug_request_function request_cb, duk_debug_detached_function detached_cb, void* udata);
void duk_debugger_detach (duk_context ctx);
void duk_debugger_cooperate (duk_context ctx);
duk_bool_t duk_debugger_notify (duk_context ctx, duk_idx_t nvalues);
void duk_debugger_pause (duk_context ctx);
}

duk_context duk_create_heap_default() {
  return duk_create_heap(null, null, null, null, null);
}

void duk_eval_string (duk_context ctx, const(char)* src) {
  duk_eval_raw((ctx), (src), 0, 1 /*args*/ | DUK_COMPILE_EVAL | DUK_COMPILE_NOSOURCE | DUK_COMPILE_STRLEN | DUK_COMPILE_NOFILENAME);
}

import std.traits;

void push(T)(duk_context ctx, T value) if (isBoolean!T) {
  duk_push_boolean(ctx, value);
}

void push(T)(duk_context ctx, T value) if (isNumeric!T) {
  duk_push_number(ctx, value);
}

void push(T)(duk_context ctx, T value) if (isSomeChar!T) {
  duk_push_lstring(ctx, &value, 1);
}

void push(T)(duk_context ctx, T value) if (is(T == string)) {
  duk_push_lstring(ctx, value.ptr, value.length);
}
