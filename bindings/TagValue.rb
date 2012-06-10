require 'ffi'
require './bindings/DTypes.rb'

class TagValueU < FFI::Union
    layout :A, :char,
           :c, :int8,
           :C, :uint8,
           :s, :int16,
           :S, :uint16,
           :i, :int32,
           :I, :uint32,
           :f, :float,
           :array, DArray # for Z, H, B:?
end

class TagValue < FFI::Struct
    layout :u,          TagValueU,
           :bam_typeid, :char,
           :tag,        :uint8


    def initialize(ptr, owned=false)
        super(ptr)

        unless owned
            @ptr = ptr
            ObjectSpace.define_finalizer(self, TagValue.finalize(ptr))
        end
    end

    def to_ruby_value
        if nothing? then
            nil
        elsif string? then
            arr = self[:u][:array]
            arr[:ptr].read_string(arr[:length])
        elsif array? then
            arr = self[:u][:array]
            arr[:ptr].send @@bam_typeid_to_method[self[:bam_typeid].chr], 
                           arr[:length]
        else
            self[:u][self[:bam_typeid].chr.to_sym]
        end
    end

    def string?
        self[:tag] & 0b11 == 0b11
    end

    def array?
        self[:tag] & 0b11 == 0b01
    end
    
    def nothing?
        self[:tag] == 0b10
    end

    private
    @@bam_typeid_to_method = { 'c' => :read_array_of_int8,
                               'C' => :read_array_of_uint8,
                               's' => :read_array_of_int16,
                               'S' => :read_array_of_uint16,
                               'i' => :read_array_of_int32,
                               'I' => :read_array_of_uint32,
                               'f' => :read_array_of_float
                             }

    def self.finalize(ptr)
        proc { LibBAM.tag_value_destroy ptr }
    end
end
