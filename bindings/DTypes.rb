require 'ffi'

class DArray < FFI::Struct 
	layout :length, :size_t,
		   :ptr, :pointer
end

DString = DArray

class DHash < FFI::Struct
    layout :keys_length, :size_t,
           :keys_ptr,    :pointer,
           :vals_length, :size_t,
           :vals_ptr,    :pointer

    def initialize(ptr)
        super(ptr)
        @ptr = ptr
        ObjectSpace.define_finalizer @ptr, DHash.finalize(ptr)
    end
   
    def to_ruby_value
        hash = {}
        i = 0
        len = self[:keys_length]
        kptr = self[:keys_ptr]
        vptr = self[:vals_ptr]
        while i < len do
            dkey = DString.new(kptr + i * DString.size)
            key = dkey[:ptr].read_string(dkey[:length])
                                                        # owned
            dval = TagValue.new(vptr + i * TagValue.size, true)
            hash[key] = dval.to_ruby_value
            i += 1
        end
        hash
    end

    private

    def self.finalize(ptr)
        proc { LibBAM.dhash_destroy ptr }
    end

end
