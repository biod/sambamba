module msgpack.exception;

@trusted:

/**
 * $(D MessagePackException) is a root Exception for MessagePack related operation.
 */
class MessagePackException : Exception
{
    pure this(string message)
    {
        super(message);
    }
}


/**
 * $(D UnpackException) is thrown on deserialization failure
 */
class UnpackException : MessagePackException
{
    this(string message)
    {
        super(message);
    }
}

