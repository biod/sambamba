// Written in the D programming language.

/**
 * User-defined class sample
 */

import std.stdio, std.math;

import msgpack;


enum Attr : ubyte
{
    A, B
}


struct User
{
    string name;
    uint   age;
    Attr   attr;
}


void main()
{
    User user = User("Foo", 20, Attr.B), other;

    unpack(pack(user), other);

    writeln("name: ", other.name, "(", other.age, ", ", other.attr, ")");

    // Complex data-structure

    auto     btree  = new BTree!(int, string);
    int[]    keys   = [3,     6,     8,     10,     1,      5,      20];
    string[] values = ["Foo", "Baz", "Bar", "Hoge", "Fuga", "Piyo", "ham"];

    foreach (i, key; keys)
        btree.insert(key, values[i]);
    btree.print();

    auto result = new BTree!(int, string);

    unpack(pack(btree), result);

    result.print();
}


/**
 * Simple B-Tree.
 */
class BTree(Key, Data)
{
  private:
    immutable uint MAX_CHILD;   // sub-tree size(m-tree)
    immutable uint HALF_CHILD;  // min number((m + 1) / 2)

    enum NodeType
    {
        Internal,
        Leaf
    }

    static class Node
    {
        NodeType type;

        union
        {
            struct  // internal
            {
                uint   childs; // child number
                Key[]  low;    // min element of each sub-tree
                Node[] child;  // sub-tree size(m)
            }

            struct  // leaf
            {
                Key  key;
                Data data;
            }
        }

        this(in uint num) { low.length = child.length = num; }

        void toMsgpack(Packer)(ref Packer packer) const
        {
            if (type == NodeType.Internal)
                packer.packArray(childs, low.length, low, child);
            else
                packer.packArray(key, data);
        }

        void fromMsgpack(ref Unpacker unpacker)
        {
            if (unpacker.beginArray() == 4) {
                uint max;

                unpacker.unpack(childs, max, low);
                child.length = unpacker.beginArray();
                foreach (i; 0..child.length)
                    unpacker.unpack(child[i], max);
                low.length = child.length = max;
            } else {
                type = NodeType.Leaf;
                unpacker.unpack(key, data);
            }
        }
    }

    Node root;


  public:
    /**
     * Params:
     *  num = node size.
     */
    this(in uint num = 5)
    in
    {
        assert(num >= 2, "The prerequisite of m-tree(size must be larger than 2)");
    }
    body
    {
        MAX_CHILD  = num;
        HALF_CHILD = cast(uint)ceil(cast(real)MAX_CHILD / 2);
    }

    /**
     * Params:
     *  key  = key that represents a data to insert.
     *  data = data to insert.
     *
     * Returns:
     *  node that inserted data.
     */
    Node insert(in Key key, in Data data)
    {
        /*
         * Params:
         *  n      = node to insert element
         *  nNode  = new node(null if not create).
         *  lowest = min element in $(D_PARAM nNode)
         *
         * Returns:
         *  node that element was inserted.
         */
        Node _insert(ref Node n, out Node nNode, out Key lowest)
        {
            Node node = n; // for updating data

            if (node.type == NodeType.Leaf) {
                if (node.key == key) {
                    return null;
                } else {
                    Node elem = allocNode();
                    elem.type = NodeType.Leaf;
                    elem.key  = key;
                    elem.data = data;

                    if (elem.key < node.key) {
                        n      = elem;
                        nNode  = node;
                        lowest = node.key;
                    } else {
                        nNode  = elem;
                        lowest = elem.key;
                    }

                    return elem;
                }
            } else {
                int  i, j, pos;  // pos = position to insert.
                Key  xLowest;    // lowest for recursion.
                Node xNode, ret; // nNode for recursion.

                pos = locateSubTree(node, key);
                ret = _insert(node.child[pos], xNode, xLowest);

                // Doesn't create and patition.
                if (xNode is null)
                    return ret;

                if (node.childs < MAX_CHILD) {
                    for (i = node.childs - 1; i > pos; i--) {
                        node.child[i+1] = node.child[i];
                        node.low[i+1]   = node.low[i];
                    }
                    node.child[pos+1] = xNode;
                    node.low[pos+1]   = xLowest;
                    node.childs++;
                    return ret;
                } else {
                    Node elem = allocNode();
                    elem.type = NodeType.Internal;

                    // insert to node side or elem side?
                    if (pos < HALF_CHILD - 1) {
                        for (i = HALF_CHILD - 1, j = 0; i < MAX_CHILD; i++, j++) {
                            elem.child[j] = node.child[i];
                            elem.low[j]   = node.low[i];
                        }

                        for (i = HALF_CHILD - 2; i > pos; i--) {
                            node.child[i+1] = node.child[i];
                            node.low[i+1]   = node.low[i];
                        }
                        node.child[pos+1] = xNode;
                        node.low[pos+1]   = xLowest;
                    } else {
                        for (i = MAX_CHILD - 1, j = MAX_CHILD - HALF_CHILD; i >= HALF_CHILD; i--) {
                            if (i == pos) {
                                elem.child[j] = xNode;
                                elem.low[j--] = xLowest;
                            }
                            elem.child[j] = node.child[i];
                            elem.low[j--] = node.low[i];
                        }

                        if (pos < HALF_CHILD) {
                            elem.child[0] = xNode;
                            elem.low[0]   = xLowest;
                        }
                    }

                    node.childs = HALF_CHILD;
                    elem.childs = MAX_CHILD+1 - HALF_CHILD;

                    nNode  = elem;
                    lowest = elem.low[0];

                    return ret;
                }
            }
        }

        if (root is null) {
            root      = allocNode();
            root.type = NodeType.Leaf;
            root.key  = key;
            root.data = data;
            return root;
        } else {
            Key  lowest;
            Node ret, newNode;

            ret = _insert(root, newNode, lowest);

            // new node and growl height if patitioned.
            if (newNode !is null) {
                Node elem     = allocNode();
                elem.type     = NodeType.Internal;
                elem.childs   = 2;
                elem.child[0] = root;
                elem.child[1] = newNode;
                elem.low[1]   = lowest;
                root          = elem;
            }

            return ret;
        }
    }

    /**
     * Params:
     *  key = key to delete.
     *
     * Returns:
     *  true if deleted.
     */
    bool remove(in Key key)
    {
        enum State
        {
            Nothing,
            Removed,
            Reconfig
        }

        /*
         * Params:
         *  n      = node to remove.
         *  result = node change because of an remove.
         *
         * Returns:
         *  true if removed.
         */
        bool _remove(ref Node n, out State result)
        {
            if (n.type == NodeType.Leaf) {
                if (n.key == key) {
                    result = State.Removed;
                    delete n;
                    return true;
                } else {
                    return false;
                }
            } else {
                int   pos, sub;    // sub-tree position to remove, for restructure.
                bool  ret, merged; // delete?, merge sub-tree?
                State state;       // sub-tree state.

                pos = locateSubTree(n, key);
                ret = _remove(n.child[pos], state);

                if (state == State.Nothing)
                    return ret;

                if (state == State.Reconfig) {
                    sub    = pos == 0 ? 0 : pos - 1;
                    merged = revisionNodes(n, sub);

                    if (merged)
                        pos = sub+1;
                }

                if (state == State.Removed || merged) {
                    // sub-tree compaction
                    for (int i = pos; i < n.childs - 1; i++) {
                        n.child[i] = n.child[i+1];
                        n.low[i]   = n.low[i+1];
                    }

                    if (--n.childs < HALF_CHILD)
                        result = State.Reconfig;
                }

                return ret;
            }
        }

        if (root is null) {
            return false;
        } else {
            State result;
            bool  ret = _remove(root, result);

            if (result == State.Removed) {
                root = null;
            } else if (result == State.Reconfig && root.childs == 1) {
                Node n = root;
                root   = root.child[0];
                delete n;
            }

            return ret;
        }
    }

    /**
     * Params:
     *  key = key to search.
     *
     * Returns:
     *  finded node.
     */
    Node search(in Key key)
    {
        if (root is null) {
            return null;
        } else {
            int  i;
            Node n = root;

            // searches internal node until find leaf.
            while (n.type == NodeType.Internal) {
                i = locateSubTree(n, key);
                n = n.child[i];
            }

            if (key == n.key)
                return n;
            else
                return null;
        }
    }

    void print()
    {
        void _print(ref Node n)
        {
            if (n.type == NodeType.Leaf) {
                writefln("[%x] Leaf : %s Data : %s", &n, n.key, n.data);
            } else {
                writef("[%x] Childs %d [%x], ", &n, n.childs, &n.child[0]);

                foreach (i; 0..MAX_CHILD)
                    writef("%d[%x] ", n.low[i], &n.child[i]);
                writeln();

                foreach (i; 0..n.childs)
                    _print(n.child[i]);
            }
        }

        if (root is null)
            writefln("Element is nothing");
        else
            _print(root);
    }

    // for MessagePack

    void toMsgpack(Packer)(ref Packer packer) const
    {
        packer.pack(root);
    }

    void fromMsgpack(ref Unpacker unpacker)
    {
        unpacker.unpack(root, MAX_CHILD);
    }


  private:
    /*
     * Returns:
     *  new node.
     */
    Node allocNode()
    {
        return new Node(MAX_CHILD);
    }

    /*
     * searches $(D_PARAM key) element in sub-tree of $(D_PARAM node).
     *
     * Params:
     *  node = node to search.
     *  key  = key to search in $(D_PARAM node).
     *
     * Returns:
     *  finded position.
     */
    int locateSubTree(ref Node node, in Key key) const
    {
        for (int i = node.childs - 1; i > 0; i--)
            if (key >= node.low[i])
                return i;
        return 0;
    }

    /*
     * Params:
     *  n = revision node.
     *  x = position to sub-tree of revision node.
     *
     * Returns:
     *  true if merged.
     */
    bool revisionNodes(ref Node n, in uint x)
    {
        int  i;
        Node a,  b;   // each sub-tree.
        uint an, bn;  // child number of each sub-tree.

        a  = n.child[x];
        b  = n.child[x+1];
        an = a.childs;
        bn = b.childs;
        b.low[0] = n.low[x+1];

        if (an + bn <= MAX_CHILD) {  // merge
            for (i = 0; i < bn; i++) {
                a.child[an+i] = b.child[i];
                a.low[an+i]   = b.low[i];
            }

            a.childs += bn;
            delete b;

            return true;
        } else {  // partition
            uint pivot = (an + bn) / 2; // pivot to patition.
            uint move;                  // element number to copy.

            if (an > pivot) {
                move = an - pivot;

                for (i = bn - 1; i >= 0; i--) {
                    b.child[move+i] = b.child[i];
                    b.low[move+i]   = b.low[i];
                }
                // copy element a to b.
                for (i = 0; i < move; i++) {
                    b.child[i] = a.child[pivot+i];
                    b.low[i]   = a.low[pivot+i];
                }
            } else {
                move = pivot - an;

                // copy element b to a.
                for (i = 0; i < move; i++) {
                    a.child[an+i] = b.child[i];
                    a.low[an+i]   = b.low[i];
                }
                for (i = 0; i < bn - move; i++) {
                    b.child[i] = b.child[move+i];
                    b.low[i]   = b.low[move+i];
                }
            }

            a.childs   = pivot;
            b.childs   = an + bn - pivot;
            n.low[x+1] = b.low[0];

            return false;
        }
    }
}
