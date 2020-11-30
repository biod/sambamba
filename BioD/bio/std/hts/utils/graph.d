/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.std.hts.utils.graph;

import std.exception;
import std.array;
import std.algorithm;
import std.range;

/// This class is used only in utils/samheadermerger.d for merging
/// reference sequence dictionaries, thus not much functionality 
/// is presented here.
class DirectedGraph {

    void addEdge(string a, string b) {
        auto indA = addNode(a);
        auto indB = addNode(b);
        _edges[indA] ~= indB;
    }

    /// Returns: unique integer identifier of node
    size_t addNode(string a) {
        if (a !in _indices) {
            _nodes ~= a;
            _edges.length = _nodes.length;
            return _indices[a] = _nodes.length - 1;
        } else {
            return _indices[a];
        }
    }

    string[] topologicalSort() {
        auto predecessor_count = new size_t[_nodes.length];
        foreach (node_neighbours; _edges) {
            foreach (neighbour; node_neighbours) 
                predecessor_count[neighbour] += 1;
        }

        size_t[] queue;
        queue.reserve(_nodes.length);
        for (size_t i = 0; i < _nodes.length; i++)
            if (predecessor_count[i] == 0)
                queue ~= i;
      
        string[] result;
        result.reserve(_nodes.length);

        while (!queue.empty) {
            auto front = queue[0];
            result ~= _nodes[front];
            queue = queue[1 .. $];
            foreach (successor; _edges[front]) {
                predecessor_count[successor] -= 1;
                if (predecessor_count[successor] == 0) {
                    queue ~= successor;
                }
            }
        }
        
        if (result.length < _nodes.length) {
            throw new Exception("graph contains cycles");
        }

        return result;
    }

private:
    string[] _nodes;
    size_t[string] _indices;
    size_t[][] _edges;
}
