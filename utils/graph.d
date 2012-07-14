module utils.graph;

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
        auto queue = array(
            filter!((size_t i) { return predecessor_count[i] == 0; })
                   (iota(_nodes.length))
        );
      
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
