// Directed graph with Tarjan's SCC algorithm.
// Simplified port of the C++ DiGraph — only implements features used by seq.rs.

use std::collections::HashMap;

/// Node index type
pub type NodeId = usize;

/// Directed graph
pub struct DiGraph {
    pub nodes: Vec<Node>,
    pub arcs: Vec<Arc>,
}

/// Node in the directed graph
pub struct Node {
    pub id: NodeId,
    pub out_arcs: Vec<usize>,  // indices into DiGraph::arcs
    pub in_arcs: Vec<usize>,
    pub scc: Option<NodeId>,   // SCC representative
    pub order_dfs: usize,      // 0 = not visited
    pub reachable: bool,
    // For SCC computation
    in_stack: bool,
    low_link: usize,
}

/// Arc (directed edge) in the graph
pub struct Arc {
    pub from: NodeId,
    pub to: NodeId,
}

impl DiGraph {
    pub fn new() -> Self {
        DiGraph {
            nodes: Vec::new(),
            arcs: Vec::new(),
        }
    }

    pub fn add_node(&mut self) -> NodeId {
        let id = self.nodes.len();
        self.nodes.push(Node {
            id,
            out_arcs: Vec::new(),
            in_arcs: Vec::new(),
            scc: None,
            order_dfs: 0,
            reachable: false,
            in_stack: false,
            low_link: 0,
        });
        id
    }

    pub fn add_arc(&mut self, from: NodeId, to: NodeId) -> usize {
        let arc_id = self.arcs.len();
        self.arcs.push(Arc { from, to });
        self.nodes[from].out_arcs.push(arc_id);
        self.nodes[to].in_arcs.push(arc_id);
        arc_id
    }

    /// Compute strongly connected components using Tarjan's algorithm.
    /// Sets node.scc to the representative node of each SCC.
    pub fn scc(&mut self) {
        let mut visited_num: usize = 0;
        let mut stack: Vec<NodeId> = Vec::new();

        // Need to collect node ids first to avoid borrow issues
        let node_ids: Vec<NodeId> = (0..self.nodes.len()).collect();

        for &node_id in &node_ids {
            if self.nodes[node_id].order_dfs == 0 {
                self.scc_visit(node_id, &mut visited_num, &mut stack);
            }
        }
    }

    fn scc_visit(
        &mut self,
        node_id: NodeId,
        visited_num: &mut usize,
        stack: &mut Vec<NodeId>,
    ) {
        *visited_num += 1;
        let order = *visited_num;
        self.nodes[node_id].order_dfs = order;
        self.nodes[node_id].low_link = order;
        self.nodes[node_id].in_stack = true;
        stack.push(node_id);

        // Get outgoing neighbors
        let out_arcs: Vec<usize> = self.nodes[node_id].out_arcs.clone();

        for arc_id in out_arcs {
            let target = self.arcs[arc_id].to;
            if self.nodes[target].order_dfs == 0 {
                self.scc_visit(target, visited_num, stack);
                let target_low = self.nodes[target].low_link;
                if target_low < self.nodes[node_id].low_link {
                    self.nodes[node_id].low_link = target_low;
                }
            } else if self.nodes[target].in_stack {
                let target_order = self.nodes[target].order_dfs;
                if target_order < self.nodes[node_id].low_link {
                    self.nodes[node_id].low_link = target_order;
                }
            }
        }

        // Root of SCC
        if self.nodes[node_id].low_link == self.nodes[node_id].order_dfs {
            loop {
                let w = stack.pop().unwrap();
                self.nodes[w].in_stack = false;
                self.nodes[w].scc = Some(node_id);
                if w == node_id {
                    break;
                }
            }
        }
    }

    /// Get connected components using union-find
    pub fn connected_components(&mut self) -> HashMap<NodeId, Vec<NodeId>> {
        let mut parent: Vec<NodeId> = (0..self.nodes.len()).collect();
        let mut rank: Vec<usize> = vec![0; self.nodes.len()];

        fn find(parent: &mut [NodeId], x: NodeId) -> NodeId {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }

        fn union(parent: &mut [NodeId], rank: &mut [usize], x: NodeId, y: NodeId) {
            let rx = find(parent, x);
            let ry = find(parent, y);
            if rx == ry {
                return;
            }
            if rank[rx] < rank[ry] {
                parent[rx] = ry;
            } else if rank[rx] > rank[ry] {
                parent[ry] = rx;
            } else {
                parent[ry] = rx;
                rank[rx] += 1;
            }
        }

        for arc in &self.arcs {
            union(&mut parent, &mut rank, arc.from, arc.to);
        }

        let mut components: HashMap<NodeId, Vec<NodeId>> = HashMap::new();
        for i in 0..self.nodes.len() {
            let root = find(&mut parent, i);
            components.entry(root).or_default().push(i);
        }

        components
    }

    /// Clear reachability flags
    pub fn clear_reachable(&mut self) {
        for node in &mut self.nodes {
            node.reachable = false;
        }
    }

    /// Mark nodes reachable from the given node (following outgoing arcs)
    pub fn set_reachable(&mut self, start: NodeId) {
        self.nodes[start].reachable = true;
        let out_arcs: Vec<usize> = self.nodes[start].out_arcs.clone();
        for arc_id in out_arcs {
            let target = self.arcs[arc_id].to;
            if !self.nodes[target].reachable {
                self.set_reachable(target);
            }
        }
    }

    /// Get root nodes (nodes with no incoming arcs)
    pub fn roots(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|n| n.in_arcs.is_empty())
            .map(|n| n.id)
            .collect()
    }

    /// Get leaf nodes (nodes with no outgoing arcs)
    pub fn leaves(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|n| n.out_arcs.is_empty())
            .map(|n| n.id)
            .collect()
    }
}

impl Default for DiGraph {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scc_simple() {
        let mut g = DiGraph::new();
        let a = g.add_node();
        let b = g.add_node();
        let c = g.add_node();

        // a -> b -> c -> a (one SCC)
        g.add_arc(a, b);
        g.add_arc(b, c);
        g.add_arc(c, a);

        g.scc();

        // All three should have the same SCC representative
        assert_eq!(g.nodes[a].scc, g.nodes[b].scc);
        assert_eq!(g.nodes[b].scc, g.nodes[c].scc);
    }

    #[test]
    fn test_scc_two_components() {
        let mut g = DiGraph::new();
        let a = g.add_node();
        let b = g.add_node();
        let c = g.add_node();
        let d = g.add_node();

        // a -> b -> a (one SCC)
        g.add_arc(a, b);
        g.add_arc(b, a);
        // c -> d (two separate SCCs)
        g.add_arc(c, d);

        g.scc();

        assert_eq!(g.nodes[a].scc, g.nodes[b].scc);
        assert_ne!(g.nodes[c].scc, g.nodes[d].scc);
    }

    #[test]
    fn test_reachable() {
        let mut g = DiGraph::new();
        let a = g.add_node();
        let b = g.add_node();
        let c = g.add_node();
        let d = g.add_node();

        g.add_arc(a, b);
        g.add_arc(b, c);
        // d is isolated

        g.set_reachable(a);

        assert!(g.nodes[a].reachable);
        assert!(g.nodes[b].reachable);
        assert!(g.nodes[c].reachable);
        assert!(!g.nodes[d].reachable);
    }

    #[test]
    fn test_roots_leaves() {
        let mut g = DiGraph::new();
        let a = g.add_node();
        let b = g.add_node();
        let c = g.add_node();

        g.add_arc(a, b);
        g.add_arc(b, c);

        assert_eq!(g.roots(), vec![a]);
        assert_eq!(g.leaves(), vec![c]);
    }
}
