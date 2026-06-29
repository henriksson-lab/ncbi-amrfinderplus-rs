// Directed graph with Tarjan's SCC algorithm.
// Port of amr/graph.{hpp,cpp} DiGraph's core operations.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::fmt::Write;

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
    pub out_arcs: Vec<usize>, // indices into DiGraph::arcs
    pub in_arcs: Vec<usize>,
    pub parent_dc: NodeId,
    pub scc: Option<NodeId>, // SCC representative
    pub order_dfs: usize,    // 0 = not visited
    pub reachable: bool,
    // For SCC computation
    in_stack: bool,
}

/// Arc (directed edge) in the graph
pub struct Arc {
    pub from: NodeId,
    pub to: NodeId,
    pub active: bool,
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
            parent_dc: id,
            scc: None,
            order_dfs: 0,
            reachable: false,
            in_stack: false,
        });
        id
    }

    pub fn add_arc(&mut self, from: NodeId, to: NodeId) -> usize {
        let arc_id = self.arcs.len();
        self.arcs.push(Arc {
            from,
            to,
            active: true,
        });
        self.nodes[from].out_arcs.push(arc_id);
        self.nodes[to].in_arcs.push(arc_id);
        arc_id
    }

    pub fn empty(&self) -> bool {
        self.nodes.is_empty()
    }

    pub fn clear(&mut self) {
        self.delete_nodes();
    }

    pub fn qc(&self) -> Result<(), String> {
        for (node_id, node) in self.nodes.iter().enumerate() {
            if node.id != node_id {
                return Err(format!("node {node_id}: stored id {}", node.id));
            }
            if node.parent_dc >= self.nodes.len() {
                return Err(format!("node {node_id}: bad disjoint-cluster parent"));
            }
            if let Some(scc) = node.scc {
                if scc >= self.nodes.len() {
                    return Err(format!("node {node_id}: bad SCC"));
                }
            }
            for &arc_id in &node.out_arcs {
                let Some(arc) = self.arcs.get(arc_id) else {
                    return Err(format!("node {node_id}: bad outgoing arc"));
                };
                if arc.active && arc.from != node_id {
                    return Err(format!("node {node_id}: outgoing arc endpoint mismatch"));
                }
            }
            for &arc_id in &node.in_arcs {
                let Some(arc) = self.arcs.get(arc_id) else {
                    return Err(format!("node {node_id}: bad incoming arc"));
                };
                if arc.active && arc.to != node_id {
                    return Err(format!("node {node_id}: incoming arc endpoint mismatch"));
                }
            }
        }
        for (arc_id, arc) in self.arcs.iter().enumerate() {
            if !arc.active {
                continue;
            }
            if arc.from >= self.nodes.len() || arc.to >= self.nodes.len() {
                return Err(format!("arc {arc_id}: bad endpoint"));
            }
            if !self.nodes[arc.from].out_arcs.contains(&arc_id) {
                return Err(format!("arc {arc_id}: missing from outgoing list"));
            }
            if !self.nodes[arc.to].in_arcs.contains(&arc_id) {
                return Err(format!("arc {arc_id}: missing from incoming list"));
            }
        }
        Ok(())
    }

    pub fn save_text(&self) -> String {
        let mut s = String::new();
        for node in &self.nodes {
            let _ = writeln!(s, "{}", self.node_save_text(node.id));
        }
        s
    }

    fn node_save_text(&self, node: NodeId) -> String {
        let mut s = String::new();
        let n = &self.nodes[node];
        let _ = write!(s, "{node:p}", node = &n.id);
        if n.order_dfs != 0 {
            let _ = write!(s, "  DFS_order = {}", n.order_dfs);
        }
        if let Some(scc) = n.scc {
            let _ = write!(s, "  SCC: {scc:p}", scc = &self.nodes[scc].id);
        }
        let _ = writeln!(s);
        let _ = writeln!(s, "  In:");
        for &arc_id in &n.in_arcs {
            if self.arcs[arc_id].active {
                let _ = writeln!(s, "    {:p}: ", &self.nodes[self.arcs[arc_id].from].id);
            }
        }
        let _ = writeln!(s, "  Out:");
        for &arc_id in &n.out_arcs {
            if self.arcs[arc_id].active {
                let _ = writeln!(s, "    {:p}: ", &self.nodes[self.arcs[arc_id].to].id);
            }
        }
        s
    }

    pub fn delete_nodes(&mut self) {
        self.nodes.clear();
        self.arcs.clear();
    }

    pub fn copy(&self) -> Self {
        let mut graph = DiGraph::new();
        let mut old2new = HashMap::with_capacity(self.nodes.len());
        graph.init(self, &mut old2new);
        graph
    }

    pub fn init(&mut self, other: &DiGraph, old2new: &mut HashMap<NodeId, NodeId>) {
        debug_assert!(old2new.is_empty());
        self.delete_nodes();

        for node in &other.nodes {
            let new_node = self.add_node();
            self.nodes[new_node].parent_dc = node.parent_dc;
            self.nodes[new_node].scc = node.scc;
            self.nodes[new_node].order_dfs = node.order_dfs;
            self.nodes[new_node].reachable = node.reachable;
            self.nodes[new_node].in_stack = node.in_stack;
            old2new.insert(node.id, new_node);
        }

        for node in &mut self.nodes {
            node.parent_dc = old2new[&node.parent_dc];
            if let Some(scc) = node.scc {
                node.scc = Some(old2new[&scc]);
            }
        }

        for arc in other.arcs.iter().filter(|arc| arc.active) {
            self.add_arc(old2new[&arc.from], old2new[&arc.to]);
        }
    }

    pub fn connected_components(&mut self) {
        self.disjoint_cluster_init_all();

        let active_arcs: Vec<usize> = self
            .arcs
            .iter()
            .enumerate()
            .filter_map(|(i, arc)| arc.active.then_some(i))
            .collect();
        for arc_id in active_arcs {
            let mut root_from = self.arcs[arc_id].from;
            while self.nodes[root_from].parent_dc != root_from {
                root_from = self.nodes[root_from].parent_dc;
            }

            let mut root_to = self.arcs[arc_id].to;
            while self.nodes[root_to].parent_dc != root_to {
                root_to = self.nodes[root_to].parent_dc;
            }

            if root_from != root_to {
                self.disjoint_cluster_merge_roots(root_from, root_to);
            }
        }
    }

    pub fn disjoint_cluster_init_all(&mut self) {
        for i in 0..self.nodes.len() {
            self.nodes[i].parent_dc = i;
        }
    }

    pub fn disjoint_cluster_root(&mut self, node: NodeId) -> NodeId {
        let parent = self.nodes[node].parent_dc;
        if parent == node {
            node
        } else {
            let root = self.disjoint_cluster_root(parent);
            self.nodes[node].parent_dc = root;
            root
        }
    }

    pub fn disjoint_cluster_root_const(&self, mut node: NodeId) -> NodeId {
        while self.nodes[node].parent_dc != node {
            node = self.nodes[node].parent_dc;
        }
        node
    }

    pub fn disjoint_cluster_merge(&mut self, a: NodeId, b: NodeId) {
        let root_a = self.disjoint_cluster_root(a);
        let root_b = self.disjoint_cluster_root(b);
        if root_a != root_b {
            self.disjoint_cluster_merge_roots(root_a, root_b);
        }
    }

    fn disjoint_cluster_merge_roots(&mut self, root_a: NodeId, root_b: NodeId) {
        self.nodes[root_b].parent_dc = root_a;
    }

    /// Compute strongly connected components using Tarjan's algorithm.
    /// Sets node.scc to the representative node of each SCC.
    pub fn scc(&mut self) {
        let mut visited_num: usize = 0;
        let mut stack: Vec<NodeId> = Vec::with_capacity(self.nodes.len());

        for node in &mut self.nodes {
            node.scc = None;
            node.order_dfs = 0;
            node.in_stack = false;
        }

        for node_id in 0..self.nodes.len() {
            if self.nodes[node_id].order_dfs == 0 {
                self.set_scc(node_id, &mut visited_num, &mut stack);
            }
            debug_assert!(stack.is_empty());
        }
    }

    fn set_scc(
        &mut self,
        node_id: NodeId,
        visited_num: &mut usize,
        stack: &mut Vec<NodeId>,
    ) -> Option<NodeId> {
        if self.nodes[node_id].order_dfs != 0 {
            if self.nodes[node_id].in_stack {
                return self.nodes[node_id].scc;
            }
            return None;
        }

        *visited_num += 1;
        self.nodes[node_id].order_dfs = *visited_num;
        stack.push(node_id);
        self.nodes[node_id].in_stack = true;
        self.nodes[node_id].scc = Some(node_id);

        let out_arcs: Vec<usize> = self.nodes[node_id].out_arcs.clone();
        for arc_id in out_arcs {
            if !self.arcs[arc_id].active {
                continue;
            }
            let target = self.arcs[arc_id].to;
            if let Some(low_node) = self.set_scc(target, visited_num, stack) {
                let current_scc = self.nodes[node_id].scc.unwrap();
                if self.nodes[low_node].order_dfs < self.nodes[current_scc].order_dfs {
                    self.nodes[node_id].scc = Some(low_node);
                }
            }
        }

        let scc = self.nodes[node_id].scc.unwrap();
        if self.nodes[scc].order_dfs < self.nodes[node_id].order_dfs {
            return Some(scc);
        }

        loop {
            let n = stack.pop().unwrap();
            self.nodes[n].in_stack = false;
            if n == node_id {
                break;
            }
            self.nodes[n].scc = Some(node_id);
        }

        None
    }

    pub fn contract_scc(&mut self) {
        for node_id in 0..self.nodes.len() {
            if let Some(scc) = self.nodes[node_id].scc {
                if scc != node_id {
                    self.contract(scc, node_id);
                }
            }
        }
        for arc_id in 0..self.arcs.len() {
            if self.self_loop(arc_id) {
                self.delete_arc(arc_id);
            }
        }
    }

    pub fn reverse(old2new: &HashMap<NodeId, NodeId>) -> HashMap<NodeId, NodeId> {
        let mut new2old = HashMap::with_capacity(old2new.len());
        for (&old, &new_node) in old2new {
            new2old.insert(new_node, old);
        }
        new2old
    }

    pub fn borrow_arcs(
        &mut self,
        other: &DiGraph,
        other2this: &HashMap<NodeId, NodeId>,
        parallel_allowed: bool,
    ) {
        for (&other_node, &from) in other2this {
            for &arc_id in &other.nodes[other_node].out_arcs {
                let arc = &other.arcs[arc_id];
                if !arc.active {
                    continue;
                }
                if let Some(&to) = other2this.get(&arc.to) {
                    if parallel_allowed || !self.is_incident(from, to, true) {
                        self.add_arc(from, to);
                    }
                }
            }
        }
    }

    /// Clear reachability flags
    pub fn clear_reachable(&mut self) {
        for node in &mut self.nodes {
            node.reachable = false;
        }
    }

    /// Mark nodes reachable from the given node following arcs[false].
    pub fn set_reachable(&mut self, start: NodeId, reach_itself: bool) {
        if self.nodes[start].reachable {
            return;
        }
        if reach_itself {
            self.nodes[start].reachable = true;
        }
        let in_arcs: Vec<usize> = self.nodes[start].in_arcs.clone();
        for arc_id in in_arcs {
            if !self.arcs[arc_id].active {
                continue;
            }
            let source = self.arcs[arc_id].from;
            self.set_reachable(source, true);
        }
    }

    pub fn is_incident(&self, node: NodeId, n: NodeId, out: bool) -> bool {
        if out {
            for &arc_id in &self.nodes[node].out_arcs {
                if self.arcs[arc_id].active && self.arcs[arc_id].to == n {
                    return true;
                }
            }
        } else {
            for &arc_id in &self.nodes[node].in_arcs {
                if self.arcs[arc_id].active && self.arcs[arc_id].from == n {
                    return true;
                }
            }
        }
        false
    }

    pub fn is_incident_except(&self, node: NodeId, n: NodeId, out: bool) -> bool {
        if out {
            for &arc_id in &self.nodes[node].out_arcs {
                if self.arcs[arc_id].active && self.arcs[arc_id].to != n {
                    return true;
                }
            }
        } else {
            for &arc_id in &self.nodes[node].in_arcs {
                if self.arcs[arc_id].active && self.arcs[arc_id].from != n {
                    return true;
                }
            }
        }
        false
    }

    pub fn get_degree(&self, node: NodeId) -> usize {
        self.nodes[node]
            .in_arcs
            .iter()
            .filter(|&&arc_id| self.arcs[arc_id].active)
            .count()
            + self.nodes[node]
                .out_arcs
                .iter()
                .filter(|&&arc_id| self.arcs[arc_id].active)
                .count()
    }

    pub fn get_neighborhood(&self, node: NodeId, out: bool) -> Vec<NodeId> {
        let mut neighborhood = Vec::new();
        if out {
            neighborhood.reserve(self.nodes[node].out_arcs.len());
            for &arc_id in &self.nodes[node].out_arcs {
                if self.arcs[arc_id].active {
                    neighborhood.push(self.arcs[arc_id].to);
                }
            }
        } else {
            neighborhood.reserve(self.nodes[node].in_arcs.len());
            for &arc_id in &self.nodes[node].in_arcs {
                if self.arcs[arc_id].active {
                    neighborhood.push(self.arcs[arc_id].from);
                }
            }
        }
        neighborhood
    }

    pub fn get_neighborhood_all(&self, node: NodeId) -> Vec<NodeId> {
        let mut neighborhood = self.get_neighborhood(node, false);
        neighborhood.extend(self.get_neighborhood(node, true));
        neighborhood
    }

    pub fn get_children(&self, node: NodeId) -> Vec<NodeId> {
        self.get_neighborhood(node, false)
    }

    pub fn delete_neighborhood(&mut self, node: NodeId, out: bool) {
        let arcs = if out {
            self.nodes[node].out_arcs.clone()
        } else {
            self.nodes[node].in_arcs.clone()
        };
        for arc_id in arcs {
            self.delete_arc(arc_id);
        }
    }

    pub fn isolate(&mut self, node: NodeId) {
        self.delete_neighborhood(node, false);
        self.delete_neighborhood(node, true);
    }

    pub fn contract(&mut self, node: NodeId, from: NodeId) {
        if node == from {
            return;
        }

        for out in [false, true] {
            let node_arcs = if out {
                self.nodes[node].out_arcs.clone()
            } else {
                self.nodes[node].in_arcs.clone()
            };
            let from_arcs = if out {
                self.nodes[from].out_arcs.clone()
            } else {
                self.nodes[from].in_arcs.clone()
            };

            for arc_id in from_arcs {
                if !self.arcs[arc_id].active {
                    continue;
                }
                let neighbor = if out {
                    self.arcs[arc_id].to
                } else {
                    self.arcs[arc_id].from
                };
                let mut parallel = false;
                for node_arc_id in &node_arcs {
                    if !self.arcs[*node_arc_id].active {
                        continue;
                    }
                    let node_neighbor = if out {
                        self.arcs[*node_arc_id].to
                    } else {
                        self.arcs[*node_arc_id].from
                    };
                    if node_neighbor == neighbor {
                        parallel = true;
                        break;
                    }
                }
                if parallel {
                    self.delete_arc(arc_id);
                } else {
                    self.set_node(arc_id, node, !out);
                }
            }
        }
        self.isolate(from);
    }

    pub fn set_node(&mut self, arc_id: usize, new_node: NodeId, out: bool) {
        if !self.arcs[arc_id].active {
            return;
        }
        if out {
            let old_node = self.arcs[arc_id].to;
            if old_node == new_node {
                return;
            }
            self.nodes[old_node].in_arcs.retain(|&id| id != arc_id);
            self.arcs[arc_id].to = new_node;
            self.nodes[new_node].in_arcs.push(arc_id);
        } else {
            let old_node = self.arcs[arc_id].from;
            if old_node == new_node {
                return;
            }
            self.nodes[old_node].out_arcs.retain(|&id| id != arc_id);
            self.arcs[arc_id].from = new_node;
            self.nodes[new_node].out_arcs.push(arc_id);
        }
    }

    pub fn self_loop(&self, arc_id: usize) -> bool {
        self.arcs[arc_id].active && self.arcs[arc_id].from == self.arcs[arc_id].to
    }

    pub fn delete_arc(&mut self, arc_id: usize) {
        if !self.arcs[arc_id].active {
            return;
        }
        let from = self.arcs[arc_id].from;
        let to = self.arcs[arc_id].to;
        self.nodes[from].out_arcs.retain(|&id| id != arc_id);
        self.nodes[to].in_arcs.retain(|&id| id != arc_id);
        self.arcs[arc_id].active = false;
    }

    pub fn get_ends(&self, out: bool) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|n| {
                if out {
                    n.out_arcs.iter().all(|&arc_id| !self.arcs[arc_id].active)
                } else {
                    n.in_arcs.iter().all(|&arc_id| !self.arcs[arc_id].active)
                }
            })
            .map(|n| n.id)
            .collect()
    }

    pub fn get_root(&self, out: bool) -> Option<NodeId> {
        let ends = self.get_ends(out);
        if ends.len() == 1 {
            Some(ends[0])
        } else {
            None
        }
    }

    /// Get root nodes (nodes with no incoming arcs)
    pub fn roots(&self) -> Vec<NodeId> {
        self.get_ends(false)
    }

    /// Get leaf nodes (nodes with no outgoing arcs)
    pub fn leaves(&self) -> Vec<NodeId> {
        self.get_ends(true)
    }
}

impl Default for DiGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl Node {
    pub fn get_disjoint_cluster(&self) -> NodeId {
        self.parent_dc
    }

    pub fn qc(&self, graph: &DiGraph) -> Result<(), String> {
        if self.id >= graph.nodes.len() {
            return Err("node id is outside graph".to_string());
        }
        if self.parent_dc >= graph.nodes.len() {
            return Err("disjoint-cluster parent is outside graph".to_string());
        }
        Ok(())
    }

    pub fn save_text(&self, graph: &DiGraph) -> String {
        graph.node_save_text(self.id)
    }
}

/// Tree node index type.
pub type TreeNodeId = usize;

/// Distance-bearing tree node used by `Tree`.
#[derive(Clone, Debug)]
pub struct TreeNode {
    pub id: TreeNodeId,
    pub name: String,
    pub parent: Option<TreeNodeId>,
    pub children: Vec<TreeNodeId>,
    pub parent_distance: f64,
    pub frequent_child: bool,
    pub frequent_degree: usize,
    pub leaves: usize,
    pub detached: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TipName {
    pub name: String,
    pub depth: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct NodeDist {
    pub node: TreeNodeId,
    pub dist: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Patristic {
    pub leaf1: TreeNodeId,
    pub leaf2: TreeNodeId,
    pub distance: f64,
}

#[derive(Clone, Debug, Default)]
pub struct LcaBuffer {
    pub vec1: Vec<TreeNodeId>,
    pub vec2: Vec<TreeNodeId>,
}

impl LcaBuffer {
    pub fn clear(&mut self) {
        self.vec1.clear();
        self.vec2.clear();
    }
}

/// Rooted tree. Parent arcs follow the original `DiGraph` convention:
/// `out = false` lists children and `out = true` lists the parent.
#[derive(Clone, Debug, Default)]
pub struct Tree {
    pub nodes: Vec<TreeNode>,
    pub root: Option<TreeNodeId>,
}

impl Tree {
    pub const OBJ_NAME_SEPARATOR: char = ':';

    pub fn new() -> Self {
        Self::default()
    }

    pub fn empty(&self) -> bool {
        self.nodes.iter().all(|n| n.detached)
    }

    pub fn clear(&mut self) {
        self.nodes.clear();
        self.root = None;
    }

    pub fn add_node(&mut self, name: impl Into<String>, parent: Option<TreeNodeId>) -> TreeNodeId {
        let id = self.nodes.len();
        self.nodes.push(TreeNode {
            id,
            name: name.into(),
            parent: None,
            children: Vec::new(),
            parent_distance: -1.0,
            frequent_child: false,
            frequent_degree: 0,
            leaves: 0,
            detached: false,
        });
        self.set_parent(id, parent);
        id
    }

    pub fn qc(&self) -> Result<(), String> {
        let live: Vec<_> = self.nodes.iter().filter(|n| !n.detached).collect();
        if live.is_empty() {
            if self.root.is_some() {
                return Err("empty tree has root".to_string());
            }
            return Ok(());
        }
        let root = self
            .root
            .ok_or_else(|| "non-empty tree has no root".to_string())?;
        if root >= self.nodes.len() || self.nodes[root].detached {
            return Err("bad root".to_string());
        }
        let mut roots = 0;
        let mut names = HashSet::new();
        for node in live {
            if node.name.contains(Self::OBJ_NAME_SEPARATOR) {
                return Err(format!(
                    "{:?} contains separator {:?}",
                    node.name,
                    Self::OBJ_NAME_SEPARATOR
                ));
            }
            if node.is_leaf() && !names.insert(node.name.clone()) {
                return Err(format!("Duplicate name: {}", node.name));
            }
            match node.parent {
                Some(parent) => {
                    if parent >= self.nodes.len() || self.nodes[parent].detached {
                        return Err(format!("node {} has bad parent", node.id));
                    }
                    if !self.nodes[parent].children.contains(&node.id) {
                        return Err(format!("node {} missing from parent children", node.id));
                    }
                }
                None => roots += 1,
            }
            for &child in &node.children {
                if child >= self.nodes.len() || self.nodes[child].detached {
                    return Err(format!("node {} has bad child", node.id));
                }
                if self.nodes[child].parent != Some(node.id) {
                    return Err(format!("child {child} parent mismatch"));
                }
            }
        }
        if roots != 1 {
            return Err(format!("expected one root, found {roots}"));
        }
        if self.nodes[root].parent.is_some() {
            return Err("root has parent".to_string());
        }
        Ok(())
    }

    pub fn save_text(&self) -> String {
        let mut s = String::new();
        if let Some(root) = self.root {
            self.save_node_text(root, &mut s, 0);
        }
        s.push('\n');
        s
    }

    fn save_node_text(&self, node: TreeNodeId, s: &mut String, depth: usize) {
        if self.nodes[node].detached {
            return;
        }
        let indent = "  ".repeat(depth);
        let _ = write!(s, "{indent}{}: ", self.nodes[node].name);
        if self.nodes[node].children.is_empty() {
            let _ = writeln!(s);
        } else {
            let _ = writeln!(s);
            for &child in &self.nodes[node].children {
                self.save_node_text(child, s, depth + 1);
            }
        }
    }

    pub fn print_newick(&self, internal_names: bool, minimal_leaf_name: bool) -> String {
        let mut s = String::new();
        if let Some(root) = self.root {
            self.print_newick_node(root, internal_names, minimal_leaf_name, &mut s);
        }
        s.push_str(";\n");
        s
    }

    fn print_newick_node(
        &self,
        node: TreeNodeId,
        internal_names: bool,
        minimal_leaf_name: bool,
        s: &mut String,
    ) {
        let n = &self.nodes[node];
        if n.is_leaf() {
            s.push_str(&TreeNode::name2newick(
                &n.get_newick_name(minimal_leaf_name),
            ));
        } else {
            s.push('(');
            for (i, &child) in n.children.iter().enumerate() {
                if i != 0 {
                    s.push(',');
                }
                self.print_newick_node(child, internal_names, minimal_leaf_name, s);
            }
            s.push(')');
            if internal_names {
                s.push_str(&TreeNode::name2newick(
                    &n.get_newick_name(minimal_leaf_name),
                ));
            }
        }
        let dist = n.get_parent_distance();
        if !dist.is_nan() && dist != -1.0 {
            let _ = write!(s, ":{dist:.6}");
        }
    }

    pub fn print_arc_lengths_columns() -> &'static str {
        "<node name> <arc length> <depth length> <log(<parent arc length>/<arc length>)"
    }

    pub fn print_arc_lengths(&self) -> String {
        let mut s = String::new();
        for node in &self.nodes {
            if node.detached || Some(node.id) == self.root {
                continue;
            }
            let dist = node.get_parent_distance();
            if !dist.is_nan() && dist != -1.0 && dist > 0.0 {
                let _ = write!(
                    s,
                    "{} {} {}",
                    self.get_lca_name(node.id),
                    dist,
                    self.get_root_distance(node.id)
                );
                if let Some(parent) = node.parent {
                    if Some(parent) != self.root {
                        let dist_parent = self.nodes[parent].get_parent_distance();
                        if !dist_parent.is_nan() && dist_parent != -1.0 {
                            let _ = write!(s, " {}", dist_parent.ln() - dist.ln());
                        } else {
                            s.push_str(" ?");
                        }
                    } else {
                        s.push_str(" ?");
                    }
                } else {
                    s.push_str(" ?");
                }
                s.push('\n');
            }
        }
        s
    }

    pub fn get_length(&self) -> f64 {
        self.root.map_or(0.0, |root| self.get_subtree_length(root))
    }

    pub fn get_ave_arc_length(&self) -> f64 {
        let mut len = 0.0;
        let mut n = 0;
        for node in &self.nodes {
            if node.detached || Some(node.id) == self.root {
                continue;
            }
            len += node.get_parent_distance();
            n += 1;
        }
        len / n as f64
    }

    pub fn set_parent(&mut self, node: TreeNodeId, new_parent: Option<TreeNodeId>) {
        assert_ne!(Some(node), new_parent);
        if let Some(old_parent) = self.nodes[node].parent {
            self.nodes[old_parent]
                .children
                .retain(|&child| child != node);
        }
        self.nodes[node].parent = new_parent;
        if let Some(parent) = new_parent {
            self.nodes[parent].children.push(node);
        } else {
            self.root = Some(node);
        }
    }

    pub fn get_parent(&self, node: TreeNodeId) -> Option<TreeNodeId> {
        self.nodes[node].parent
    }

    pub fn get_children(&self, node: TreeNodeId) -> Vec<TreeNodeId> {
        self.nodes[node].children.clone()
    }

    pub fn get_ancestor(&self, mut node: TreeNodeId, height: usize) -> TreeNodeId {
        for _ in 0..height {
            if let Some(parent) = self.nodes[node].parent {
                node = parent;
            } else {
                return self.root.expect("tree root is required");
            }
        }
        node
    }

    pub fn get_tip_name(&self, node: TreeNodeId) -> TipName {
        if self.nodes[node].is_leaf() {
            return TipName {
                name: self.nodes[node].name.clone(),
                depth: 0,
            };
        }
        let mut best: Option<TipName> = None;
        for &child in &self.nodes[node].children {
            let tn = self.get_tip_name(child);
            if best.as_ref().is_none_or(|b| b.name > tn.name) {
                best = Some(tn);
            }
        }
        let mut tn = best.unwrap();
        tn.depth += 1;
        tn
    }

    pub fn get_topological_depth(&self, node: TreeNodeId) -> usize {
        self.nodes[node]
            .parent
            .map_or(0, |parent| self.get_topological_depth(parent) + 1)
    }

    pub fn get_subtree_heights(&self, node: TreeNodeId, node_heights: &mut Vec<NodeDist>) {
        if self.nodes[node].is_leaf() {
            return;
        }
        let mut height = 0.0;
        for &child in &self.nodes[node].children {
            let dist = self.nodes[child].get_parent_distance();
            self.get_subtree_heights(child, node_heights);
            let child_height = if self.nodes[child].is_leaf() {
                0.0
            } else {
                node_heights.last().map_or(0.0, |nd| nd.dist)
            };
            height = f64::max(height, dist + child_height);
        }
        node_heights.push(NodeDist { node, dist: height });
    }

    pub fn get_leaf_depths(&self, node: TreeNodeId, leaf_depths: &mut Vec<NodeDist>) {
        self.get_leaf_depths_(node, leaf_depths, true);
    }

    fn get_leaf_depths_(&self, node: TreeNodeId, leaf_depths: &mut Vec<NodeDist>, first: bool) {
        let start = leaf_depths.len();
        if self.nodes[node].is_leaf() {
            leaf_depths.push(NodeDist { node, dist: 0.0 });
        } else {
            for &child in &self.nodes[node].children {
                self.get_leaf_depths_(child, leaf_depths, false);
            }
        }
        if !first {
            let parent_distance = self.nodes[node].get_parent_distance();
            for nd in &mut leaf_depths[start..] {
                nd.dist += parent_distance;
            }
        }
    }

    pub fn get_height(&self, node: TreeNodeId) -> usize {
        self.nodes[node]
            .children
            .iter()
            .map(|&child| 1 + self.get_height(child))
            .max()
            .unwrap_or(0)
    }

    pub fn get_interior_height(&self) -> usize {
        self.root
            .filter(|&root| self.nodes[root].is_interior_type())
            .map_or(0, |root| self.get_interior_height_node(root))
    }

    pub fn get_interior_height_node(&self, node: TreeNodeId) -> usize {
        self.nodes[node]
            .children
            .iter()
            .filter(|&&child| self.nodes[child].is_interior_type())
            .map(|&child| 1 + self.get_interior_height_node(child))
            .max()
            .unwrap_or(0)
    }

    pub fn get_distance_height(&self, node: TreeNodeId) -> f64 {
        self.nodes[node]
            .children
            .iter()
            .map(|&child| self.nodes[child].get_parent_distance() + self.get_distance_height(child))
            .fold(0.0, f64::max)
    }

    pub fn get_root_distance(&self, node: TreeNodeId) -> f64 {
        self.nodes[node].parent.map_or(0.0, |parent| {
            self.nodes[node].get_parent_distance() + self.get_root_distance(parent)
        })
    }

    pub fn descendant_of(&self, node: TreeNodeId, ancestor: Option<TreeNodeId>) -> bool {
        let Some(ancestor) = ancestor else {
            return true;
        };
        if node == ancestor {
            return true;
        }
        self.nodes[node]
            .parent
            .is_some_and(|parent| self.descendant_of(parent, Some(ancestor)))
    }

    pub fn get_prev_ancestor(
        &self,
        node: TreeNodeId,
        ancestor: Option<TreeNodeId>,
    ) -> Option<TreeNodeId> {
        let ancestor = ancestor?;
        let parent = self.nodes[node].parent?;
        if parent == ancestor {
            Some(node)
        } else {
            self.get_prev_ancestor(parent, Some(ancestor))
        }
    }

    pub fn get_path_length(&self, node: TreeNodeId, ancestor: Option<TreeNodeId>) -> f64 {
        if Some(node) == ancestor {
            return 0.0;
        }
        if let Some(parent) = self.nodes[node].parent {
            return self.nodes[node].get_parent_distance() + self.get_path_length(parent, ancestor);
        }
        if ancestor.is_none() {
            0.0
        } else {
            f64::NAN
        }
    }

    pub fn get_subtree_size(&self, node: TreeNodeId, count_leaves: bool) -> usize {
        let mut n = 0;
        for &child in &self.nodes[node].children {
            if !count_leaves && self.nodes[child].children.is_empty() {
                continue;
            }
            n += 1 + self.get_subtree_size(child, count_leaves);
        }
        n
    }

    pub fn get_subtree_length(&self, node: TreeNodeId) -> f64 {
        self.nodes[node]
            .children
            .iter()
            .map(|&child| self.nodes[child].get_parent_distance() + self.get_subtree_length(child))
            .sum()
    }

    pub fn subtree_size2leaves(&mut self, node: TreeNodeId) {
        self.nodes[node].leaves = 0;
        let children = self.nodes[node].children.clone();
        for child in children {
            self.subtree_size2leaves(child);
            self.nodes[node].leaves += 1 + self.nodes[child].leaves;
        }
    }

    pub fn set_leaves(&mut self) {
        if let Some(root) = self.root {
            self.set_leaves_node(root);
        }
    }

    pub fn set_leaves_node(&mut self, node: TreeNodeId) {
        self.nodes[node].leaves = 0;
        let children = self.nodes[node].children.clone();
        for child in children {
            self.set_leaves_node(child);
            self.nodes[node].leaves += self.nodes[child].leaves;
        }
        if self.nodes[node].leaves == 0 {
            self.nodes[node].leaves = 1;
        }
    }

    pub fn set_leaves_value(&mut self, node: TreeNodeId, leaves: usize) {
        self.nodes[node].leaves = leaves;
        let children = self.nodes[node].children.clone();
        for child in children {
            self.set_leaves_value(child, leaves);
        }
    }

    pub fn get_leaves_size(&self, node: TreeNodeId) -> usize {
        usize::max(
            self.nodes[node]
                .children
                .iter()
                .map(|&child| self.get_leaves_size(child))
                .sum(),
            1,
        )
    }

    pub fn children2frequent_child(&mut self, node: TreeNodeId, rare_prob: f64) {
        assert!((0.0..0.5).contains(&rare_prob));
        assert!(self.nodes[node].frequent_child);
        let mut child_freqs: Vec<_> = self.nodes[node]
            .children
            .iter()
            .map(|&child| (child, self.nodes[child].leaves))
            .collect();
        child_freqs.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        let mut sum = 0usize;
        for (child, freq) in child_freqs {
            sum += freq;
            if (freq as f64) / (sum as f64) < rare_prob {
                continue;
            }
            self.nodes[child].frequent_child = true;
            self.children2frequent_child(child, rare_prob);
        }
    }

    pub fn get_leaves(&self, node: TreeNodeId, leaf_vec: &mut Vec<TreeNodeId>) {
        if self.nodes[node].children.is_empty() {
            leaf_vec.push(node);
        } else {
            for &child in &self.nodes[node].children {
                self.get_leaves(child, leaf_vec);
            }
        }
    }

    pub fn get_closest_leaf(&self, node: TreeNodeId, leaf_depth: &mut usize) -> TreeNodeId {
        let mut leaf = None;
        *leaf_depth = usize::MAX;
        for &child in &self.nodes[node].children {
            let mut depth1 = 0;
            let leaf1 = self.get_closest_leaf(child, &mut depth1);
            if depth1 < *leaf_depth {
                *leaf_depth = depth1;
                leaf = Some(leaf1);
            }
        }
        if let Some(leaf) = leaf {
            *leaf_depth += 1;
            leaf
        } else {
            *leaf_depth = 0;
            node
        }
    }

    pub fn get_other_child(&self, node: TreeNodeId, child: TreeNodeId) -> Option<TreeNodeId> {
        debug_assert_eq!(self.nodes[child].parent, Some(node));
        let mut other = None;
        for &n in &self.nodes[node].children {
            if n != child {
                debug_assert!(other.is_none());
                other = Some(n);
            }
        }
        other
    }

    pub fn get_different_child(&self, node: TreeNodeId, child: TreeNodeId) -> TreeNodeId {
        self.nodes[node]
            .children
            .iter()
            .copied()
            .find(|&n| n != child)
            .expect("getDifferentChild(): Transient parent")
    }

    pub fn get_first_decendant(&self, mut node: TreeNodeId) -> TreeNodeId {
        while let Some(&child) = self.nodes[node].children.first() {
            node = child;
        }
        node
    }

    pub fn get_last_decendant(&self, mut node: TreeNodeId) -> TreeNodeId {
        while let Some(&child) = self.nodes[node].children.last() {
            node = child;
        }
        node
    }

    pub fn get_lca_name(&self, node: TreeNodeId) -> String {
        let left = self.get_first_decendant(node);
        let right = self.get_last_decendant(node);
        if left == right {
            self.nodes[left].get_leaf_name()
        } else {
            format!(
                "{}{}{}",
                self.nodes[left].get_leaf_name(),
                Self::OBJ_NAME_SEPARATOR,
                self.nodes[right].get_leaf_name()
            )
        }
    }

    pub fn children_up(&mut self, node: TreeNodeId) {
        let parent = self.nodes[node].parent;
        let children = self.nodes[node].children.clone();
        for child in children {
            self.set_parent(child, parent);
        }
        self.nodes[node].children.clear();
    }

    pub fn detach_children_up(&mut self, node: TreeNodeId) {
        self.children_up(node);
        if let Some(parent) = self.nodes[node].parent {
            self.nodes[parent].children.retain(|&child| child != node);
        }
        self.nodes[node].parent = None;
        self.nodes[node].detached = true;
        if self.root == Some(node) {
            self.root = None;
        }
    }

    pub fn is_transient(&self, node: TreeNodeId) -> Option<TreeNodeId> {
        (self.nodes[node].children.len() == 1).then_some(self.nodes[node].children[0])
    }

    pub fn isolate_transient(&mut self, node: TreeNodeId) -> Option<TreeNodeId> {
        let transient = self.is_transient(node);
        if transient.is_some() {
            self.detach_children_up(node);
        }
        transient
    }

    pub fn delete_transient_node(&mut self, node: TreeNodeId) -> bool {
        if self.isolate_transient(node).is_none() {
            return false;
        }
        self.nodes[node].detached = true;
        true
    }

    pub fn delete_subtree(&mut self, node: TreeNodeId) {
        let children = self.nodes[node].children.clone();
        for child in children {
            self.delete_subtree(child);
            self.nodes[child].detached = true;
        }
        self.nodes[node].children.clear();
    }

    pub fn make_root(&mut self, node: TreeNodeId) -> TreeNodeId {
        let root_old = self.root.expect("tree root is required");
        if node == root_old {
            return node;
        }
        let mut parent_new = None;
        let mut current = Some(node);
        while let Some(n) = current {
            let parent_old = self.nodes[n].parent;
            self.set_parent(n, parent_new);
            parent_new = Some(n);
            current = parent_old;
        }
        root_old
    }

    pub fn get_area(
        &self,
        node: TreeNodeId,
        radius: u32,
        area: &mut Vec<TreeNodeId>,
        boundary: &mut Vec<TreeNodeId>,
    ) {
        self.get_area_(node, radius, None, area, boundary);
    }

    fn get_area_(
        &self,
        node: TreeNodeId,
        radius: u32,
        prev: Option<TreeNodeId>,
        area: &mut Vec<TreeNodeId>,
        boundary: &mut Vec<TreeNodeId>,
    ) {
        area.push(node);
        let mut degree = usize::from(prev.is_some());
        if radius != 0 {
            if let Some(parent) = self.nodes[node].parent {
                if Some(parent) != prev {
                    self.get_area_(parent, radius - 1, Some(node), area, boundary);
                    degree += 1;
                }
            }
            for &child in &self.nodes[node].children {
                if Some(child) != prev {
                    self.get_area_(child, radius - 1, Some(node), area, boundary);
                    degree += 1;
                }
            }
        }
        if degree <= 1 {
            boundary.push(node);
        }
    }

    pub fn get_distance_area(
        &self,
        node: TreeNodeId,
        radius: f64,
        area: &mut Vec<TreeNodeId>,
        boundary: &mut Vec<TreeNodeId>,
    ) {
        self.get_distance_area_(node, radius, None, area, boundary);
    }

    fn get_distance_area_(
        &self,
        node: TreeNodeId,
        radius: f64,
        prev: Option<TreeNodeId>,
        area: &mut Vec<TreeNodeId>,
        boundary: &mut Vec<TreeNodeId>,
    ) {
        if radius < 0.0 {
            return;
        }
        area.push(node);
        let mut degree = usize::from(prev.is_some());
        if let Some(parent) = self.nodes[node].parent {
            if Some(parent) != prev {
                self.get_distance_area_(
                    parent,
                    radius - self.nodes[node].get_parent_distance(),
                    Some(node),
                    area,
                    boundary,
                );
                degree += 1;
            }
        }
        for &child in &self.nodes[node].children {
            if Some(child) != prev {
                self.get_distance_area_(
                    child,
                    radius - self.nodes[child].get_parent_distance(),
                    Some(node),
                    area,
                    boundary,
                );
                degree += 1;
            }
        }
        if degree <= 1 {
            boundary.push(node);
        }
    }

    pub fn get_closest_leaves(
        &self,
        node: TreeNodeId,
        neighbors_max: usize,
        neighbors: &mut Vec<NodeDist>,
    ) {
        self.get_closest_leaves_(node, None, 0.0, neighbors_max, neighbors);
    }

    fn get_closest_leaves_(
        &self,
        node: TreeNodeId,
        prev: Option<TreeNodeId>,
        distance: f64,
        neighbors_max: usize,
        neighbors: &mut Vec<NodeDist>,
    ) {
        assert!(neighbors_max != 0);
        if neighbors.len() == neighbors_max && neighbors.last().is_some_and(|n| n.dist <= distance)
        {
            return;
        }
        if prev.is_some() && self.nodes[node].is_leaf() {
            neighbors.push(NodeDist {
                node,
                dist: distance,
            });
            neighbors.sort_by(|a, b| a.dist.total_cmp(&b.dist));
            if neighbors.len() > neighbors_max {
                neighbors.pop();
            }
        }
        if let Some(parent) = self.nodes[node].parent {
            if Some(parent) != prev {
                self.get_closest_leaves_(
                    parent,
                    Some(node),
                    distance + self.nodes[node].get_parent_distance(),
                    neighbors_max,
                    neighbors,
                );
            }
        }
        for &child in &self.nodes[node].children {
            if Some(child) != prev {
                self.get_closest_leaves_(
                    child,
                    Some(node),
                    distance + self.nodes[child].get_parent_distance(),
                    neighbors_max,
                    neighbors,
                );
            }
        }
    }

    pub fn get_subtree_area(
        &self,
        node: TreeNodeId,
        possible_boundary: &[TreeNodeId],
        area: &mut Vec<TreeNodeId>,
        boundary: &mut Vec<TreeNodeId>,
    ) {
        area.push(node);
        if possible_boundary.contains(&node) || self.nodes[node].children.is_empty() {
            boundary.push(node);
        } else {
            for &child in &self.nodes[node].children {
                self.get_subtree_area(child, possible_boundary, area, boundary);
            }
        }
    }

    pub fn get_leaf_distances(&self) -> Vec<Patristic> {
        let mut leaf2dist = HashMap::new();
        let Some(root) = self.root else {
            return Vec::new();
        };
        self.node2leaf_distances(root, &mut leaf2dist)
    }

    fn node2leaf_distances(
        &self,
        node: TreeNodeId,
        leaf2dist: &mut HashMap<TreeNodeId, f64>,
    ) -> Vec<Patristic> {
        let mut res = Vec::new();
        if self.nodes[node].is_leaf() {
            leaf2dist.insert(node, 0.0);
        } else {
            for &child in &self.nodes[node].children {
                let mut node_leaf2dist = HashMap::new();
                res.extend(self.node2leaf_distances(child, &mut node_leaf2dist));
                let dist = self.nodes[child].get_parent_distance();
                for d in node_leaf2dist.values_mut() {
                    *d += dist;
                }
                for (&leaf1, &dist1) in leaf2dist.iter() {
                    for (&leaf2, &dist2) in node_leaf2dist.iter() {
                        res.push(self.patristic(leaf1, leaf2, dist1 + dist2));
                    }
                }
                leaf2dist.extend(node_leaf2dist);
            }
        }
        res
    }

    fn patristic(&self, leaf1: TreeNodeId, leaf2: TreeNodeId, distance: f64) -> Patristic {
        assert_ne!(self.nodes[leaf1].name, self.nodes[leaf2].name);
        if self.nodes[leaf1].name > self.nodes[leaf2].name {
            Patristic {
                leaf1: leaf2,
                leaf2: leaf1,
                distance,
            }
        } else {
            Patristic {
                leaf1,
                leaf2,
                distance,
            }
        }
    }

    pub fn size(&self, count_leaves: bool) -> usize {
        let live = self.nodes.iter().filter(|n| !n.detached).count();
        if live <= 1 {
            usize::from(count_leaves && live == 1)
        } else {
            1 + self.get_subtree_size(self.root.expect("tree root is required"), count_leaves)
        }
    }

    pub fn count_interior_nodes(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| !n.detached && n.is_interior_type())
            .count()
    }

    pub fn is_star(&self) -> bool {
        self.count_interior_nodes() == 1
    }

    pub fn radius2boundary_size(radius: u32) -> usize {
        if radius == 0 {
            1
        } else {
            3 * 2usize.pow(radius - 1)
        }
    }

    pub fn get_bifurcating_interior_branching(&self) -> f64 {
        let Some(root) = self.root else {
            return -1.0;
        };
        if !self.nodes[root].is_interior_type() {
            return -1.0;
        }
        let mut bifurcating_interior_nodes = 0;
        let mut branches = 0;
        self.get_bifurcating_interior_branching_node(
            root,
            &mut bifurcating_interior_nodes,
            &mut branches,
        );
        branches as f64 / bifurcating_interior_nodes as f64
    }

    fn get_bifurcating_interior_branching_node(
        &self,
        node: TreeNodeId,
        bifurcating_interior_nodes: &mut usize,
        branches: &mut usize,
    ) {
        let child_count = self.nodes[node].children.len();
        assert!(child_count >= 2);
        *branches += child_count - 2;
        *bifurcating_interior_nodes += child_count - 2;
        let mut interior_arcs = 0;
        for &child in &self.nodes[node].children {
            if self.nodes[child].is_interior_type() {
                interior_arcs += 1;
                self.get_bifurcating_interior_branching_node(
                    child,
                    bifurcating_interior_nodes,
                    branches,
                );
            }
        }
        *branches += interior_arcs;
        if interior_arcs != 0 {
            *bifurcating_interior_nodes += 1;
        }
    }

    pub fn count_interior_undirected_arcs(&self) -> usize {
        let Some(root) = self.root else {
            return 0;
        };
        let mut n = self.count_interior_nodes();
        if self.nodes[root].is_interior_type() {
            n = n.saturating_sub(1);
        }
        match self.nodes[root].children.as_slice() {
            [child] if self.nodes[*child].is_interior_type() => n = n.saturating_sub(1),
            [a, b] if self.nodes[*a].is_interior_type() || self.nodes[*b].is_interior_type() => {
                n = n.saturating_sub(1)
            }
            _ => {}
        }
        n
    }

    pub fn get_lca(
        &self,
        n1: Option<TreeNodeId>,
        n2: Option<TreeNodeId>,
        buf: &mut LcaBuffer,
    ) -> Option<TreeNodeId> {
        let (n1, n2) = (n1?, n2?);
        if n1 == n2 {
            return Some(n1);
        }
        buf.clear();
        if self.get_parents_or_target(n1, n2, &mut buf.vec1) {
            return Some(n2);
        }
        if self.get_parents_or_target(n2, n1, &mut buf.vec2) {
            return Some(n1);
        }
        let mut m = None;
        let mut i1 = buf.vec1.len() - 1;
        let mut i2 = buf.vec2.len() - 1;
        while buf.vec1[i1] == buf.vec2[i2] {
            m = Some(buf.vec1[i1]);
            if i1 == 0 || i2 == 0 {
                break;
            }
            i1 -= 1;
            i2 -= 1;
        }
        m
    }

    fn get_parents_or_target(
        &self,
        mut from: TreeNodeId,
        target: TreeNodeId,
        parents: &mut Vec<TreeNodeId>,
    ) -> bool {
        while Some(from) != self.nodes[from].parent {
            if from == target {
                return true;
            }
            parents.push(from);
            let Some(parent) = self.nodes[from].parent else {
                return false;
            };
            from = parent;
        }
        false
    }

    pub fn get_lca_vec(&self, node_vec: &[TreeNodeId], buf: &mut LcaBuffer) -> Option<TreeNodeId> {
        let mut n = *node_vec.first()?;
        for &node in &node_vec[1..] {
            n = self.get_lca(Some(n), Some(node), buf)?;
        }
        Some(n)
    }

    pub fn get_parents(
        &self,
        node_vec: &[TreeNodeId],
        buf: &mut LcaBuffer,
    ) -> BTreeSet<TreeNodeId> {
        let mut s = BTreeSet::new();
        let lca = self.get_lca_vec(node_vec, buf);
        for &node in node_vec {
            let mut n = Some(node);
            while n != lca {
                let current = n.expect("node has no parent before LCA");
                s.insert(current);
                n = self.nodes[current].parent;
            }
        }
        s
    }

    pub fn get_path(
        &self,
        n1: TreeNodeId,
        n2: TreeNodeId,
        ca: Option<TreeNodeId>,
        lca: &mut Option<TreeNodeId>,
        buf: &mut LcaBuffer,
    ) -> Vec<TreeNodeId> {
        buf.clear();
        if n1 == n2 {
            *lca = Some(n1);
            return Vec::new();
        }
        let n1_init = n1;
        let mut cur1 = Some(n1);
        while cur1 != ca && cur1 != Some(n2) {
            let node = cur1.expect("first node path misses common ancestor");
            buf.vec1.push(node);
            cur1 = self.nodes[node].parent;
        }
        if cur1 == Some(n2) {
            *lca = Some(n2);
            return buf.vec1.clone();
        }
        let mut cur2 = Some(n2);
        while cur2 != ca && cur2 != Some(n1_init) {
            let node = cur2.expect("second node path misses common ancestor");
            buf.vec2.push(node);
            cur2 = self.nodes[node].parent;
        }
        if cur2 == Some(n1_init) {
            *lca = Some(n1_init);
            return buf.vec2.clone();
        }
        *lca = ca;
        let mut i1 = buf.vec1.len() - 1;
        let mut i2 = buf.vec2.len() - 1;
        while buf.vec1[i1] == buf.vec2[i2] {
            *lca = Some(buf.vec1[i1]);
            buf.vec1.pop();
            buf.vec2.pop();
            if i1 == 0 || i2 == 0 {
                break;
            }
            i1 -= 1;
            i2 -= 1;
        }
        if buf.vec1.len() < buf.vec2.len() {
            let mut path = buf.vec2.clone();
            path.extend(buf.vec1.iter().rev().copied());
            path
        } else {
            let mut path = buf.vec1.clone();
            path.extend(buf.vec2.iter().rev().copied());
            path
        }
    }

    pub fn leaves2lcas(&mut self, leaves: &BTreeSet<TreeNodeId>) -> Vec<TreeNodeId> {
        let Some(root) = self.root else {
            return Vec::new();
        };
        let mut lcas = Vec::new();
        if leaves.is_empty() {
            return lcas;
        }
        self.set_leaves_node(root);
        let mut leaves2nodes = vec![BTreeSet::new(); self.nodes[root].leaves];
        for node in &self.nodes {
            if node.detached || node.leaves == 0 {
                continue;
            }
            leaves2nodes[node.leaves - 1].insert(node.id);
        }
        for i in 0..leaves2nodes.len() {
            let tree_nodes: Vec<_> = leaves2nodes[i].iter().copied().collect();
            let mut bads = Vec::new();
            for parent in tree_nodes {
                let mut good = true;
                let children = self.nodes[parent].children.clone();
                if children.is_empty() {
                    if !leaves.contains(&parent) {
                        good = false;
                    }
                } else {
                    for &child in &children {
                        if !leaves2nodes[self.nodes[child].leaves - 1].contains(&child) {
                            good = false;
                            break;
                        }
                    }
                }
                if good {
                    for child in children {
                        let child_idx = self.nodes[child].leaves - 1;
                        leaves2nodes[child_idx].remove(&child);
                    }
                } else {
                    bads.push(parent);
                }
            }
            for bad in bads {
                leaves2nodes[i].remove(&bad);
            }
        }
        for tree_nodes in leaves2nodes {
            lcas.extend(tree_nodes);
        }
        lcas
    }

    pub fn set_frequent_child(&mut self, rare_prob: f64) {
        let Some(root) = self.root else {
            return;
        };
        self.set_leaves();
        for node in &mut self.nodes {
            node.frequent_child = false;
        }
        self.nodes[root].frequent_child = true;
        self.children2frequent_child(root, rare_prob);
        let transients: Vec<_> = self
            .nodes
            .iter()
            .filter(|tn| tn.frequent_child)
            .filter(|tn| {
                tn.children
                    .iter()
                    .filter(|&&child| self.nodes[child].frequent_child)
                    .count()
                    == 1
            })
            .map(|tn| tn.id)
            .collect();
        for tn in transients {
            self.nodes[tn].frequent_child = false;
        }
    }

    pub fn set_frequent_degree(&mut self, rare_prob: f64) {
        assert!((0.0..0.3).contains(&rare_prob));
        let Some(root) = self.root else {
            return;
        };
        self.set_leaves();
        let all_leaves = self.nodes[root].leaves;
        for node_id in 0..self.nodes.len() {
            if self.nodes[node_id].detached {
                continue;
            }
            let mut degree = 0;
            let mut sum = 0;
            for &child in &self.nodes[node_id].children {
                let leaves = self.nodes[child].leaves;
                if (leaves as f64) / (all_leaves as f64) >= rare_prob {
                    degree += 1;
                }
                sum += leaves;
            }
            if Some(node_id) != self.root {
                let leaves = all_leaves - sum;
                if (leaves as f64) / (all_leaves as f64) >= rare_prob {
                    degree += 1;
                }
            }
            self.nodes[node_id].frequent_degree = degree;
        }
        for node_id in 0..self.nodes.len() {
            if !self.nodes[node_id].is_leaf_type() {
                continue;
            }
            let mut parent = Some(node_id);
            while let Some(p) = parent {
                if self.nodes[p].frequent_degree != 1 {
                    break;
                }
                parent = self.nodes[p].parent;
            }
            if parent.is_some_and(|p| self.nodes[p].frequent_degree < 3) {
                self.nodes[node_id].frequent_degree = 0;
            }
        }
    }

    pub fn set_root(&mut self) {
        self.root = None;
        for node_id in 0..self.nodes.len() {
            if !self.nodes[node_id].detached && self.nodes[node_id].parent.is_none() {
                assert!(self.root.is_none());
                self.root = Some(node_id);
            }
        }
    }

    pub fn delete_transients(&mut self) -> usize {
        let mut n = 0;
        for node_id in 0..self.nodes.len() {
            if !self.nodes[node_id].detached && self.delete_transient_node(node_id) {
                n += 1;
            }
        }
        n
    }

    pub fn restrict_leaves(
        &mut self,
        leaf_names: &[String],
        delete_transient_ancestor: bool,
    ) -> usize {
        let mut n = 0;
        let leaf_names: BTreeSet<&str> = leaf_names.iter().map(String::as_str).collect();
        let node_vec: Vec<_> = self.nodes.iter().map(|node| node.id).collect();
        for node in node_vec {
            if self.nodes[node].detached || !self.nodes[node].is_leaf_type() {
                continue;
            }
            if !leaf_names.contains(self.nodes[node].name.as_str()) {
                self.delete_leaf(node, delete_transient_ancestor);
                n += 1;
            }
        }
        n
    }

    pub fn delete_leaf(&mut self, leaf: TreeNodeId, delete_transient_ancestor: bool) {
        let parent = self.nodes[leaf].parent;
        if let Some(parent) = parent {
            self.nodes[parent].children.retain(|&child| child != leaf);
        }
        self.nodes[leaf].detached = true;
        self.nodes[leaf].parent = None;
        if delete_transient_ancestor {
            if let Some(parent) = parent {
                self.delete_transient_node(parent);
            }
        }
    }

    pub fn sort(&mut self) {
        self.set_leaves();
        if let Some(root) = self.root {
            self.sort_node(root);
        }
    }

    fn sort_node(&mut self, node: TreeNodeId) {
        let children = self.nodes[node].children.clone();
        for child in children {
            self.sort_node(child);
        }
        let mut sorted = self.nodes[node].children.clone();
        sorted.sort_by(|&a, &b| {
            self.nodes[b]
                .leaves
                .cmp(&self.nodes[a].leaves)
                .then_with(|| {
                    let an = self.nodes[self.get_first_decendant(a)].name.as_str();
                    let bn = self.nodes[self.get_first_decendant(b)].name.as_str();
                    an.cmp(bn)
                })
        });
        self.nodes[node].children = sorted;
    }
}

impl TreeNode {
    pub fn qc(&self, tree: &Tree) -> Result<(), String> {
        if self.id >= tree.nodes.len() {
            return Err("tree node id is outside tree".to_string());
        }
        if self.name.contains(Tree::OBJ_NAME_SEPARATOR) {
            return Err(format!(
                "{:?} contains separator {:?}",
                self.name,
                Tree::OBJ_NAME_SEPARATOR
            ));
        }
        Ok(())
    }

    pub fn save_text(&self, tree: &Tree) -> String {
        let mut s = String::new();
        tree.save_node_text(self.id, &mut s, 0);
        s
    }

    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    pub fn is_leaf_type(&self) -> bool {
        self.is_leaf()
    }

    pub fn is_interior_type(&self) -> bool {
        !self.is_leaf()
    }

    pub fn get_parent_distance(&self) -> f64 {
        self.parent_distance
    }

    pub fn get_newick_name(&self, _minimal: bool) -> String {
        self.name.clone()
    }

    pub fn get_leaf_name(&self) -> String {
        self.name.clone()
    }

    pub fn get_human_name(&self, tree: &Tree) -> String {
        tree.get_lca_name(self.id)
    }

    pub fn get_save_subtree_p(&self) -> bool {
        true
    }

    pub fn name2newick(s: &str) -> String {
        url_escape(s)
    }
}

fn url_escape(s: &str) -> String {
    let mut out = String::new();
    for b in s.bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~') {
            out.push(b as char);
        } else {
            let _ = write!(out, "%{b:02X}");
        }
    }
    out
}

pub struct TypeNode<'a, T> {
    pub node: NodeId,
    pub t: &'a T,
}

pub fn topological_sort<T: Ord>(vec: &[T]) -> Vec<&T> {
    let mut gr = DiGraph::new();
    let mut nodes: Vec<TypeNode<'_, T>> = Vec::with_capacity(vec.len());
    for t in vec {
        let n = gr.add_node();
        for prev in &nodes {
            if prev.t < t {
                gr.add_arc(prev.node, n);
            } else if t < prev.t {
                gr.add_arc(n, prev.node);
            }
        }
        nodes.push(TypeNode { node: n, t });
    }

    fn process<'a, T: Ord>(
        gr: &mut DiGraph,
        nodes: &[TypeNode<'a, T>],
        node: NodeId,
        res: &mut Vec<&'a T>,
    ) {
        if gr.nodes[node].reachable {
            return;
        }
        gr.nodes[node].reachable = true;
        let in_arcs = gr.nodes[node].in_arcs.clone();
        for arc_id in in_arcs {
            if gr.arcs[arc_id].active {
                process(gr, nodes, gr.arcs[arc_id].from, res);
            }
        }
        res.push(nodes[node].t);
    }

    let mut res = Vec::with_capacity(vec.len());
    for node in 0..nodes.len() {
        process(&mut gr, &nodes, node, &mut res);
    }
    res
}
