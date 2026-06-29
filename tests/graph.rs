use amrfinder::graph::{topological_sort, DiGraph, LcaBuffer, Tree};
use std::collections::BTreeSet;

#[test]
fn scc_simple() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    let c = g.add_node();

    g.add_arc(a, b);
    g.add_arc(b, c);
    g.add_arc(c, a);

    g.scc();

    assert_eq!(g.nodes[a].scc, g.nodes[b].scc);
    assert_eq!(g.nodes[b].scc, g.nodes[c].scc);
}

#[test]
fn scc_two_components() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    let c = g.add_node();
    let d = g.add_node();

    g.add_arc(a, b);
    g.add_arc(b, a);
    g.add_arc(c, d);

    g.scc();

    assert_eq!(g.nodes[a].scc, g.nodes[b].scc);
    assert_ne!(g.nodes[c].scc, g.nodes[d].scc);
}

#[test]
fn set_reachable_follows_original_arcs_false_direction() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    let c = g.add_node();
    let d = g.add_node();

    g.add_arc(a, b);
    g.add_arc(b, c);

    g.set_reachable(c, true);

    assert!(g.nodes[a].reachable);
    assert!(g.nodes[b].reachable);
    assert!(g.nodes[c].reachable);
    assert!(!g.nodes[d].reachable);
}

#[test]
fn roots_leaves_and_single_root() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    let c = g.add_node();

    g.add_arc(a, b);
    g.add_arc(b, c);

    assert_eq!(g.get_ends(false), vec![a]);
    assert_eq!(g.get_ends(true), vec![c]);
    assert_eq!(g.roots(), vec![a]);
    assert_eq!(g.leaves(), vec![c]);
    assert_eq!(g.get_root(false), Some(a));
}

#[test]
fn incident_neighborhood_degree_and_self_loop() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    let c = g.add_node();

    let ab = g.add_arc(a, b);
    g.add_arc(c, b);
    let cc = g.add_arc(c, c);

    assert!(g.is_incident(a, b, true));
    assert!(g.is_incident(b, a, false));
    assert!(g.is_incident_except(b, a, false));
    assert_eq!(g.get_degree(b), 2);
    assert_eq!(g.get_neighborhood(a, true), vec![b]);
    assert_eq!(g.get_neighborhood(b, false), vec![a, c]);
    assert_eq!(g.get_children(b), vec![a, c]);
    assert!(!g.self_loop(ab));
    assert!(g.self_loop(cc));
}

#[test]
fn set_node_updates_arc_endpoints_and_adjacency_lists() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    let c = g.add_node();
    let d = g.add_node();

    let arc = g.add_arc(a, b);
    g.set_node(arc, c, true);
    assert_eq!(g.arcs[arc].to, c);
    assert!(!g.nodes[b].in_arcs.contains(&arc));
    assert!(g.nodes[c].in_arcs.contains(&arc));

    g.set_node(arc, d, false);
    assert_eq!(g.arcs[arc].from, d);
    assert!(!g.nodes[a].out_arcs.contains(&arc));
    assert!(g.nodes[d].out_arcs.contains(&arc));
}

#[test]
fn connected_components_sets_disjoint_cluster_parent_state() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    let c = g.add_node();

    g.add_arc(a, b);
    g.connected_components();

    assert_eq!(g.nodes[b].parent_dc, a);
    assert_eq!(g.nodes[c].parent_dc, c);
}

#[test]
fn contract_scc_merges_cycle_and_removes_duplicate_arc() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    let c = g.add_node();

    g.add_arc(a, b);
    g.add_arc(b, a);
    g.add_arc(a, c);
    g.add_arc(b, c);

    g.scc();
    let root = g.nodes[a].scc.unwrap();
    g.contract_scc();

    assert_eq!(g.get_neighborhood(root, true), vec![c]);
    assert_eq!(g.get_degree(if root == a { b } else { a }), 0);
}

#[test]
fn copy_reverse_and_borrow_arcs_preserve_original_mapping_behavior() {
    let mut other = DiGraph::new();
    let a = other.add_node();
    let b = other.add_node();
    let c = other.add_node();
    other.add_arc(a, b);
    other.add_arc(b, c);

    let copied = other.copy();
    assert_eq!(copied.nodes.len(), 3);
    assert_eq!(copied.get_neighborhood(0, true), vec![1]);
    assert_eq!(copied.get_neighborhood(1, true), vec![2]);

    let mut mapping = std::collections::HashMap::new();
    mapping.insert(a, 10);
    mapping.insert(b, 11);
    let reversed = DiGraph::reverse(&mapping);
    assert_eq!(reversed[&10], a);
    assert_eq!(reversed[&11], b);

    let mut target = DiGraph::new();
    let ta = target.add_node();
    let tb = target.add_node();
    let tc = target.add_node();
    target.add_arc(ta, tb);

    let mut other2this = std::collections::HashMap::new();
    other2this.insert(a, ta);
    other2this.insert(b, tb);
    other2this.insert(c, tc);
    target.borrow_arcs(&other, &other2this, false);

    assert_eq!(target.get_neighborhood(ta, true), vec![tb]);
    assert_eq!(target.get_neighborhood(tb, true), vec![tc]);
}

#[test]
fn graph_qc_and_save_text_cover_original_hooks() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    g.add_arc(a, b);
    g.scc();

    assert!(g.qc().is_ok());
    let text = g.save_text();
    assert!(text.contains("Out:"));
    assert!(text.contains("In:"));
    assert!(g.nodes[a].qc(&g).is_ok());
    assert!(g.nodes[a].save_text(&g).contains("DFS_order"));
}

#[test]
fn disjoint_cluster_api_merges_and_finds_roots() {
    let mut g = DiGraph::new();
    let a = g.add_node();
    let b = g.add_node();
    let c = g.add_node();

    g.disjoint_cluster_init_all();
    g.disjoint_cluster_merge(a, b);

    assert_eq!(g.disjoint_cluster_root_const(b), a);
    assert_eq!(g.disjoint_cluster_root(c), c);
    assert_eq!(g.nodes[b].get_disjoint_cluster(), a);
}

#[test]
fn tree_lca_path_and_parent_sets_match_upstream_shape() {
    let mut tree = Tree::new();
    let root = tree.add_node("root", None);
    let left = tree.add_node("left", Some(root));
    tree.nodes[left].parent_distance = 1.0;
    let right = tree.add_node("right", Some(root));
    tree.nodes[right].parent_distance = 2.0;
    let ll = tree.add_node("aa", Some(left));
    tree.nodes[ll].parent_distance = 3.0;
    let lr = tree.add_node("bb", Some(left));
    tree.nodes[lr].parent_distance = 4.0;

    assert!(tree.qc().is_ok());
    assert_eq!(tree.root, Some(root));
    assert_eq!(tree.get_parent(left), Some(root));
    assert_eq!(tree.get_children(left), vec![ll, lr]);
    assert_eq!(tree.get_ancestor(ll, 2), root);
    assert_eq!(tree.get_topological_depth(ll), 2);

    let mut buf = LcaBuffer::default();
    assert_eq!(tree.get_lca(Some(ll), Some(lr), &mut buf), Some(left));
    assert_eq!(tree.get_lca_vec(&[ll, lr, right], &mut buf), Some(root));

    let parents = tree.get_parents(&[ll, lr], &mut buf);
    assert_eq!(parents, BTreeSet::from([ll, lr]));

    let mut lca = None;
    let path = tree.get_path(ll, right, Some(root), &mut lca, &mut buf);
    assert_eq!(lca, Some(root));
    assert_eq!(path, vec![ll, left, right]);
}

#[test]
fn tree_distances_leaves_and_text_outputs() {
    let mut tree = Tree::new();
    let root = tree.add_node("root", None);
    let a = tree.add_node("a leaf", Some(root));
    tree.nodes[a].parent_distance = 1.5;
    let b = tree.add_node("b", Some(root));
    tree.nodes[b].parent_distance = 2.5;

    tree.set_leaves();
    assert_eq!(tree.nodes[root].leaves, 2);
    assert_eq!(tree.get_length(), 4.0);
    assert_eq!(tree.get_ave_arc_length(), 2.0);
    assert_eq!(tree.get_height(root), 1);
    assert_eq!(tree.get_distance_height(root), 2.5);
    assert_eq!(tree.get_lca_name(root), "a leaf:b");
    assert_eq!(
        tree.print_newick(false, false),
        "(a%20leaf:1.500000,b:2.500000);\n"
    );
    assert!(tree.save_text().contains("root:"));

    let distances = tree.get_leaf_distances();
    assert_eq!(distances.len(), 1);
    assert_eq!(distances[0].leaf1, a);
    assert_eq!(distances[0].leaf2, b);
    assert_eq!(distances[0].distance, 4.0);
}

#[test]
fn tree_area_closest_leaves_and_restrict_leaves() {
    let mut tree = Tree::new();
    let root = tree.add_node("root", None);
    let a = tree.add_node("a", Some(root));
    tree.nodes[a].parent_distance = 1.0;
    let b = tree.add_node("b", Some(root));
    tree.nodes[b].parent_distance = 2.0;
    let c = tree.add_node("c", Some(b));
    tree.nodes[c].parent_distance = 3.0;

    let mut area = Vec::new();
    let mut boundary = Vec::new();
    tree.get_area(root, 1, &mut area, &mut boundary);
    assert_eq!(area, vec![root, a, b]);
    assert_eq!(boundary, vec![a, b]);

    let mut closest = Vec::new();
    tree.get_closest_leaves(root, 2, &mut closest);
    assert_eq!(closest[0].node, a);
    assert_eq!(closest[1].node, c);

    let keep = vec!["c".to_string()];
    assert_eq!(tree.restrict_leaves(&keep, false), 1);
    assert!(tree.nodes[a].detached);
}

#[test]
fn leaves2lcas_and_frequent_metrics_are_available() {
    let mut tree = Tree::new();
    let root = tree.add_node("root", None);
    let left = tree.add_node("left", Some(root));
    let right = tree.add_node("right", Some(root));
    let a = tree.add_node("a", Some(left));
    let b = tree.add_node("b", Some(left));
    let c = tree.add_node("c", Some(right));

    let lcas = tree.leaves2lcas(&BTreeSet::from([a, b, c]));
    assert_eq!(lcas, vec![root]);

    tree.set_frequent_child(0.1);
    assert!(tree.nodes[root].frequent_child);
    tree.set_frequent_degree(0.1);
    assert!(tree.nodes[root].frequent_degree >= 2);
    assert_eq!(Tree::radius2boundary_size(0), 1);
    assert_eq!(Tree::radius2boundary_size(2), 6);
}

#[test]
fn topological_sort_returns_ascending_order_using_original_algorithm() {
    let values = vec![3, 1, 2, 1];
    let sorted: Vec<_> = topological_sort(&values).into_iter().copied().collect();
    assert_eq!(sorted, vec![1, 1, 2, 3]);
}
