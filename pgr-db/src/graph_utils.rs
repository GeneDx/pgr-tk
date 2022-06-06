use petgraph::visit::{GraphRef, IntoNeighbors, VisitMap, Visitable};
use rustc_hash::FxHashMap;
use std::hash::Hash;

/// Visit nodes of a graph in a depth-first-search (DFS) emitting nodes in
/// preorder (when they are first discovered).
///
/// The traversal starts at a given node and only traverses nodes reachable
/// from it.
///
/// `Dfs` is not recursive.
///
/// `Dfs` does not itself borrow the graph, and because of this you can run
/// a traversal over a graph while still retaining mutable access to it, if you
/// use it like the following example:
///
/// ```
/// use petgraph::Graph;
/// use petgraph::visit::Dfs;
///
/// let mut graph = Graph::<_,()>::new();
/// let a = graph.add_node(0);
///
/// let mut dfs = Dfs::new(&graph, a);
/// while let Some(nx) = dfs.next(&graph) {
///     // we can access `graph` mutably here still
///     graph[nx] += 1;
/// }
///
/// assert_eq!(graph[a], 1);
/// ```
///
/// **Note:** The algorithm may not behave correctly if nodes are removed
/// during iteration. It may not necessarily visit added nodes or edges.
#[derive(Clone, Debug)]
pub struct WeightedDfs<'a, N, VM> {
    /// The stack of nodes to visit
    pub stack: Vec<N>,
    /// The map of discovered nodes
    pub discovered: VM,
    pub node_score: Option<&'a FxHashMap<N, u32>>,
}

impl<'a, N, VM> Default for WeightedDfs<'a, N, VM>
where
    VM: Default,
{
    fn default() -> Self {
        WeightedDfs {
            stack: Vec::new(),
            discovered: VM::default(),
            node_score: None,
        }
    }
}

impl<'a, N, VM> WeightedDfs<'a, N, VM>
where
    N: Copy + PartialEq + Eq + Hash,
    VM: VisitMap<N>,
{
    /// Create a new **Dfs**, using the graph's visitor map, and put **start**
    /// in the stack of nodes to visit.
    pub fn new<G>(graph: G, start: N, node_score: &'a FxHashMap<N, u32>) -> Self
    where
        G: GraphRef + Visitable<NodeId = N, Map = VM>,
    {
        let mut dfs = WeightedDfs::empty(graph, node_score);
        dfs.move_to(start);
        dfs
    }

    /// Create a `Dfs` from a vector and a visit map
    pub fn from_parts(stack: Vec<N>, discovered: VM, node_score: &'a FxHashMap<N, u32>) -> Self {
        WeightedDfs {
            stack,
            discovered,
            node_score: Some(node_score),
        }
    }

    /// Clear the visit state
    pub fn reset<G>(&mut self, graph: G)
    where
        G: GraphRef + Visitable<NodeId = N, Map = VM>,
    {
        graph.reset_map(&mut self.discovered);
        self.stack.clear();
    }

    /// Create a new **Dfs** using the graph's visitor map, and no stack.
    pub fn empty<G>(graph: G, node_score: &'a FxHashMap<N, u32>) -> Self
    where
        G: GraphRef + Visitable<NodeId = N, Map = VM>,
    {
        WeightedDfs {
            stack: Vec::new(),
            discovered: graph.visit_map(),
            node_score: Some(node_score),
        }
    }

    /// Keep the discovered map, but clear the visit stack and restart
    /// the dfs from a particular node.
    pub fn move_to(&mut self, start: N) {
        self.stack.clear();
        self.stack.push(start);
    }

    /// Return the next node in the dfs, or **None** if the traversal is done.
    pub fn next<G>(&mut self, graph: G) -> Option<(N, bool)>
    where
        G: IntoNeighbors<NodeId = N>,
    {
        while let Some(node) = self.stack.pop() {
            if self.discovered.visit(node) {
                let mut out_count = 0_usize;
                let mut succs = Vec::<(u32, N)>::new();
                for succ in graph.neighbors(node) {
                    if !self.discovered.is_visited(&succ) {
                        out_count += 1;
                        let s = self.node_score.unwrap().get(&succ).unwrap();
                        succs.push((*s, succ));
                    }
                    succs.sort_by(|a, b| a.0.cmp(&b.0));
                    succs.iter().for_each(|(_s, succ)| {
                        self.stack.push(*succ);
                    });
                }
                let mut is_leaf = true;
                if out_count == 0 {
                    is_leaf = false;
                    // dead end
                    let mut best_score_node: Option<N> = None;
                    let mut best_score = std::u32::MIN;
                    self.stack.iter().for_each(|node| {
                        if let Some(score) = self.node_score.unwrap().get(&node) {
                            if *score >= best_score {
                                best_score = *score;
                                best_score_node = Some(*node);
                            }
                        }
                    });
                    self.stack.clear();
                    if let Some(node) = best_score_node {
                        self.stack.push(node);
                    }
                }
                return Some((node, is_leaf));
            }
        }
        None
    }
}
