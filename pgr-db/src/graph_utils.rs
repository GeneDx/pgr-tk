use core::cmp::Ord;
use core::fmt::Debug;
use petgraph::visit::{GraphRef, IntoNeighbors, VisitMap, Visitable};
use rustc_hash::FxHashMap;
use std::collections::BinaryHeap;
use std::hash::Hash;

// A struct to support node weight prioritized path traversal.
#[derive(Copy, Clone)]
pub struct WeightedNode<N>(u32, N);

impl<N> Ord for WeightedNode<N> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

impl<N> PartialOrd for WeightedNode<N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<N> PartialEq for WeightedNode<N> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<N> Eq for WeightedNode<N> {}

impl<N> Debug for WeightedNode<N>
where
    N: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {:?})", self.0, self.1)
    }
}

// To facilitate specical graph traversal for skew symmetrical graphs.
pub trait Node {
    fn reverse(&self) -> Self;
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
pub struct SNode (pub u64, pub u64, pub u8);

impl Node for SNode {
    fn reverse(&self) -> Self {
        SNode(self.0, self.1, 1-self.2)
    }
}
    


/// Code adapte from Petgraph's DFS
/// 
#[derive(Clone, Debug)]
pub struct WeightedDfs<'a, N, VM>
where
    N: Ord,
{
    /// The stack of nodes to visit
    pub priority_queue: BinaryHeap<WeightedNode<N>>,
    /// The map of discovered nodes
    pub discovered: VM,
    pub next_node: Option<WeightedNode<N>>,
    pub node_score: Option<&'a FxHashMap<N, u32>>,
}

impl<'a, N, VM> Default for WeightedDfs<'a, N, VM>
where
    VM: Default,
    N: Ord,
{
    fn default() -> Self {
        WeightedDfs {
            priority_queue: BinaryHeap::<WeightedNode<N>>::new(),
            discovered: VM::default(),
            next_node: None,
            node_score: None,
        }
    }
}

impl<'a, N, VM> WeightedDfs<'a, N, VM>
where
    N: Copy + PartialEq + Eq + Hash + Ord + Node,
    VM: VisitMap<N>,
{
    /// Create a new **Dfs**, using the graph's visitor map, and put **start**
    /// in the stack of nodes to visit.
    pub fn new<G>(graph: G, start: N, node_score: &'a FxHashMap<N, u32>) -> Self
    where
        G: GraphRef + Visitable<NodeId = N, Map = VM>,
    {
        let mut dfs = WeightedDfs::empty(graph, node_score);
        let s = node_score.get(&start).expect("Node not found");
        dfs.move_to(start);
        dfs.next_node = Some(WeightedNode(*s, start));
        dfs
    }

    /// Create a `Dfs` from a vector and a visit map
    pub fn from_parts(
        stack: BinaryHeap<WeightedNode<N>>,
        discovered: VM,
        node_score: &'a FxHashMap<N, u32>,
    ) -> Self {
        WeightedDfs {
            priority_queue: stack,
            discovered,
            next_node: None,
            node_score: Some(node_score),
        }
    }

    /// Clear the visit state
    pub fn reset<G>(&mut self, graph: G)
    where
        G: GraphRef + Visitable<NodeId = N, Map = VM>,
    {
        graph.reset_map(&mut self.discovered);
        self.priority_queue.clear();
    }

    /// Create a new **Dfs** using the graph's visitor map, and no stack.
    pub fn empty<G>(graph: G, node_score: &'a FxHashMap<N, u32>) -> Self
    where
        G: GraphRef + Visitable<NodeId = N, Map = VM>,
    {
        WeightedDfs {
            priority_queue: BinaryHeap::<WeightedNode<N>>::new(),
            next_node: None,
            discovered: graph.visit_map(),
            node_score: Some(node_score),
        }
    }

    /// Keep the discovered map, but clear the visit stack and restart
    /// the dfs from a particular node.
    pub fn move_to(&mut self, start: N) {
        let s = self.node_score.unwrap().get(&start).unwrap();
        let wn = WeightedNode(*s, start);
        self.priority_queue.clear();
        self.priority_queue.push(wn);
        self.next_node = Some(WeightedNode(*s, start));
    }

    /// Return the next node in the dfs, or **None** if the traversal is done.
    pub fn next<G>(&mut self, graph: G) -> Option<(N, bool)>
    where
        G: IntoNeighbors<NodeId = N>
    {
        loop {
            let node;
            if let Some(n) = self.next_node {
                node = n;
            } else {
                if self.priority_queue.is_empty() {
                    return None;
                }
                node = self.priority_queue.pop().unwrap();
            }

            if self.discovered.visit(node.1) {
                let rnode = node.1.reverse();
                self.discovered.visit(rnode);
                //println!("{:?} {:?}", node.1, rnode);

                let mut out_count = 0_usize;
                let mut succ_list = Vec::<WeightedNode<N>>::new();
                for succ in graph.neighbors(node.1) {
                    if !self.discovered.is_visited(&succ) {
                        out_count += 1;
                        let s = self.node_score.unwrap().get(&succ).unwrap();
                        succ_list.push(WeightedNode(*s, succ));
                    }
                }

                let mut is_leaf = true;
                if out_count == 0 {
                    is_leaf = false;
                    self.next_node = None;
                } else {
                    succ_list.sort();
                    self.next_node = succ_list.pop();
                    succ_list.iter().for_each(|s| self.priority_queue.push(*s));
                }
                return Some((node.1, is_leaf));
            }
        }
    }
}
