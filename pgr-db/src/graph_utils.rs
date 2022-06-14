use core::cmp::Ord;
use petgraph::visit::{GraphRef, IntoNeighbors, IntoNeighborsDirected, VisitMap, Visitable};
use petgraph::EdgeDirection::{Incoming, Outgoing};
use rustc_hash::FxHashMap;
use std::collections::BinaryHeap;
use std::fmt::Debug;
use std::hash::Hash;

// A struct to support node weight prioritized path traversal.
#[derive(Copy, Clone)]
pub struct WeightedNode<N>(pub u32, pub N);

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
pub trait BiDiNode {
    fn reverse(&self) -> Self;
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
pub struct ShmmrGraphNode(pub u64, pub u64, pub u8);

impl BiDiNode for ShmmrGraphNode {
    fn reverse(&self) -> Self {
        ShmmrGraphNode(self.0, self.1, 1 - self.2)
    }
}

/// Code adapte from Petgraph's DFS
///
#[derive(Clone, Debug)]
pub struct BiDiGraphWeightedDfs<'a, N, VM>
where
    N: Ord + Debug,
{
    /// The stack of nodes to visit
    pub priority_queue: BinaryHeap<WeightedNode<N>>,
    /// The map of discovered nodes
    pub discovered: VM,
    pub next_node: Option<WeightedNode<N>>,
    current_branch: u32,
    branch_rank: u32,
    global_rank: FxHashMap<N, u32>,
    pub node_score: Option<&'a FxHashMap<N, u32>>,
}

impl<'a, N, VM> Default for BiDiGraphWeightedDfs<'a, N, VM>
where
    VM: Default,
    N: Ord + Debug,
{
    fn default() -> Self {
        BiDiGraphWeightedDfs {
            priority_queue: BinaryHeap::<WeightedNode<N>>::new(),
            discovered: VM::default(),
            next_node: None,
            current_branch: 0_u32,
            branch_rank: 0_u32,
            global_rank: FxHashMap::<N, u32>::default(),
            node_score: None,
        }
    }
}

impl<'a, N, VM> BiDiGraphWeightedDfs<'a, N, VM>
where
    N: Copy + PartialEq + Eq + Hash + Ord + BiDiNode + Debug,
    VM: VisitMap<N>,
{
    /// Create a new **Dfs**, using the graph's visitor map, and put **start**
    /// in the stack of nodes to visit.
    pub fn new<G>(graph: G, start: N, node_score: &'a FxHashMap<N, u32>) -> Self
    where
        G: GraphRef + Visitable<NodeId = N, Map = VM>,
    {
        let mut dfs = BiDiGraphWeightedDfs::empty(graph, node_score);
        let s = node_score.get(&start).expect("Node not found");
        dfs.move_to(start);
        dfs.next_node = Some(WeightedNode(*s, start));
        dfs.global_rank.insert(start, 0);
        dfs
    }

    /// Create a `Dfs` from a vector and a visit map
    pub fn from_parts(
        stack: BinaryHeap<WeightedNode<N>>,
        discovered: VM,
        node_score: &'a FxHashMap<N, u32>,
    ) -> Self {
        BiDiGraphWeightedDfs {
            priority_queue: stack,
            discovered,
            next_node: None,
            current_branch: 0_u32,
            branch_rank: 0_u32,
            global_rank: FxHashMap::<N, u32>::default(),
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
        self.global_rank.clear();
    }

    /// Create a new **Dfs** using the graph's visitor map, and no stack.
    pub fn empty<G>(graph: G, node_score: &'a FxHashMap<N, u32>) -> Self
    where
        G: GraphRef + Visitable<NodeId = N, Map = VM>,
    {
        BiDiGraphWeightedDfs {
            priority_queue: BinaryHeap::<WeightedNode<N>>::new(),
            next_node: None,
            current_branch: 0_u32,
            branch_rank: 0_u32,
            global_rank: FxHashMap::<N, u32>::default(),
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
        self.global_rank.insert(start, 0);
    }

    /// Return the next node in the dfs, or **None** if the traversal is done.
    pub fn next<G>(&mut self, graph: G) -> Option<(N, Option<N>, bool, u32, u32, u32)>
    where
        G: IntoNeighbors<NodeId = N> + IntoNeighborsDirected<NodeId = N>,
    {
        let mut branch_rank;
        let global_rank = &mut self.global_rank;
        let mut branch = self.current_branch;
        loop {
            let node;
            if let Some(n) = self.next_node { // the next_node is the prioritized node 
                node = n;
                branch_rank = self.branch_rank;
            } else {
                if self.priority_queue.is_empty() {
                    return None;
                }
                node = self.priority_queue.pop().unwrap();
                self.branch_rank = 0;
                branch_rank = 0;
                self.current_branch += 1;
                branch = self.current_branch;
            }
            //println!("DBG: current node: {:?}", node);

            if self.discovered.visit(node.1) {
                //the node is not visited before
                let rnode = node.1.reverse();
                self.discovered.visit(rnode);
                //println!("DBG, visited: {:?}, {:?}", node, rnode);

                let mut out_count = 0_usize;
                let mut succ_list_f = Vec::<WeightedNode<N>>::new();
                for succ in graph.neighbors_directed(node.1, Outgoing) {
                    //println!("DBG: succ: {:?} {:?}", node.1, succ);
                    if !self.discovered.is_visited(&succ) {
                        //println!("DBG: pushing0: {:?}", succ);
                        out_count += 1;
                        let s = self.node_score.unwrap().get(&succ).unwrap();
                        succ_list_f.push(WeightedNode(*s, succ));
                    }
                }

                let mut succ_list_r = Vec::<WeightedNode<N>>::new();
                for succ in graph.neighbors_directed(node.1.reverse(), Outgoing) {
                    //println!("DBG: succ: {:?} {:?}", node.1, succ);
                    if !self.discovered.is_visited(&succ) {
                        //println!("DBG: pushing0: {:?}", succ);
                        let s = self.node_score.unwrap().get(&succ).unwrap();
                        succ_list_r.push(WeightedNode(*s, succ));
                    }
                }

                let mut is_leaf = false;
                if out_count == 0 {
                    is_leaf = true;
                    self.next_node = None;
                } else {
                    if succ_list_f.len() > 0 {
                        // we prefer the same direction first
                        succ_list_f.sort();
                        self.next_node = succ_list_f.pop();
                        succ_list_f.iter().for_each(|s| {
                            //println!("DBG, pushing1: {:?}", s);
                            self.priority_queue.push(*s);
                        });
                    } 

                    succ_list_r.sort();
                    succ_list_r.iter().for_each(|s| {
                        //println!("DBG, pushing1: {:?}", s);
                        self.priority_queue.push(*s);
                    });
                }
                //println!("DBG: next node: {:?}", self.next_node);

                let mut node_rank = u32::MAX;
                let mut p_node: Option<N> = None;
                graph.neighbors_directed(node.1, Incoming).for_each(|n| {
                    if let Some(r) = global_rank.get(&n) {
                        if *r < node_rank {
                            node_rank = *r;
                            p_node = Some(n)
                        }
                    }
                });

                graph
                    .neighbors_directed(node.1.reverse(), Incoming)
                    .for_each(|n| {
                        if let Some(r) = global_rank.get(&n) {
                            if *r < node_rank {
                                node_rank = *r;
                                p_node = Some(n);
                            }
                        }
                    });
                
                if node_rank == u32::MAX {
                    node_rank = 0;
                }
                node_rank += 1;
                global_rank.insert(node.1, node_rank);
                global_rank.insert(node.1.reverse(), node_rank);

                self.branch_rank += 1;
                //println!("DBG: out {:?}", node.1);
                return Some((node.1, p_node, is_leaf, node_rank, branch, branch_rank));
            } // else continue the loop
        }
    }
}
