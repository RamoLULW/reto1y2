from dataclasses import dataclass
from collections import defaultdict, deque
import math
import heapq
import sys
import os

BASE_PATH = os.path.dirname(os.path.abspath(__file__))

@dataclass
class Node:
    nid: int
    x: float
    y: float
    is_source: bool = False

class Graph:
    def __init__(self):
        self.nodes: dict[int, Node] = {}
        self.sources: set[int] = set()
        self.office: int | None = None

        self.edges: list[tuple[int, int, float, float]] = []

        self.adj: dict[int, list[tuple[int, float, float]]] = defaultdict(list)

        self.new_specs: list[tuple[float, float, float]] = []

    def add_node(self, nid: int, x: float, y: float, is_source: bool):
        self.nodes[nid] = Node(nid, x, y, is_source)
        if is_source:
            self.sources.add(nid)

    def add_edge(self, u: int, v: int, capacity: float):
        nu, nv = self.nodes[u], self.nodes[v]
        length = math.hypot(nu.x - nv.x, nu.y - nv.y)
        self.edges.append((u, v, capacity, length))

        self.adj[u].append((v, length, capacity))
        self.adj[v].append((u, length, capacity))

def read_graph(path: str) -> Graph:
    g = Graph()
    with open(path, "r") as f:
        _ = f.readline()

        for line in f:
            if line.strip() == "[NODES]":
                break

        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == "[EDGES]":
                break
            nid, x, y, flag = line.split()
            g.add_node(int(nid), float(x), float(y), int(flag) == 1)

        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == "[OFFICE]":
                office_line = f.readline().strip()
                g.office = int(office_line)
                break
            u, v, cap = line.split()
            g.add_edge(int(u), int(v), float(cap))

        for line in f:
            if line.strip() == "[NEW]":
                 break
        else:
            return g
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            x, y, cap = parts
            g.new_specs.append((float(x), float(y), float(cap)))

    return g

def dijkstra_multi(g: Graph):
    dist = {nid: math.inf for nid in g.nodes}
    source_of: dict[int, int] = {}
    pq: list[tuple[float, int, int]] =[]

    for s in g.sources:
        dist[s] = 0.0
        source_of[s] = s
        heapq.heappush(pq, (0.0, s, s))

    while pq:
        d, u, src = heapq.heappop(pq)
        if d != dist[u]:
            continue
        for v, length, _cap in g.adj[u]:
            nd = d + length
            if nd < dist[v]:
                dist[v] = nd
                source_of[v] = src
                heapq.heappush(pq, (nd, v, src))

    return dist, source_of

def sectorization(g: Graph):
    dist, source_of = dijkstra_multi(g)

    sectors: dict[int, list[int]] = defaultdict(list)
    for v, s in source_of.items():
        sectors[s].append(v)

    closed_edges: set[tuple[int, int]] = set()
    open_edges: list[tuple[int, int, float, float]] = []
    for u, v, cap, length in g.edges:
        if source_of[u] != source_of[v]:
            closed_edges.add((min(u, v), max(u, v)))
        else:
            open_edges.append((u, v, cap, length))

    return sectors, closed_edges, open_edges, dist

def farthest_nodes_per_sector(sectors, dist):
    farthest: dict[int, tuple[int, float]] = {}
    for src, nodes in sectors.items():
        best_node = None
        best_dist = -1.0
        for v in nodes:
            if dist[v] > best_dist:
                best_dist = dist[v]
                best_node = v
        farthest[src] = (best_node, best_dist)
    return farthest

def build_capacity_graph_4sector(g: Graph, sector_nodes: list[int], closed_edges: set[tuple[int, int]],):
    allowed = set(sector_nodes)
    banned = {tuple(sorted(e)) for e in closed_edges}
    cap: dict[int, dict[int, float]] = defaultdict(lambda: defaultdict(float))

    for u, v, c, _length in g.edges:
        if u in allowed and v in allowed and tuple(sorted((u, v))) not in banned: 
            cap[u][v] += c
            cap[v][u] += c

    return cap

def edmons_karp(cap: dict[int, dict[int, float]], s: int, t: int) -> float:
    res: dict[int, dict[int, float]] = defaultdict(dict)
    for u in cap:
        for v, c in cap[u].items():
            res[u][v] = c

    max_flow = 0.0
    eps = 1e-9

    while True:
        parent: dict[int, int | None] = {s: None}
        q = deque([s])

        while q and t not in parent:
            u = q.popleft()
            for v, c in res.get(u, {}).items():
                if v not in parent and c > eps:
                    parent[v] = u
                    q.append(v)

        if t not in parent:
            break

        bottleneck = float("inf")
        v = t
        while parent[v] is not None:
            u = parent[v]
            bottleneck = min(bottleneck, res[u][v])
            v = u

        v = t
        while parent[v] is not None:
            u = parent[v]
            res[u][v] -= bottleneck
            res.setdefault(v, {}).setdefault(u, 0.0)
            res[v][u] += bottleneck
            v = u

        max_flow += bottleneck
    
    return max_flow

def max_flow_per_sector(g: Graph, sectors, closed_edges, farthest):
    results = {}
    for src, nodes in sectors.items():
        t, _ = farthest[src]
        cap = build_capacity_graph_4sector(g, nodes, closed_edges)
        flow = edmons_karp(cap, src, t)
        results[src] = (t, flow)
    return results

def prim_mst(g: Graph):
    start = g.office
    if start is None: 
        raise ValueError("Office node not defined")
    
    visited = {start}
    mst_adj: dict[int, list[tuple[int, float]]] = defaultdict(list)
    heap: list[tuple[float, int, int]] = []

    for v, length, _cap in g.adj[start]:
        heapq.heappush(heap, (length, start, v))

    while heap and len(visited) < len(g.nodes):
        length, u, v = heapq.heappop(heap)
        if v in visited:
            continue
        visited.add(v)
        mst_adj[u].append((v, length))
        mst_adj[v].append((u, length))
        for w, l, _c in g.adj[v]:
            if w not in visited:
                heapq.heappush(heap, (l, v, w))

    return mst_adj

def tsp_tour_from_mst(g: Graph, mst_adj):
    start = g.office
    if start is None:
        raise ValueError("Office node not defined")
    
    walk = []
    visited = set()

    def dfs(u: int):
        visited.add(u)
        walk.append(u)
        for v, _ in mst_adj[u]:
            if v not in visited:
                dfs(v)
                walk.append(u)

    dfs(start)

    seen = set()
    tour: list[int] = []
    for v in walk:
        if v not in seen:
            seen.add(v)
            tour.append(v)
    tour.append(start)

    return tour

def tour_length(g: Graph, tour: list[int]) -> float:
    total = 0.0
    for a, b in zip(tour, tour[1:]):
        na, nb = g.nodes[a], g.nodes[b]
        total += math.hypot(na.x - nb.x, na.y - nb.y)
    return total

def add_new_node(g: Graph, spec: tuple[float, float, float]):
    x, y, capacity = spec
    new_id = max(g.nodes) + 1

    best_node = None
    best_dist = float("inf")
    for nid, node in g.nodes.items():
        if node.is_source:
            continue
        d = math.hypot(node.x - x, node.y - y)
        if d < best_dist:
            best_dist = d
            best_node = nid

    if best_node is None:
        raise ValueError("No non-source node found to connect the new node")
    
    g.add_node(new_id, x, y, is_source=False)
    g.add_edge(new_id, best_node, capacity)
    return new_id, best_node, best_dist

def check_file(path: str):
    g = read_graph(path)
    print("1. Instances")
    print(f"Nodes: {len(g.nodes)} and Edges: {len(g.edges)}")
    print(f"Sources: {sorted(g.sources)}")
    print(f"Office: {g.office}\n")

    print("2. Pipe lengths")
    print("Computed length for all pipes\n")

    sectors, closed_edges, open_edges, dist = sectorization(g)
    print("3. Sectorization")
    for src, nodes in sectors.items():
        print(f"Source {src}: {len(nodes)} nodes")
    print(f"Closed pipes: {len(closed_edges)}\n")

    print("4. Water Freshness (slowest node)")
    farthest = farthest_nodes_per_sector(sectors, dist)
    for src, (node, d) in farthest.items():
        print(f"Source {src} \u2192 Node {node} (distance = {d:.2f})")
    print()

    print("5. Maximum Flow per Sector")
    flows = max_flow_per_sector(g, sectors, closed_edges, farthest)
    for src, (t, flow) in flows.items():
        print(f"Source {src} \u2192 Node {t}: max flow = {flow:.2f}")
    print()

    print("6. Water quality sampling (TSP approx)")
    mst_adj = prim_mst(g) 
    tour = tsp_tour_from_mst(g, mst_adj)
    L = tour_length(g, tour)
    print(f"Tour length \u2248 {L:.2f}")
    print(f"Nodes visited: {len(tour) - 1}\n")

    print("7. Network expansion")
    if g.new_specs:
        new_id, attach_to, d = add_new_node(g, g.new_specs[0])
        print(f"Added new node {new_id}")
        print(f"Connected to node {attach_to} (distance = {d:.2f})\n")
    else:
        print("No new nodes provided\n")



if __name__ == "__main__":
    if len(sys.argv) >= 2:
        filepath = sys.argv[1]
    else:
        filepath = os.path.join(BASE_PATH, "HAN.txt")
    check_file(filepath)