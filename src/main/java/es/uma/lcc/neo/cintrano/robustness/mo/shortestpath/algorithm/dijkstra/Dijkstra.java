package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;

import java.util.*;
import java.util.Map.Entry;

/**
 * Dijkstra algorithm
 */
public class Dijkstra {

    private GraphTable graph;
    private Map<Long, Float> distances;
    private Map<Long, Long> predecessors;
    private HashSet<Long> unvisited;
    private Node from;
    private Node to;

    Dijkstra() {
        distances = new HashMap<Long, Float>();
        predecessors = new HashMap<Long, Long>();
        unvisited = new HashSet<Long>();
    }

    GraphTable getGraph() {
        return graph;
    }

    public void setGraph(GraphTable graph) {
        this.graph = graph;
    }

    public List<Node> getPath(Long from, Long to) {
        //System.out.println("From:" + from.getId() + " To:" + to.getId());
        System.out.println("From:" + from + " To:" + to);
        if (from == to) { // path with a single node
            List<Node> path = new LinkedList<Node>();
            path.add(graph.getIntersections().get(from));
            return path;
        }

        setShortestDistance(from, 0f);
        unvisited.add(from);

//        System.out.println();
        while (!unvisited.isEmpty()) {
//            System.out.print(".");
            Long node = minDistance(unvisited);
            if (node == to) return computePath(to);
            unvisited.remove(node);
            visitNeighbors(node);
        }
//System.out.println(predecessors);
        System.out.println("___________ B: ");
        return computePath(to);
    }

    private void visitNeighbors(Long source) {
        // System.out.println("Visiting source " + source);
        Map<Long, Long> neighbors = getNeighbors(source);
        for (Entry<Long, Long> neighbor : neighbors.entrySet()) {
            // System.out.println("Neighbor: " + neighbor.getKey());
            if (getShortestDistance(neighbor.getKey()) > getShortestDistance(source) + getEdgeDistance(source, neighbor)) {
                setShortestDistance(neighbor.getKey(), getShortestDistance(source) + getEdgeDistance(source, neighbor));
                predecessors.put(neighbor.getKey(), source);
                unvisited.add(neighbor.getKey());
            }
        }
    }

    protected float getEdgeDistance(Long source, Entry<Long, Long> target) {
        return graph.getWeightsMatrix().get(target.getValue(), 0L);
    }

    private float getShortestDistance(long destId) {
        if (!distances.containsKey(destId)) {
            distances.put(destId, Float.MAX_VALUE);
        }
        return distances.get(destId);
    }

    private void setShortestDistance(long destId, float distance) {
        distances.put(destId, distance);
    }

    protected Map<Long, Long> getNeighbors(Long sourceId) {
        // System.out.println("Id" + sourceId + " cols:" +
        // graph.getAdjacencyMatrix().containsColumn(sourceId) + " rows:"
        // + graph.getAdjacencyMatrix().containsColumn(sourceId));
        return graph.getAdjacencyMatrix().row(sourceId);
    }

    private Long minDistance(HashSet<Long> nodes) {
        Long minimum = null;
        for (Long vertex : nodes) {
            if (minimum == null) {
                minimum = vertex;
            } else if (getShortestDistance(vertex) < getShortestDistance(minimum)) {
                minimum = vertex;
            }
        }
        return minimum;
    }

    private List<Node> computePath(Long target) {
        List<Node> path = new LinkedList<Node>();
        Long step = target;
        // check if a path exists
        if (predecessors.get(step) == null) {
            System.out.print(from + " ");
            System.out.println(to);
            if (getGraph().getAdjacencyMatrix().get(from.getId(), to.getId()) != null) {
            //if (getGraph().getAdjacencyMatrix().get(getGraph().getMapping().get(from.getId())
            //        , getGraph().getMapping().get(to.getId())) != null) {
                path.add(from);
                path.add(to);
                return path;
            }
//            System.out.println("No path found");
            return path;
        }
        path.add(graph.getIntersections().get(step));
        while (predecessors.get(step) != null) {
            step = predecessors.get(step);
            path.add(graph.getIntersections().get(step));
//            System.out.println("-------" + step + " "+ graph.getIntersections().get(step).getId());
        }
        // Put it into the correct order
        Collections.reverse(path);
        return path;
    }

}
