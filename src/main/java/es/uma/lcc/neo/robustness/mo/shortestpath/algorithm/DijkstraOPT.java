package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm;

import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;

import java.util.*;
import java.util.Map.Entry;

/**
 * Dijkstra algorithm
 */
public class DijkstraOPT implements RoutingAlgorithm {

    private GraphTable graph;
    // Suponemos que los arcos van numerados desde 1...N
    private float[] distances;
    //private Map<Long, Float> distances;
    private Map<Long, Long> predecessors;
    private HashSet<Long> unvisited;

    DijkstraOPT() {
        //distances = new HashMap<Long, Float>();
        predecessors = new HashMap<Long, Long>();
        unvisited = new HashSet<Long>();
    }

    GraphTable getGraph() {
        return graph;
    }

    public void setGraph(GraphTable graph) {
        this.graph = graph;
        distances = new float[graph.getIntersections().keySet().size() + 1]; // el 0 nunca se usa
        //for (Long node : graph.getIntersections().keySet()) {
        //    distances.put(node, Float.MAX_VALUE);
        //}
        for (int i = 0; i < distances.length; i++) {
            distances[i] = Float.MAX_VALUE;
        }
    }

    public List<Node> getPath(Node from, Node to) {
        //distances.put(from.getId(), 0f);
        distances[(int) from.getId()] = 0f;
        unvisited.add(from.getId());

        while (!unvisited.isEmpty()) {
            Long node = minDistance(unvisited);
            unvisited.remove(node);
            visitNeighbors(node);
        }

        return computePath(to.getId());
    }

    private void visitNeighbors(Long source) {
        Map<Long, Long> neighbors = getNeighbors(source);
        float value;
        for (Entry<Long, Long> neighbor : neighbors.entrySet()) {
            value = distances[source.intValue()] + getEdgeDistance(source, neighbor);
            if (distances[neighbor.getKey().intValue()] > value) {
                distances[neighbor.getKey().intValue()] = 0f;
                //distances.put(neighbor.getKey(), value);
                predecessors.put(neighbor.getKey(), source);
                unvisited.add(neighbor.getKey());
            }
        }
    }

    protected float getEdgeDistance(Long source, Entry<Long, Long> target) {
        return graph.getWeightsMatrix().get(target.getValue(), 0L);
    }

    private Map<Long, Long> getNeighbors(Long sourceId) {
        return graph.getAdjacencyMatrix().row(sourceId);
    }

    private Long minDistance(HashSet<Long> nodes) {
        Long minimum = null;
        for (Long vertex : nodes) {
            if (minimum == null || distances[vertex.intValue()] < distances[minimum.intValue()]) {
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
            System.out.println("No path found");
            return path;
        }
        path.add(graph.getIntersections().get(step));
        while (predecessors.get(step) != null) {
            step = predecessors.get(step);
            path.add(graph.getIntersections().get(step));
        }
        // Put it into the correct order
        Collections.reverse(path);
        return path;
    }

}
