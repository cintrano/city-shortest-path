package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;

import java.util.*;
import java.util.Map.Entry;

/**
 * Dijkstra algorithm
 */
public class DijkstraSimple {

    private GraphTable graph;
    private Map<Long, Float> distances;
    private Map<Long, Long> predecessors;
    private Set<Long> unvisited;

    DijkstraSimple() {
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

    public List<Long> getPath(Long from, Long to) {
        //System.out.println("From:" + from.getId() + " To:" + to.getId());

        setShortestDistance(from, 0f);
        unvisited.add(from);

        while (!unvisited.isEmpty()) {
            Long node = minDistance(unvisited);
            if (node.equals(to)) return computePath(to);
            unvisited.remove(node);
            visitNeighbors(node);
        }
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

    private Long minDistance(Set<Long> nodes) {
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

    private List<Long> computePath(Long target) {
        List<Long> path = new ArrayList<Long>();
        Long step = target;
        // check if a path exists
        if (predecessors.get(step) == null) {
            System.out.println("No path found");
            return path;
        }
        path.add(step);
        while (predecessors.get(step) != null) {
            step = predecessors.get(step);
            path.add(step);
        }
        // Put it into the correct order
        Collections.reverse(path);
        return path;
    }

}
