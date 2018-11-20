package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;

import java.util.*;
import java.util.Map.Entry;

/**
 * Dijkstra algorithm
 */
public class DijkstraLite {

    private GraphTable graph;
    private Float[] distances;
    private Long[] predecessors;
    private HashSet<Long> unvisited;
    private Long from;
    private Long to;

    DijkstraLite() {
        distances = new Float[graph.getIntersections().keySet().size()];
        for (int i = 0; i < distances.length; i++) {
            distances[i] = Float.MAX_VALUE;
        }
        predecessors = new Long[graph.getIntersections().keySet().size()];
        for (int i = 0; i < distances.length; i++) {
            distances[i] = null;
        }
        unvisited = new HashSet<>();
    }

    GraphTable getGraph() {
        return graph;
    }

    public void setGraph(GraphTable graph) {
        this.graph = graph;
    }

    public List<Node> getPath(Long from, Long to) {
        this.from = from;
        this.to = to;
        //System.out.println("From:" + from.getId() + " To:" + to.getId());

        distances[from.intValue()] = 0f;
        unvisited.add(from);

        while (!unvisited.isEmpty()) {
            Long node = minDistance(unvisited);
            if (Objects.equals(node, to)) return computePath(to);
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
            Float newValue = distances[source.intValue()] + getEdgeDistance(source, neighbor);
            int neighborId = neighbor.getKey().intValue();
            if (distances[neighborId] > newValue) {
                distances[neighborId] = newValue;
                predecessors[neighborId] = source;
                unvisited.add(neighbor.getKey());
            }
        }
    }

    protected float getEdgeDistance(Long source, Entry<Long, Long> target) {
        return graph.getWeightsMatrix().get(target.getValue(), 0L);
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
        if (predecessors[step.intValue()] == null) {
            if (getGraph().getAdjacencyMatrix().get(from, to) != null) {
                path.add(graph.getIntersections().get(from));
                path.add(graph.getIntersections().get(to));
                return path;
            }
            System.out.println("No path found");
            return path;
        }
        path.add(graph.getIntersections().get(step));
        while (predecessors[step.intValue()] != null) {
            step = predecessors[step.intValue()];
            path.add(graph.getIntersections().get(step));
        }
        // Put it into the correct order
        Collections.reverse(path);
        return path;
    }

}
