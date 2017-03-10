package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.dijkstra;

import es.uma.lcc.neo.robustness.mo.shortestpath.utilities.MyPair;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.RoutingAlgorithm;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;

import java.util.*;


/**
 * Dijkstra algorithm
 */
public class DijkstraBounded implements RoutingAlgorithm {

    private GraphTable graph;
    // Suponemos que los arcos van numerados desde 1...N
    private float[] distances;
    //private Map<Long, Float> distances;
    private Map<Long, Long> predecessors;
    private TreeSet<MyPair> unvisited;

    private Long[] labs;
    private float[] weights;

    private float[] boundaries;

    public float[] getBoundaries() {
        return boundaries;
    }

    public void setBoundaries(float[] boundaries) {
        this.boundaries = boundaries;
    }

    DijkstraBounded() {
        //distances = new HashMap<Long, Float>();
        predecessors = new HashMap<Long, Long>();
        unvisited = new TreeSet<MyPair>();
    }

    public DijkstraBounded(Long[] labs) {
        this();
        this.labs = labs;
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
        unvisited.add(new MyPair(0f, from.getId()));

        while (!unvisited.isEmpty()) {
            MyPair top = unvisited.first();
            Long node = top.getV();
            unvisited.remove(top);

            if (node == to.getId()) {
                return computePath(to.getId());
            }
            // visit neightbors
            Map<Long, Long> neighbors = graph.getAdjacencyMatrix().row(node);
            float value;
            float edgeDistance;
            for (Map.Entry<Long, Long> neighbor : neighbors.entrySet()) {
                edgeDistance = getEdgeDistance(node, neighbor);
                value = distances[node.intValue()] + edgeDistance;
                if (distances[neighbor.getKey().intValue()] > value) {
                    distances[neighbor.getKey().intValue()] = value;
                    //distances.put(neighbor.getKey(), value);
                    predecessors.put(neighbor.getKey(), node);
                    unvisited.add(new MyPair(value, neighbor.getKey()));
                }
            }
        }

        return computePath(to.getId());
    }

    private float getEdgeDistance(Long source, Map.Entry<Long, Long> target) {
        float sum = 0f;
        for (int i = 0; i < labs.length; i++) {
            sum += (weights[i] * getGraph().getWeightsMatrix().get(target.getValue(), labs[i]));
        }
        return sum;
    }


    private List<Node> computePath(Long target) {
        List<Node> path = new ArrayList<Node>();
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

    public void setWeights(float[] weights) {
        this.weights = weights;
    }

}
