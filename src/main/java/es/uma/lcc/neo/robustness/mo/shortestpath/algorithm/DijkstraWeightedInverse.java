package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm;

import java.util.Map;
import java.util.Map.Entry;

public class DijkstraWeightedInverse extends Dijkstra {

    private Long lab1;
    private Long lab2;
    private float weighing;

    public DijkstraWeightedInverse(Long lab1, Long lab2, float weighing) {
        super();
        this.lab1 = lab1;
        this.lab2 = lab2;
        this.weighing = weighing;
    }

    protected float getEdgeDistance(Long source, Entry<Long, Long> target) {
        Float w1 = getGraph().getWeightsMatrix().get(target.getValue(), lab1);
        Float w2 = getGraph().getWeightsMatrix().get(target.getValue(), lab2);
        return (w1 * weighing) + (w2 * (1- weighing));
    }

    protected Map<Long, Long> getNeighbors(Long sourceId) {
        // System.out.println("Id" + sourceId + " cols:" +
        // graph.getAdjacencyMatrix().containsColumn(sourceId) + " rows:"
        // + graph.getAdjacencyMatrix().containsColumn(sourceId));
        return getGraph().getAdjacencyMatrix().column(sourceId);
    }
}
