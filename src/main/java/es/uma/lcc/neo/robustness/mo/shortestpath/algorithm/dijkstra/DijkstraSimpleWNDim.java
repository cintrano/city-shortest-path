package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.dijkstra;

import java.util.Map.Entry;

public class DijkstraSimpleWNDim extends DijkstraSimple {

    private Long[] labs;
    private float[] weights;

    public DijkstraSimpleWNDim(Long[] labs) {
        super();
        this.labs = labs;
    }

    public void setWeights(float[] weights) {
        this.weights = weights;
    }

    protected float getEdgeDistance(Long source, Entry<Long, Long> target) {
        float sum = 0f;
        for (int i = 0; i < labs.length; i++) {
            sum += (weights[i] * getGraph().getWeightsMatrix().get(target.getValue(), labs[i]));
        }
        return sum;
    }

}
