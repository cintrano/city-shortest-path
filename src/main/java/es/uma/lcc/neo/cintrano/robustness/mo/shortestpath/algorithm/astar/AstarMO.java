package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.astar;


/**
 * Created by Christian Cintrano on 14/03/17.
 * A* algorithm version to 4-objectives
 */
public class AstarMO extends Astar{

    private float[] weights;

    public void setWeights(float[] weights) {
        this.weights = weights;
    }

    protected Float distBetween(Long current, Long neighbor) {
        return (getGraph().getWeightsMatrix().get(getGraph().getAdjacencyMatrix().get(current, neighbor), 0L) * weights[0]) +
                (getGraph().getWeightsMatrix().get(getGraph().getAdjacencyMatrix().get(current, neighbor), 1L) * weights[1]) +
                (getGraph().getWeightsMatrix().get(getGraph().getAdjacencyMatrix().get(current, neighbor), 2L) * weights[2]) +
                (getGraph().getWeightsMatrix().get(getGraph().getAdjacencyMatrix().get(current, neighbor), 3L) * weights[3]);
    }


}
