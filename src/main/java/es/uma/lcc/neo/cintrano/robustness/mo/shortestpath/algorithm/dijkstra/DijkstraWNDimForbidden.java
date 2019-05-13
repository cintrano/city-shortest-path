package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;

import java.util.List;
import java.util.Map.Entry;
import java.util.stream.Collectors;

public class DijkstraWNDimForbidden extends Dijkstra {

    private Long[] labs;
    private float[] weights;
    private List<Long> forbidden;

    public DijkstraWNDimForbidden(Long[] labs) {
        super();
        this.labs = labs;
    }

    public void setWeights(float[] weights) {
        this.weights = weights;
    }

    protected float getEdgeDistance(Long source, Entry<Long, Long> target) {
        if (forbidden.contains(target.getKey())) {
            return Float.MAX_VALUE;
        }
        float sum = 0f;
        for (int i = 0; i < labs.length; i++) {
            sum += (weights[i] * getGraph().getWeightsMatrix().get(target.getValue(), labs[i]));
        }
        return sum;
    }

    public void setForbidden(List<Node> fNodes) {
        //List<String> names = cars.stream().map( Car::getName ).collect( Collectors.toList() );
        //forbidden = forbiddenList;
        forbidden = fNodes.stream().map(Node::getId).collect(Collectors.toList());
    }
}
