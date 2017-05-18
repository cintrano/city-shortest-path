package es.uma.lcc.neo.robustness.mo.shortestpath.eaparameters.utilities;

import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDim;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;

import java.util.List;

/**
 * Created by Christian Cintrano on 18/05/17.
 * Dummy class to the individual of the evolutive algorithm
 */
public class EASolution {
    // Variables of the solution
    private int iteration;
    private float[] objectives;

    // Variables auxiliares
    private GraphTable graph;
    private Long start, end; // start -> end

    // Variables of the path
    private float[] fitness;
    private List<Node> variables;

    public EASolution(int iteration, GraphTable graph, Long start, Long end) {
        this.iteration = iteration;
        this.graph = graph;
        this.start = start;
        this.end = end;
    }

    public int getIteration() {
        return iteration;
    }

    public void setIteration(int iteration) {
        this.iteration = iteration;
    }

    public float[] getObjectives() {
        return objectives;
    }

    public void setObjectives(float[] objectives) {
        this.objectives = objectives;
    }

    public List<Node> getVariables() {
        return variables;
    }

    public void setVariables(List<Node> variables) {
        this.variables = variables;
    }

    public float[] getFitness() {
        return fitness;
    }

    public void setFitness(float[] fitness) {
        this.fitness = fitness;
    }

    public NodePathSolution getNodePathSolution() {
        return new NodePathSolution(fitness, variables);
    }

    public void evaluate() {
        DijkstraWNDim dj = new DijkstraWNDim(new Long[]{0L, 1L, 2L, 3L});
        dj.setWeights(objectives);
        dj.setGraph(graph);
        List<Node> path = dj.getPath(start, end);
        variables = path;
        fitness = graph.getFitness(path, "");
    }
}
