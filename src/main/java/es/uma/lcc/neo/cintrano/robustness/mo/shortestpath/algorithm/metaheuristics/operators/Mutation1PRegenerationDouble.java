package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.operators;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.RunTLMain;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDimForbidden;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.MOShortestPathProblem;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.MOShortestPathProblemDouble;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.MyDoubleSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.NodePathSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import org.apache.commons.lang3.ArrayUtils;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Created by Christian Cintrano on 14/12/16.
 * Mutation operation: swap and new individual
 */
public class Mutation1PRegenerationDouble implements MutationOperator<MyDoubleSolution> {
    private double mutationProbability;

    private JMetalRandom randomGenerator ;

    private MOShortestPathProblemDouble problem;

    /** Constructor */
    public Mutation1PRegenerationDouble(MOShortestPathProblemDouble problem, double mutationProbability, double distributionIndex) {
        if (mutationProbability < 0) {
            throw new JMetalException("Crossover probability is negative: " + mutationProbability) ;
        } else if (distributionIndex < 0) {
            throw new JMetalException("Distribution index is negative: " + distributionIndex);
        }

        this.mutationProbability =mutationProbability ;

        this.problem = problem;

        randomGenerator = JMetalRandom.getInstance();
    }

    public MyDoubleSolution execute(MyDoubleSolution pathSolution) {
        if (randomGenerator.nextDouble() < mutationProbability) {
            Set<Long> neighbours;
            int point;
            do {
                point = randomGenerator.nextInt(0, pathSolution.getNumberOfVariables() - 1);
                neighbours = problem.getGraph().getAdjacencyMatrix().row(pathSolution.getVariableValue(point).longValue()).keySet();
            } while (neighbours.size() < 2);
            //List<Node> newSolution = new ArrayList<>(s.subList(0, point + 1));
            Double[] head = Arrays.copyOfRange(pathSolution.getVariables(), 0, point + 1);

            // Dijkstra algorithm helper
            // TODO: Change the implementation of Dijkstra to take into account a set of forbidden nodes
            DijkstraWNDimForbidden dj = new DijkstraWNDimForbidden(RunTLMain.objectives);
            dj.setWeights(new float[]{1, 1, 1, 1});//dj.setWeights(weights);
            dj.setGraph(problem.getGraph());

            //dj.setForbiddenList(Arrays.stream(head).mapToLong(Double::longValue).boxed().collect(Collectors.toList()));

            // Select the new node
            List<Long> neighboursList = new ArrayList<>(neighbours);
            Collections.shuffle(neighboursList);
            int j = 0;
            while (j < neighboursList.size() && neighboursList.get(j).equals(pathSolution.getVariableValue(point + 1).longValue())) {
                j++;
            }
            List<Node> subPath = dj.getPath(neighboursList.get(j), pathSolution.getVariableValue(pathSolution.getNumberOfVariables() - 1).longValue());
            Double[] aux = new Double[head.length + subPath.size()];
            System.arraycopy(head, 0, aux, 0, head.length);
            for (int i = 0; i < subPath.size(); i++) {
                aux[i + head.length] = (double) problem.getGraph().getMapping().get(subPath.get(i).getId());
            }
            pathSolution.setVariables(aux);
        }
        return pathSolution;
    }
}
