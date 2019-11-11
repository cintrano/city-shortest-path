package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.operators;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.MOShortestPathProblemDouble;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.MyDoubleSolution;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.util.ArrayList;

/**
 * Created by Christian Cintrano on 14/12/16.
 * Mutation operation: swap and new individual
 */
public class Mutation1PChangeD implements MutationOperator<MyDoubleSolution> {
    private double mutationProbability;

    private JMetalRandom randomGenerator ;

    private MOShortestPathProblemDouble problem;

    /** Constructor */
    public Mutation1PChangeD(MOShortestPathProblemDouble problem, double mutationProbability, double distributionIndex) {
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
            int index = randomGenerator.nextInt(0, problem.getGraph().getIntersections().keySet().size()-1);
            Long nextNode = new ArrayList<>(problem.getGraph().getIntersections().keySet()).get(index);
            Long[] head = problem.getShortestPathBetween( (pathSolution.getVariableValue(0)).longValue(), nextNode);
            Long[] tail = problem.getShortestPathBetween(nextNode, pathSolution.getVariableValue(pathSolution.getNumberOfVariables() - 1).longValue());

            Double[] out;
            if (tail.length == 0 || head.length == 0) {
                out = new Double[pathSolution.getNumberOfVariables()];
                for (int i = 0; i < out.length; i++) {
                    out[i] = pathSolution.getVariableValue(i);
                }
            } else {
                out = new Double[head.length + tail.length - 1];
                for (int i = 0; i < head.length; i++) {
                    out[i] = head[i].doubleValue();
                }
                for (int i = 1; i < tail.length; i++) {
                    out[i - 1 + head.length] = tail[i].doubleValue();
                }
            }
            pathSolution.setVariables(out);
        }

        return pathSolution;
    }


}
