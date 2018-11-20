package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.operators;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MOShortestPathProblem;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

/**
 * Created by Christian Cintrano on 14/12/16.
 * Mutation operation: swap and new individual
 */
public class Mutation1PRegeneration implements MutationOperator<NodePathSolution> {
    private double mutationProbability;

    private JMetalRandom randomGenerator ;

    private MOShortestPathProblem problem;

    /** Constructor */
    public Mutation1PRegeneration(MOShortestPathProblem problem, double mutationProbability, double distributionIndex) {
        if (mutationProbability < 0) {
            throw new JMetalException("Crossover probability is negative: " + mutationProbability) ;
        } else if (distributionIndex < 0) {
            throw new JMetalException("Distribution index is negative: " + distributionIndex);
        }

        this.mutationProbability =mutationProbability ;

        this.problem = problem;

        randomGenerator = JMetalRandom.getInstance();
    }

    public NodePathSolution execute(NodePathSolution pathSolution) {
        /*
        if (randomGenerator.nextDouble() < mutationProbability) {
            int i = randomGenerator.nextInt(0, pathSolution.getNumberOfVariables() - 1); // NOT we make sure that the new individual is different
            Long t = pathSolution.getVariableValue(i);
            Long[] partialPath = problem.getShortestPartialPath(t);

            Long[] first = Arrays.copyOfRange(pathSolution.getVariables(),0, i);
            Long[] variables = ArrayUtils.addAll(first, partialPath);

            pathSolution.setVariables(variables);
        }
        */
        return pathSolution;
    }


}
