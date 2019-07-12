package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.operators;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MOShortestPathProblem;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MOShortestPathProblemDouble;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDoubleSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.util.ArrayList;

/**
 * Created by Christian Cintrano on 14/12/16.
 * Mutation operation: swap and new individual
 */
public class Mutation1PChangeD implements MutationOperator<DoubleSolution> {
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

    public DoubleSolution execute(DoubleSolution pathSolution) {

        if (randomGenerator.nextDouble() < mutationProbability) {
            ///////////////////////////

            int index = randomGenerator.nextInt(0, problem.getGraph().getIntersections().keySet().size()-1);
            Long nextNode = new ArrayList<Long>(problem.getGraph().getIntersections().keySet()).get(index);
            Long[] head = problem.getShortestPathBetween( (pathSolution.getVariableValue(0)).longValue(), nextNode);
            Long[] tail = problem.getShortestPathBetween(nextNode, pathSolution.getVariableValue(pathSolution.getNumberOfVariables() - 1).longValue());

            Long[] out;
            if (tail.length == 0 || head.length == 0) {
                out = new Long[pathSolution.getNumberOfVariables()];
                for (int i = 0; i < out.length; i++) {
                    out[i] = pathSolution.getVariableValue(i).longValue();
                }
                //out = pathSolution.getVariables();
            } else {
                out = new Long[head.length + tail.length - 1];
                for (int i = 0; i < head.length; i++) {
                    out[i] = head[i];
                }
                for (int i = 1; i < tail.length; i++) {
                    out[i + head.length - 1] = tail[i];
                }
            }

            ((MyDoubleSolution) pathSolution).setVariables(new Double[out.length]);
            for (int i = 0; i < out.length; i++) {
                pathSolution.setVariableValue(i, out[i].doubleValue());
            }

            ///////////////////////////
            /*
            int index = randomGenerator.nextInt(0, pathSolution.getVariables().length - 2);
            Long nextNode = problem.findRandomNextNode(pathSolution.getVariableValue(index));
            Long[] tail = problem.getShortestPathWithMiddlePoint(nextNode);

            Long[] out;
            if (tail.length == 0) {
                out = pathSolution.getVariables();
            } else {
                out = new Long[index + tail.length + 1];
                for (int i = 0; i <= index; i++) {
                    out[i] = pathSolution.getVariableValue(i);
                }
                for (int i = 0; i < tail.length; i++) {
                    out[i + index + 1] = tail[i];
                }
            }
            pathSolution.setVariables(out);
            */
        }

        return pathSolution;
    }


}
