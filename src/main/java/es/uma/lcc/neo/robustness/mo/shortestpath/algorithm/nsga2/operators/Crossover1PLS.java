package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.operators;

import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution;
import org.apache.commons.lang3.ArrayUtils;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by cintrano on 22/12/16.
 * Crossover one point + local search
 */
public class Crossover1PLS implements CrossoverOperator<NodePathSolution> {
    private double crossoverProbability;

    private JMetalRandom randomGenerator;


    /** Constructor */
    public Crossover1PLS(double crossoverProbability, double distributionIndex) {
        if (crossoverProbability < 0) {
            throw new JMetalException("Crossover probability is negative: " + crossoverProbability);
        } else if (distributionIndex < 0) {
            throw new JMetalException("Distribution index is negative: " + distributionIndex);
        }

        this.crossoverProbability = crossoverProbability;
        randomGenerator = JMetalRandom.getInstance();
    }
    public List<NodePathSolution> execute(List<NodePathSolution> solutions){
        if (null == solutions) {
            throw new JMetalException("Null parameter") ;
        } else if (solutions.size() != 2) {
            throw new JMetalException("There must be two parents instead of " + solutions.size()) ;
        }

        return doCrossover(crossoverProbability, solutions.get(0), solutions.get(1)) ;
    }

    private List<NodePathSolution> doCrossover(double crossoverProbability, NodePathSolution solution1, NodePathSolution solution2) {
        // Crossover
        List<NodePathSolution> list = new ArrayList<NodePathSolution>();
        if (randomGenerator.nextDouble() < crossoverProbability) {
            List<Long> points = getPoints(solution1, solution2);
            int index = randomGenerator.nextInt(0, points.size()-1);
            int i1 = -1, i2 = -1;
            for (int j = 0; j < solution1.getVariables().length; j++) {
                if (points.get(index).equals(solution1.getVariableValue(j))) {
                    i1 = j;
                    break;
                }
            }
            for (int j = 0; j < solution2.getVariables().length; j++) {
                if (points.get(index).equals(solution2.getVariableValue(j))) {
                    i2 = j;
                    break;
                }
            }

            Long[] first = new Long[i1 + solution2.getVariables().length - i2];
            for (int j = 0; j < i1; j++) {
                first[j] = solution1.getVariableValue(j);
            }
            for (int j = i2; j < solution2.getVariables().length; j++) {
                first[j-i2+i1] = solution2.getVariableValue(j);
            }

            Long[] second = new Long[i2 + solution1.getVariables().length - i1];
            for (int j = 0; j < i2; j++) {
                second[j] = solution2.getVariableValue(j);
            }
            for (int j = i1; j < solution1.getVariables().length; j++) {
                second[j-i1+i2] = solution1.getVariableValue(j);
            }


            NodePathSolution child1 = new NodePathSolution(new double[solution1.getNumberOfObjectives()], first);
            NodePathSolution child2 = new NodePathSolution(new double[solution1.getNumberOfObjectives()], second);



            // Local Search


            // END
            list.add(child1);
            list.add(child2);
        } else {
            list.add(solution1);
            list.add(solution2);
        }
        return list;
    }

    private static List<Long> getPoints(NodePathSolution solution1, NodePathSolution solution2) {
        List<Long> solutions = new ArrayList<Long>();

        for (int i = 0; i < solution1.getNumberOfVariables(); i++) {
            //System.out.println(solution1.getVariableValue(i) + " " + Arrays.asList(solution2.getVariables()).contains(solution1.getVariableValue(i)));
            if (Arrays.asList(solution2.getVariables()).contains(solution1.getVariableValue(i))) {
                solutions.add(solution1.getVariableValue(i));
            }
        }
        return solutions;
    }
}
