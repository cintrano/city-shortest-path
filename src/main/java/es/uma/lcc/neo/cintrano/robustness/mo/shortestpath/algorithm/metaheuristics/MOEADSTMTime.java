package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics;

import org.uma.jmetal.algorithm.multiobjective.moead.MOEADSTM;
import org.uma.jmetal.algorithm.multiobjective.moead.util.MOEADUtils;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;

import java.util.List;


/**
 * This is a modification of the JMetal Class MOEADSTM
 * with max time stopping criteria
 *
 * @author Christian Cintrano
 * @version 1.0
 */
@SuppressWarnings("serial")

public class MOEADSTMTime extends MOEADSTM {

    private long initComputingTime;
    private long thresholdComputingTime;

    public MOEADSTMTime(Problem<DoubleSolution> problem, long maxComputingTime, int populationSize, int resultPopulationSize,
                    int maxEvaluations,
                    MutationOperator<DoubleSolution> mutation, CrossoverOperator<DoubleSolution> crossover,
                    FunctionType functionType, String dataDirectory, double neighborhoodSelectionProbability,
                    int maximumNumberOfReplacedSolutions, int neighborSize) {
        super(problem, populationSize, resultPopulationSize, maxEvaluations, mutation, crossover,
                functionType, dataDirectory, neighborhoodSelectionProbability,
                maximumNumberOfReplacedSolutions, neighborSize);
        initComputingTime = System.currentTimeMillis() ;
        thresholdComputingTime = maxComputingTime ;
    }

    @Override
    public void run() {
        initializePopulation();
        initializeUniformWeight();
        initializeNeighborhood();
        idealPoint.update(population);
        nadirPoint.update(population);

        int generation = 0;
        evaluations = populationSize;
        do {
            int[] permutation = new int[populationSize];
            MOEADUtils.randomPermutation(permutation, populationSize);
            offspringPopulation.clear();

            for (int i = 0; i < populationSize; i++) {
                int subProblemId = permutation[i];
                frequency[subProblemId]++;

                NeighborType neighborType = chooseNeighborType();
                List<DoubleSolution> parents = parentSelection(subProblemId, neighborType);

                differentialEvolutionCrossover.setCurrentSolution(population.get(subProblemId));
                List<DoubleSolution> children = differentialEvolutionCrossover.execute(parents);

                DoubleSolution child = children.get(0);
                mutationOperator.execute(child);
                problem.evaluate(child);

                evaluations++;

                idealPoint.update(population);
                nadirPoint.update(population);
                updateNeighborhood(child, subProblemId, neighborType);

                offspringPopulation.add(child);
            }

            // Combine the parent and the current offspring populations
            jointPopulation.clear();
            jointPopulation.addAll(population);
            jointPopulation.addAll(offspringPopulation);

            // selection process
            stmSelection();

            generation++;
            if (generation % 30 == 0) {
                utilityFunction();
            }

        } while (!isStoppingConditionReached());

    }

    private boolean isStoppingConditionReached() {
        long currentComputingTime = System.currentTimeMillis() - initComputingTime ;
        return currentComputingTime > thresholdComputingTime ;
    }

    @Override
    public String getName() {
        return "MOEADSTMTimer";
    }

    @Override
    public String getDescription() {
        return "Multi-Objective Evolutionary Algorithm based on Decomposition. Version with Timer";
    }
}
