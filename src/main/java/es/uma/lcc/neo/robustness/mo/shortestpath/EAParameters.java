package es.uma.lcc.neo.robustness.mo.shortestpath;

import es.uma.lcc.neo.robustness.mo.shortestpath.eaparameters.utilities.EASolution;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;

import java.util.List;
import java.util.Random;

/**
 * Created by Christian Cintrano on 18/05/17.
 * Evolutive Algorithm to calculate the best weights to transform the multi- to mono- objective
 */
public class EAParameters {
    private final static int MAX_ITERATIONS = 625; // 5^4 = (0, .25, .5, .75, 1)^objectives
    private final static int OBJECTIVE_SIZE = 4;

    private GraphTable graph;
    private Long start, end;
    private float mutationProbability;

    public void runEA(List<EASolution> population, Random random) {
        int it = 0;
        population = generatePopulation();
        population = evaluatePopulation(population);
        while (it < MAX_ITERATIONS) {
            EASolution[] parents = selection(population, random);
            EASolution child = crossover(parents, it, random);
            child = mutation(child, random);
            child.evaluate();
            replace(population, child);

            it++;
        }
    }

    /**
     * Random Selection
     * @param population Current population
     * @param random Random
     * @return Array with the two parents
     */
    private EASolution[] selection(List<EASolution> population, Random random) {
        int j, i = random.nextInt(population.size());
        do {
            j = random.nextInt(population.size());
        } while (j == i);
        return new EASolution[]{population.get(i), population.get(j)};
    }

    /**
     * Crossover 1-point cut
     * @param parents Array with the two parents
     * @param random Random
     * @return New individual
     */
    private EASolution crossover(EASolution[] parents, int iteration, Random random) {
        int cut = random.nextInt(OBJECTIVE_SIZE - 1) + 1; // cut = {1,2,3}
        EASolution individual = new EASolution(iteration, graph, start, end);
        float[] objectives = new float[OBJECTIVE_SIZE];
        System.arraycopy(parents[1].getFitness(), 0, objectives, 0, cut);
        System.arraycopy(parents[2].getFitness(), cut, objectives, cut, OBJECTIVE_SIZE - cut);
        individual.setObjectives(objectives);
        return individual;
    }

    /**
     * Mutation one objective
     * @param child individual
     * @param random Random
     * @return Individual mutated
     */
    private EASolution mutation(EASolution child, Random random) {
        float prob = random.nextFloat();
        if (prob > mutationProbability) {
            int i = random.nextInt(OBJECTIVE_SIZE);
            child.getFitness()[i] = random.nextFloat();
        }
        return child;
    }

    private void replace(List<EASolution> population, EASolution child) {
        
    }

}
