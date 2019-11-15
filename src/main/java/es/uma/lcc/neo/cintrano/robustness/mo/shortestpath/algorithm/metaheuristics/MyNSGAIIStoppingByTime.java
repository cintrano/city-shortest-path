package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics;

import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAII;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.problem.ConstrainedProblem;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;

import java.util.Comparator;
import java.util.List;

/**
 * This class shows a version of NSGA-II having a stopping condition depending on run-time
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
@SuppressWarnings("serial")
public class MyNSGAIIStoppingByTime <S extends Solution<?>> extends NSGAII<S> {
    private long initComputingTime;
    private long thresholdComputingTime;
    private boolean stoppingCondition;
    /**
     * Constructor
     */
    public MyNSGAIIStoppingByTime(Problem<S> problem, int populationSize,
                                  long maxComputingTime,
                                  CrossoverOperator<S> crossoverOperator, MutationOperator<S> mutationOperator,
                                  SelectionOperator<List<S>, S> selectionOperator, Comparator<S> dominanceComparator,
                                  SolutionListEvaluator<S> evaluator) {
        super(problem, 0, populationSize, crossoverOperator, mutationOperator,
                selectionOperator, dominanceComparator, evaluator);

        initComputingTime = System.currentTimeMillis();
        stoppingCondition = false;
        thresholdComputingTime = maxComputingTime;
    }

//    @Override
//    protected void initProgress() {
//        evaluations = getMaxPopulationSize();
//    }

//    @Override
//    protected List<S> evaluatePopulation(List<S> population) {
//        int index = 0;
//
//        while ((index < population.size()) && !stoppingCondition) {
//            if (getProblem() instanceof ConstrainedProblem) {
//                getProblem().evaluate(population.get(index));
//                ((ConstrainedProblem<S>) getProblem()).evaluateConstraints(population.get(index));
//            } else {
//                getProblem().evaluate(population.get(index));
//            }
//
//            if ((System.currentTimeMillis() - initComputingTime) > thresholdComputingTime) {
//                stoppingCondition = true;
//            } else {
//                evaluations++;
//                index++;
//            }
//        }
//
//        return population;
//    }

    @Override protected boolean isStoppingConditionReached() {
        return (System.currentTimeMillis() - initComputingTime) > thresholdComputingTime;
//        return stoppingCondition;
    }

    @Override public String getName() {
        return "MyNSGAII by time";
    }

    @Override public String getDescription() {
        return "MyNSGAII by time";
    }
}
