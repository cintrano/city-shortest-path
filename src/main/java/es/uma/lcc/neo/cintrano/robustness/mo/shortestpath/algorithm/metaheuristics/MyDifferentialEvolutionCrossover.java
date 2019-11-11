package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics;

import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.pseudorandom.BoundedRandomGenerator;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;
import org.uma.jmetal.util.pseudorandom.RandomGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MyDifferentialEvolutionCrossover extends DifferentialEvolutionCrossover {
    private static final double DEFAULT_CR = 0.5;
    private static final double DEFAULT_F = 0.5;
    private static final double DEFAULT_K = 0.5;
    private static final String DEFAULT_DE_VARIANT = "rand/1/bin";

    private double cr;
    private double f;
    private double k;
    // DE variant (rand/1/bin, rand/1/exp, etc.)
    private String variant;

    private DoubleSolution currentSolution ;

    private BoundedRandomGenerator<Integer> jRandomGenerator ;
    private BoundedRandomGenerator<Double> crRandomGenerator ;

    /** Constructor */
    public MyDifferentialEvolutionCrossover() {
        this(DEFAULT_CR, DEFAULT_F, DEFAULT_K, DEFAULT_DE_VARIANT) ;
    }

    /**
     * Constructor
     * @param cr
     * @param f
     * @param variant
     */
    public MyDifferentialEvolutionCrossover(double cr, double f, String variant) {
        this(cr, f, variant, (a, b) -> JMetalRandom.getInstance().nextInt(a, b), (a, b) -> JMetalRandom.getInstance().nextDouble(a, b));
    }


    public MyDifferentialEvolutionCrossover(double cr, double f, String variant, RandomGenerator<Double> randomGenerator) {
        this(cr, f, variant, BoundedRandomGenerator.fromDoubleToInteger(randomGenerator), BoundedRandomGenerator.bound(randomGenerator));
    }

    /**
     * Constructor
     * @param cr
     * @param f
     * @param variant
     * @param jRandomGenerator
     * @param crRandomGenerator
     */
    public MyDifferentialEvolutionCrossover(double cr, double f, String variant, BoundedRandomGenerator<Integer> jRandomGenerator, BoundedRandomGenerator<Double> crRandomGenerator) {
        this.cr = cr;
        this.f = f;
        this.k = DEFAULT_K ;
        this.variant = variant ;

        this.jRandomGenerator = jRandomGenerator;
        this.crRandomGenerator = crRandomGenerator ;
    }

    /** Constructor */
    public MyDifferentialEvolutionCrossover(double cr, double f, double k, String variant) {
        this(cr, f, variant) ;
        this.k = k ;
    }

    /* Getters */
    public double getCr() {
        return cr;
    }

    public double getF() {
        return f;
    }

    public double getK() {
        return k;
    }

    public String getVariant() {
        return variant;
    }

    /* Setters */
    public void setCurrentSolution(DoubleSolution current) {
        this.currentSolution = current ;
    }

    public void setCr(double cr) {
        this.cr = cr;
    }

    public void setF(double f) {
        this.f = f;
    }

    public void setK(double k) {
        this.k = k;
    }

    /** Execute() method */
    @Override
    public List<DoubleSolution> execute(List<DoubleSolution> parentSolutions) {
        // Crossover
        DoubleSolution solution1 = parentSolutions.get(0);
        DoubleSolution solution2 = parentSolutions.get(1);
        List<DoubleSolution> list = new ArrayList<>();
        if (JMetalRandom.getInstance().nextDouble() < 0.9) {
            List<Double> points = getPoints((MyDoubleSolution) solution1, (MyDoubleSolution) solution2);//getPoints(solution1, solution2);
            int index = JMetalRandom.getInstance().nextInt(0, points.size()-1);
            int i1 = -1, i2 = -1; // index in the solutions
            for (int j = 0; j < solution1.getNumberOfVariables(); j++) {
                if (points.get(index).equals(solution1.getVariableValue(j))) {
                    i1 = j;
                    break;
                }
            }
            for (int j = 0; j < solution2.getNumberOfVariables(); j++) {
                if (points.get(index).equals(solution2.getVariableValue(j))) {
                    i2 = j;
                    break;
                }
            }

            Double[] first = new Double[i1 + solution2.getNumberOfVariables() - i2];
            for (int j = 0; j < i1; j++) {
                first[j] = solution1.getVariableValue(j);
            }
            for (int j = i2; j < solution2.getNumberOfVariables(); j++) {
                first[j-i2+i1] = solution2.getVariableValue(j);
            }

            Double[] second = new Double[i2 + solution1.getNumberOfVariables() - i1];
            for (int j = 0; j < i2; j++) {
                second[j] = solution2.getVariableValue(j);
            }
            for (int j = i1; j < solution1.getNumberOfVariables(); j++) {
                second[j-i1+i2] = solution1.getVariableValue(j);
            }


            MyDoubleSolution child1 = new MyDoubleSolution(new double[solution1.getNumberOfObjectives()], first);
            MyDoubleSolution child2 = new MyDoubleSolution(new double[solution1.getNumberOfObjectives()], second);

            // Local Search


            // END
            list.add(child1);
            //list.add(child2);
        } else {
            //list.add(solution1);
            //list.add(solution2);
            list.add((MyDoubleSolution) parentSolutions.get(2).copy());
        }
        return list;
    }

    private static List<Double> getPoints(MyDoubleSolution solution1, MyDoubleSolution solution2) {
        List<Double> solutions = new ArrayList<>();

        for (int i = 0; i < solution1.getNumberOfVariables(); i++) {
            if (Arrays.asList((solution2).getVariables()).contains(solution1.getVariableValue(i))) {
                solutions.add(solution1.getVariableValue(i));
            }
        }
        return solutions;
    }

    public int getNumberOfRequiredParents() {
        return 3 ;
    }

    public int getNumberOfGeneratedChildren() {
        return 1 ;
    }
}
