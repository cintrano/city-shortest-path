package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2;

import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.JMetalLogger;
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
        List<DoubleSolution> list = new ArrayList<DoubleSolution>();
        if (JMetalRandom.getInstance().nextDouble() < 0.9) {
            List<Double> points = getPoints((MyDoubleSolution) solution1, (MyDoubleSolution) solution2);
            int index = JMetalRandom.getInstance().nextInt(0, points.size()-1);
            int i1 = -1, i2 = -1;
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
            list.add(parentSolutions.get(2));
        }
        return list;
        /*
        DoubleSolution child;

        int jrand;

        child = (DoubleSolution)currentSolution.copy() ;

        int numberOfVariables = parentSolutions.get(0).getNumberOfVariables();
        jrand = jRandomGenerator.getRandomValue(0, numberOfVariables - 1);

        // STEP 4. Checking the DE variant
        if ((DEFAULT_DE_VARIANT.equals(variant)) ||
                "best/1/bin".equals(variant)) {
            for (int j = 0; j < numberOfVariables; j++) {
                if (crRandomGenerator.getRandomValue(0.0, 1.0) < cr || j == jrand) {
                    double value;
                    value = parentSolutions.get(2).getVariableValue(j) + f * (parentSolutions.get(0).getVariableValue(
                            j) -
                            parentSolutions.get(1).getVariableValue(j));

                    if (value < child.getLowerBound(j)) {
                        value = child.getLowerBound(j);
                    }
                    if (value > child.getUpperBound(j)) {
                        value = child.getUpperBound(j);
                    }
                    child.setVariableValue(j, value);
                } else {
                    double value;
                    value = currentSolution.getVariableValue(j);
                    child.setVariableValue(j, value);
                }
            }
        } else if ("rand/1/exp".equals(variant) ||
                "best/1/exp".equals(variant)) {
            for (int j = 0; j < numberOfVariables; j++) {
                if (crRandomGenerator.getRandomValue(0.0, 1.0) < cr || j == jrand) {
                    double value;
                    value = parentSolutions.get(2).getVariableValue(j) + f * (parentSolutions.get(0).getVariableValue(j) -
                            parentSolutions.get(1).getVariableValue(j));

                    if (value < child.getLowerBound(j)) {
                        value = child.getLowerBound(j);
                    }
                    if (value > child.getUpperBound(j)) {
                        value = child.getUpperBound(j);
                    }

                    child.setVariableValue(j, value);
                } else {
                    cr = 0.0;
                    double value;
                    value = currentSolution.getVariableValue(j);
                    child.setVariableValue(j, value);
                }
            }
        } else if ("current-to-rand/1".equals(variant) ||
                "current-to-best/1".equals(variant)) {
            for (int j = 0; j < numberOfVariables; j++) {
                double value;
                value = currentSolution.getVariableValue(j) + k * (parentSolutions.get(2).getVariableValue(j) -
                        currentSolution.getVariableValue(j)) +
                        f * (parentSolutions.get(0).getVariableValue(j) - parentSolutions.get(1).getVariableValue(j));

                if (value < child.getLowerBound(j)) {
                    value = child.getLowerBound(j);
                }
                if (value > child.getUpperBound(j)) {
                    value = child.getUpperBound(j);
                }

                child.setVariableValue(j, value);
            }
        } else if ("current-to-rand/1/bin".equals(variant) ||
                "current-to-best/1/bin".equals(variant)) {
            for (int j = 0; j < numberOfVariables; j++) {
                if (crRandomGenerator.getRandomValue(0.0, 1.0) < cr || j == jrand) {
                    double value;
                    value = currentSolution.getVariableValue(j) + k * (parentSolutions.get(2).getVariableValue(j) -
                            currentSolution.getVariableValue(j)) +
                            f * (parentSolutions.get(0).getVariableValue(j) - parentSolutions.get(1).getVariableValue(j));

                    if (value < child.getLowerBound(j)) {
                        value = child.getLowerBound(j);
                    }
                    if (value > child.getUpperBound(j)) {
                        value = child.getUpperBound(j);
                    }

                    child.setVariableValue(j, value);
                } else {
                    double value;
                    value = currentSolution.getVariableValue(j);
                    child.setVariableValue(j, value);
                }
            }
        } else if ("current-to-rand/1/exp".equals(variant) ||
                "current-to-best/1/exp".equals(variant)) {
            for (int j = 0; j < numberOfVariables; j++) {
                if (crRandomGenerator.getRandomValue(0.0, 1.0) < cr || j == jrand) {
                    double value;
                    value = currentSolution.getVariableValue(j) + k * (parentSolutions.get(2).getVariableValue(j) -
                            currentSolution.getVariableValue(j)) +
                            f * (parentSolutions.get(0).getVariableValue(j) - parentSolutions.get(1).getVariableValue(j));

                    if (value < child.getLowerBound(j)) {
                        value = child.getLowerBound(j);
                    }
                    if (value > child.getUpperBound(j)) {
                        value = child.getUpperBound(j);
                    }

                    child.setVariableValue(j, value);
                } else {
                    cr = 0.0;
                    double value;
                    value = currentSolution.getVariableValue(j);
                    child.setVariableValue(j, value);
                }
            }
        } else {
            JMetalLogger.logger.severe("DifferentialEvolutionCrossover.execute: " +
                    " unknown DE variant (" + variant + ")");
            Class<String> cls = String.class;
            String name = cls.getName();
            throw new JMetalException("Exception in " + name + ".execute()");
        }

        List<DoubleSolution> result = new ArrayList<>(1) ;
        result.add(child) ;
        return result;
        */
    }

    private static List<Double> getPoints(MyDoubleSolution solution1, MyDoubleSolution solution2) {
        List<Double> solutions = new ArrayList<Double>();

        for (int i = 0; i < solution1.getNumberOfVariables(); i++) {
            //System.out.println(solution1.getVariableValue(i) + " " + Arrays.asList(solution2.getVariables()).contains(solution1.getVariableValue(i)));
            if (Arrays.asList(solution2.getVariables()).contains(solution1.getVariableValue(i))) {
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
