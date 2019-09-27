package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.astar.AstarMO;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDim;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.iterated.IteratedLS;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.iterated.MOIteratedLS;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MOShortestPathProblem;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MOShortestPathProblemDouble;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDifferentialEvolutionCrossover;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.operators.Crossover1PLS;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.operators.Mutation1PChange;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.operators.Mutation1PChangeD;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.pulse.Pulse;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.pulse.PulseMO;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.utilities.ProcessGraph;
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEAD;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEADSTM;
import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAII;
import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIMeasures;
import org.uma.jmetal.measure.MeasureListener;
import org.uma.jmetal.measure.MeasureManager;
import org.uma.jmetal.measure.PushMeasure;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.comparator.DominanceComparator;
import org.uma.jmetal.util.comparator.RankingAndCrowdingDistanceComparator;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static org.uma.jmetal.util.AbstractAlgorithmRunner.printFinalSolutionSet;

/**
 * Created by Christian Cintrano on 22/01/17.
 * Main Class
 */
public class RunTLMain {

    private final static int[] seed = new int[]{0, 3, 23, 29, 67, 79, 97, 101, 139, 199, 211, 269, 311, 347, 389,
            431, 461, 503, 569, 601, 607, 653, 691, 701, 733, 739, 761, 809, 811, 997};

    public final static Long[] objectives = new Long[]{0L, 1L, 2L, 3L};
    private final static float[] weights = new float[]{0.0001f, 0.25f, 0.5f, 0.75f, 1.0f};

    public static void main (String[] args) {
        System.out.println("=== START EXPERIMENTS ===");
        for (String s : args) {
            System.out.print(s + " ");
        }
        System.out.println();

        if (args.length >= 3) {
            GraphTable graph = null;
            // Map
            switch (args[0]) {
                case "Malaga":
                    System.out.println("Loading graph...");
                    graph = ProcessGraph.prepareGraph("g_CC.xml", "w_filter.xml", "malaga.osm-mapping.txt", "MAL");
                    ProcessGraph.fixVertexIndex(graph);
                    ProcessGraph.readTlLogics(graph, "malaga_tls.json");
                    break;
                case "Colorado":
                    graph = ProcessGraph.prepareGraph("hbefa-col-graph.xml", "col_weights_time-hbefa-MOD.xml", "mapping-col.txt", "COL");
                    break;
                case "NY":
                    graph = ProcessGraph.prepareGraph("hbefa-ny-graph.xml", "ny_weights_time-hbefa.xml", "mapping-ny.txt", "NY");
                    break;
                case "Graph":
                    //graph = ProcessGraph.prepareGraph("malga.osm-tlgraph.xml", "malaga.osm-tlweight.xml", "malaga.osm-mapping.txt", "MAL");
                    graph = ProcessGraph.prepareGraph("g_filter.xml", "w_filter.xml", "malaga.osm-mapping.txt", "MAL");
                    ProcessGraph.printMapping(graph);
                    ProcessGraph.getConnectedComponents(graph);
                    ProcessGraph.printSCCs(graph);
                    ProcessGraph.getMaxConnectedComponent(graph, "g_map.txt");
                    ProcessGraph.printGraph(graph, "g-CC.xml");
                    ProcessGraph.getConnectedComponents(graph);
                    ProcessGraph.printSCCs(graph);
                    break;
            }

            // Points
            System.out.print("Reading points...");
            Long[] points;
            if (!args[2].equals("default")) {
                points = MyUtility.readPoints(Integer.parseInt(args[2]), args[0]);
            } else {
                points = new Long[]{8L, 720L};
            }

            // Algorithm
            if (points != null) {
                System.out.println(points[0] + " " + points[1]);
                //boolean robust = Boolean.parseBoolean(args[4]);
                //boolean tlOption = Boolean.parseBoolean(args[5]);
                boolean robust = args[4].equals("Robust");
                boolean tlOption = args[5].equals("TL");
                System.out.println("FLAGs:: Robust=" + robust + ", TL=" + tlOption);

                switch (args[1]) {
                    case "PULSE4":
                        System.out.println("Run PULSE4...");
                        rPulse(graph, points);
                        break;
                    case "PULSE2":
                        System.out.println("Run PULSE2...");
                        rPulse2Experiments(graph, points);
                        break;
                    case "Dijkstra":
                        System.out.println("Run Dijkstra...");
                        rDijkstraExperiments(graph, points);
                        break;
                    case "Astar":
                        System.out.println("Run A*...");
                        rAstarExperiments(graph, points);
                        break;
                    case "NSGAII":
                        System.out.println("Run NSGA-II...");
                        rNSGAII(graph, points, seed[Integer.parseInt(args[3])], robust, tlOption);
                        break;
                    case "MOEAD":
                        System.out.println("Run MOEA-D...");
                        rMOEAD(graph, points, seed[Integer.parseInt(args[3])], robust, tlOption);
                        break;
                    case "NSGAIII":
                        System.out.println("Run NSGA-III...");
                        rNSGAIII(graph, points, seed[Integer.parseInt(args[3])]);
                        break;
                    case "Iterated":
                        System.out.println("Run Iterated...");
                        rIterated(graph, points, seed[Integer.parseInt(args[3])]);
                        break;
                    case "MOIterated":
                        System.out.println("Run MOIterated...");
                        rMOIterated(graph, points, seed[Integer.parseInt(args[3])], robust, tlOption);
                        break;
                }
            }
        } else if (args[0].equals("Mapping")) {
            GraphTable graph = ProcessGraph.parserFile(args[1]);
            assert graph != null;
            ProcessGraph.printMapping(graph);
        }

        System.out.println("=== END EXPERIMENTS ===");
    }


    private static void rIterated(GraphTable graph, Long[] points, int seed) {
        long timeInit = System.currentTimeMillis();
        IteratedLS algorithm = new IteratedLS(seed);
        algorithm.setGraph(graph);
        long timeStart = System.currentTimeMillis();
        List<Node> path = algorithm.getPath(points[0], points[1], new float[]{1, 1, 1, 1});
        long timeEnd = System.currentTimeMillis();
        String wTag = "1-1-1-1";
        MyUtility.printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "Iterated", wTag, seed);
    }


    private static void rMOIterated(GraphTable graph, Long[] points, int seed, boolean robust, boolean tlOption) {
        long timeInit = System.currentTimeMillis();
        MOIteratedLS algorithm = new MOIteratedLS(seed);
        algorithm.setGraph(graph);
        algorithm.setRobustFlag(robust);
        algorithm.setTLFlag(tlOption);
        long timeStart = System.currentTimeMillis();
        Set<NodePathSolution> paths = algorithm.getPath(points[0], points[1], new float[]{1, 1, 1, 1});
        long timeEnd = System.currentTimeMillis();
        String wTag = "1-1-1-1";
        MyUtility.printSolutions(graph, paths, timeStart - timeInit, timeEnd - timeStart, "MOIterated", wTag, seed);

    }

    private static void rDijkstraExperiments(GraphTable graph, Long[] randomPoints) {
        for (float w0 : weights) {
            for (float w1 : weights) {
                for (float w2 : weights) {
                    for (float w3 : weights) {
                        rDijkstraWD(graph, randomPoints, new float[]{w0, w1, w2, w3});
                    }
                }
            }
        }
    }

    private static void rDijkstraWD(GraphTable graph, Long[] randomPoints, float[] weights) {
        long timeInit = System.currentTimeMillis();
        DijkstraWNDim dj = new DijkstraWNDim(RunTLMain.objectives);
        dj.setWeights(weights);
        dj.setGraph(graph);
        long timeStart = System.currentTimeMillis();
        List<Node> path = dj.getPath(randomPoints[0], randomPoints[1]);
        long timeEnd = System.currentTimeMillis();
        String wTag = weights[0] + "-" +weights[1] + "-" +weights[2] + "-" +weights[3];
        MyUtility.printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "Dijkstra", wTag);
    }

    private static void rAstarExperiments(GraphTable graph, Long[] randomPoints) {
        float[] weights = new float[]{0.0001f, 0.25f, 0.5f, 0.75f, 1.0f};

        for (float w0 : weights) {
            for (float w1 : weights) {
                for (float w2 : weights) {
                    for (float w3 : weights) {
                        rAstarMO(graph, randomPoints, new float[]{w0, w1, w2, w3});
                    }
                }
            }
        }
    }

    private static void rAstarMO(GraphTable graph, Long[] randomPoints, float[] weights) {
        long timeInit = System.currentTimeMillis();
        AstarMO astar = new AstarMO();
        astar.setWeights(weights);
        astar.setGraph(graph);
        astar.setTarget(1F);
        long timeStart = System.currentTimeMillis();
        List<Node> path = astar.getPath(randomPoints[0], randomPoints[1]);
        long timeEnd = System.currentTimeMillis();
        String wTag = weights[0] + "-" +weights[1] + "-" +weights[2] + "-" +weights[3];
        MyUtility.printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "A*", wTag);
    }

    private static void rPulse(GraphTable graph, Long[] randomPoints) {
        long timeInit = System.currentTimeMillis();
        PulseMO pulse = new PulseMO(objectives.length);
        pulse.setGraph(graph);

        long timeStart = System.currentTimeMillis();
        Set<NodePathSolution> solutions = pulse.pulseAlgorithm(randomPoints[0], randomPoints[1]);
        long timeEnd = System.currentTimeMillis();

        MyUtility.printSummary(timeInit, timeStart, solutions, timeEnd);
        MyUtility.printSolutions(solutions, (timeStart - timeInit), (timeEnd - timeStart), "PULSE4");
    }

    private static void rPulse2Experiments(GraphTable graph, Long[] randomPoints) {
        for (int i = 0; i < objectives.length - 1; i++) {
            for (int j = i + 1; j < objectives.length; j++) {
                System.out.println(objectives[i] + " " + objectives[j]);
                rPulse2(graph, randomPoints, objectives[i], objectives[j]);
            }
        }
    }

    private static void rPulse2(GraphTable graph, Long[] randomPoints, Long obj1, Long obj2) {
        long timeInit = System.currentTimeMillis();
        Pulse pulse = new Pulse(obj1, obj2);
        pulse.setGraph(graph);

        long timeStart = System.currentTimeMillis();
        Set<NodePathSolution> solutions = pulse.pulseAlgorithm(randomPoints[0], randomPoints[1]);
        long timeEnd = System.currentTimeMillis();

        MyUtility.printSummary(timeInit, timeStart, solutions, timeEnd);
        MyUtility.printSolutions(solutions, (timeStart - timeInit), (timeEnd - timeStart),
                "PULSE2", obj1 + "-" + obj2);
    }



    private static Map<Long, List<double[]>> fitnessAll = new HashMap<>();
    private static Map<Long, List<Double[]>> solutions = new HashMap<>();
    private static void rNSGAII(GraphTable graph, Long[] points, long seed, boolean robust, boolean tlOption) {
        executeAlgorithm(
                new MOShortestPathProblem(graph, points[0], points[1], 10, robust, tlOption),
                seed,
                0.9, //crossoverProbability,
                0.1, //mutationProbability,
                10000,//numIterations,
                16,//populationSize,
                "NSGAII"
        );
    }

    private static void rMOEAD(GraphTable graph, Long[] points, long seed, boolean robust, boolean tloption) {
        executeAlgorithm(
                new MOShortestPathProblemDouble(graph, points[0], points[1], 10, robust, tloption),
                seed,
                0.9, //crossoverProbability,
                0.1, //mutationProbability,
                10000,//numIterations,
                16,//populationSize,
                "MOEAD"
        );
    }

    private static void rNSGAIII(GraphTable graph, Long[] points, long seed) {
        executeAlgorithm(
                new MOShortestPathProblem(graph, points[0], points[1], 10),
                seed,
                0.9, //crossoverProbability,
                0.1, //mutationProbability,
                10000,//numIterations,
                10,//populationSize,
                "NSGAIII"
        );
    }

    private static void executeAlgorithm(Problem problem, long seed, double crossoverProbability,
                                         double mutationProbability, int numIterations, int populationSize,
                                         String metaheuristic) {
        System.out.println("Problem: " + problem);
        System.out.println("Random seed: " + seed);
        System.out.println("Population size: " + populationSize);
        System.out.println("Number of iterations: " + numIterations);
        JMetalRandom.getInstance().setSeed(seed);

        //Algorithm<List<NodePathSolution>> algorithm = null;
        CrossoverOperator crossover;
        MutationOperator mutation;

        double crossoverDistributionIndex = 20.0 ;
        crossover = new Crossover1PLS(crossoverProbability, crossoverDistributionIndex);


        String label = "Colorado_" + seed + "_";
        Algorithm algorithm = null;

        if (metaheuristic.equals("MOEAD")) {
            double mutationDistributionIndex = 20.0 ;
            mutation = new Mutation1PChangeD((MOShortestPathProblemDouble) problem, mutationProbability, mutationDistributionIndex);
            label += "MOEAD";

            algorithm = new MOEADSTM(problem, populationSize, populationSize, numIterations, mutation, new MyDifferentialEvolutionCrossover(0.5, 0.5, "rand/1/bin"), MOEAD.FunctionType.TCHE, "es/uma/lcc/neo/cintrano/robustness/mo/shortestpath", 0.9, 5, 5);//selection, evaluator);

            InputStream in = algorithm.getClass().getResourceAsStream("/" + "" + "/" + "/W4D_10.dat");
            System.out.println(in);
            System.out.println(algorithm.getClass().getName());


            AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm).execute();
            List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDoubleSolution> population =
                    (List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDoubleSolution>) algorithm.getResult();
            long computingTime = algorithmRunner.getComputingTime();

            JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");

            fillLog(population, 111L);
            printStaticFinalLog(algorithm, computingTime, label);

            printFinalSolutionSet(population);

        }

        if (metaheuristic.equals("NSGAII")) {
            double mutationDistributionIndex = 20.0;
            mutation = new Mutation1PChange((MOShortestPathProblem) problem, mutationProbability, mutationDistributionIndex);
            label += "NSGAII";
            final BinaryTournamentSelection<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution> selection = new BinaryTournamentSelection<>(new RankingAndCrowdingDistanceComparator<>());
            SolutionListEvaluator<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution> evaluator = new SequentialSolutionListEvaluator<>();

            algorithm = new NSGAII<>(problem, numIterations, populationSize, crossover, mutation, selection, new DominanceComparator<>(), evaluator);
//
//            final MeasureManager measures = ((NSGAIIMeasures) algorithm).getMeasureManager();
//            PushMeasure<Object> pushMeasure = measures.getPushMeasure("currentIteration");
//            final Algorithm<List<NodePathSolution>> finalAlgorithm = algorithm;
//            //final Algorithm finalAlgorithm = algorithm;
//            System.err.print("-------> : ");
//            System.err.println((List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>) ((NSGAIIMeasures) finalAlgorithm).getPopulation());
//            System.err.println(pushMeasure);
//            pushMeasure.register(o -> {
//                fillLogN((List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>) ((NSGAIIMeasures) finalAlgorithm).getPopulation(), o);
//            });


            AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm).execute();
            List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution> population =
                    (List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>) algorithm.getResult();
            long computingTime = algorithmRunner.getComputingTime();

            JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");

            fillLogN(population, 111L);
            printStaticFinalLog(algorithm, computingTime, label);

            printFinalSolutionSet(population);
        }
//        if (algorithm != null) {
//            AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm).execute();
//            List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDoubleSolution> population =
//                    (List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDoubleSolution>) algorithm.getResult();
//            long computingTime = algorithmRunner.getComputingTime();
//
//            JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");
//
//            fillLog(population, 111L);
//            printStaticFinalLog(algorithm, computingTime, label);
//
//            printFinalSolutionSet(population);
//        }
    }
    private static void fillLog(List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDoubleSolution> population, Object o) {
        List<double[]> f = new ArrayList<>();
        List<Double[]> v = new ArrayList<>();
        for (es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDoubleSolution solution : population) {
            f.add(solution.getObjectives());
            v.add(solution.getVariables());
        }
        fitnessAll.put((Long) o, f);
        solutions.put((Long) o, v);
    }

    private static void fillLogN(List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution> population, Object o) {
        List<double[]> f = new ArrayList<>();
        List<Double[]> v = new ArrayList<>();
        for (es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution solution : population) {
            //lista.add(solution.getObjective(0));
            f.add(solution.getObjectives());

            v.add(solution.getVariablesDouble());
        }
        fitnessAll.put((Long) o, f);
        solutions.put((Long) o, v);
    }

    static void printStaticFinalLog(Algorithm<List<NodePathSolution>> algorithm, long time, String label) {
        Writer writer = null;
        String filename = "result_F.csv";
        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filename), StandardCharsets.UTF_8));
            int i;
            for (Long iteration: fitnessAll.keySet()){
                i = 0;
                for (double[] listA : fitnessAll.get(iteration)) {
                    writer.write(iteration + " " + i + " ");
                    for (double e : listA) {
                        writer.write(e + " ");
                    }
                    i++;
                    writer.write("\n");
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {
                assert writer != null;
                writer.close();} catch (Exception ex) {ex.printStackTrace();}
        }
        filename = "result_V.csv";
        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filename), StandardCharsets.UTF_8));
            int i;
            for (Long iteration: solutions.keySet()){
                i = 0;
                for (Double[] listA : solutions.get(iteration)) {
                    writer.write(iteration + " " + i + " ");
                    for (Double e : listA) {
                        writer.write(e + " ");
                    }
                    i++;
                    writer.write("\n");
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {writer.close();} catch (Exception ex) {ex.printStackTrace();}
        }


    }

}
