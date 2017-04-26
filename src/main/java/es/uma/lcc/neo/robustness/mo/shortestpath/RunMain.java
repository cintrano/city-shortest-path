package es.uma.lcc.neo.robustness.mo.shortestpath;

import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.astar.Astar;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.astar.AstarMO;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDim;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWeighted;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.MOShortestPathProblem;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.operators.Crossover1PLS;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.operators.Mutation1PChange;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.pulse.Pulse;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.pulse.PulseMO;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.robustness.mo.shortestpath.utilities.ProcessGraph;
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.mocell.MOCellBuilder;
import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIMeasures;
import org.uma.jmetal.algorithm.multiobjective.spea2.SPEA2Builder;
import org.uma.jmetal.measure.MeasureListener;
import org.uma.jmetal.measure.MeasureManager;
import org.uma.jmetal.measure.PushMeasure;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.runner.multiobjective.NSGAIIIntegerRunner;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.comparator.RankingAndCrowdingDistanceComparator;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.io.*;
import java.lang.management.ManagementFactory;
import java.util.*;

import static org.uma.jmetal.runner.AbstractAlgorithmRunner.printFinalSolutionSet;

/**
 * Created by Christian Cintrano on 22/01/17.
 * Main Class
 */
public class RunMain {

    private final static Long[] points = new Long[]{
            352681962L,
            2397542792L,
            304830449L,
            189419939L,
            2334032985L,
            630737144L,
            1402836812L,
            4025969288L,
            2180436378L,
            913155755L,
            765640946L,
            729655979L,
            737470528L,
            2009776762L
    };

    private final static int[] seed = new int[]{0, 3, 23, 29, 67, 79, 97, 101, 139, 199, 211, 269, 311, 347, 389,
            431, 461, 503, 569, 601, 607, 653, 691, 701, 733, 739, 761, 809, 811, 997};


    public static void main (String[] args) throws IOException, InterruptedException {
        System.out.println("=== START EXPERIMENTS ===");
        for (String s : args) {
            System.out.print(s + " ");
        }
        System.out.println();

        if (args[0].equals("PrepareGraph")) {
            GraphTable graph = prepareGraph();
            ProcessGraph.printRandomWeights(graph, "wVar0.xml", 0L, 0.9f, 1.1f, 2);
            ProcessGraph.printRandomWeights(graph, "wVar1.xml", 1L, 0.9f, 1.1f, 3);
            ProcessGraph.printMapping(graph);
        }

        if (args.length >= 3) {
            GraphTable graph = null;
            // Map
            if (args[0].equals("Malaga")) {
                System.out.println("Loading graph...");
                graph = prepareGraph();
                graph = ProcessGraph.fixVertexIndex(graph);
            } else if (args[0].equals("Colorado")) {
                graph = prepareSimpleColoradoGraph();
            }

            // Points
            System.out.print("Reading points...");
            Long[] points = null;
            if (!args[2].equals("default")) {
                points = readPoints(args[2]);
            } else if (args[0].equals("Colorado") && args[2].equals("default")) {
                points = new Long[]{8L, 720L};
            } else if (args[0].equals("Colorado") && args[2].equals("default2")) {
                points = new Long[]{8L, 7324L};
            }

            // Algorithm
            if (points != null) {
                System.out.println(points[0] + " " + points[1]);
                if (Objects.equals(args[1], "PULSE4")) {
                    System.out.println("Run PULSE4...");
                    rPulse(graph, points);
                }
                if (Objects.equals(args[1], "PULSE2")) {
                    System.out.println("Run PULSE2...");
                    rPulse2Experiments(graph, points);
                }
                if (Objects.equals(args[1], "Dijkstra")) {
                    System.out.println("Run Dijkstra...");
                    rDijkstraExperiments(graph, points);
                }
                if (Objects.equals(args[1], "Astar")) {
                    System.out.println("Run A*...");
                    rAstarExperiments(graph, points);
                }
                if (Objects.equals(args[1], "NSGAII")) {
                    System.out.println("Run NSGA-II...");
                    rNSGAII(graph, points, seed[Integer.parseInt(args[3])]);
                }
                if (Objects.equals(args[1], "NSGAIII")) {
                    System.out.println("Run NSGA-III...");
                    rNSGAIII(graph, points, seed[Integer.parseInt(args[3])]);
                }
            }
        }

        if (args[0].equals("Google")) {
            GraphTable graph = prepareGraph();

            for (int aSeed : seed) {
                for (int i = 0; i < points.length - 1; i++) {
                    for (int j = i + 1; j < points.length; j++) {
                        runCommandGoogle(graph, aSeed, points[i], points[j]);
                        runCommandGoogle(graph, aSeed, points[j], points[i]);
                    }
                }
            }
        }
        if (args[0].equals("GoogleMO")) {
            GraphTable graph = prepareGraph();

            for (int aSeed : seed) {
                for (int i = 0; i < points.length - 1; i++) {
                    for (int j = i + 1; j < points.length; j++) {
                        runCommandGoogleMO(graph, aSeed, points[i], points[j]);
                        runCommandGoogleMO(graph, aSeed, points[j], points[i]);
                    }
                }
            }
        }
        if (args[0].equals("ColoradoN")) {
            System.out.println("Loading graph...");
            GraphTable graph = prepareColoradoGraph();
            //List<Long> puntos = getPoints(graph);
            long seed = Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0] + System.currentTimeMillis());//0;
            /*
            seed = Integer.parseInt(ManagementFactory.getRuntimeMXBean().getName().split("@")[0] );
            System.out.println("S " + seed);
            seed = Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0]) + System.currentTimeMillis();
            System.out.println("S " + seed);
            seed = Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0] + System.currentTimeMillis());
            System.out.println("S " + seed);
*/
            System.out.println("Start algorithm...");
/*
            es.uma.lcc.neo.cintrano.RobustShortestPathMain.executeAlgorithm(
                    new MOShortestPathProblem(graph, 8L, 720L, 10),
                    seed,
                    0.9,
                    0.1,
                    1000,
                    10,
                    "NSGAII"
            );
*/
        }
        if (args[0].equals("ColoradoD")) {
            System.out.println("Loading graph...");
            GraphTable graph = prepareColoradoGraph();
            //List<Long> puntos = getPoints(graph);
            long seed = Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0] + System.currentTimeMillis());//0;
            /*
            seed = Integer.parseInt(ManagementFactory.getRuntimeMXBean().getName().split("@")[0] );
            System.out.println("S " + seed);
            seed = Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0]) + System.currentTimeMillis();
            System.out.println("S " + seed);
            seed = Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0] + System.currentTimeMillis());
            System.out.println("S " + seed);
*/
            System.out.println("Start algorithm...");

        }
        if (args[0].equals("ColoradoAstar")) {
            System.out.println("Loading graph...");
            GraphTable graph = prepareColoradoGraph();
            //List<Long> puntos = getPoints(graph);
            long seed = Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0] + System.currentTimeMillis());//0;
            /*
            seed = Integer.parseInt(ManagementFactory.getRuntimeMXBean().getName().split("@")[0] );
            System.out.println("S " + seed);
            seed = Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0]) + System.currentTimeMillis();
            System.out.println("S " + seed);
            seed = Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0] + System.currentTimeMillis());
            System.out.println("S " + seed);
*/
            System.out.println("Start algorithm...");

            /*
            es.uma.lcc.neo.cintrano.RobustShortestPathMain.executeAlgorithm(
                    new MOShortestPathProblem(graph, 8L, 720L, 10),
                    seed,
                    0.9,
                    0.1,
                    1000,
                    10,
                    "NSGAII"
            );
*/
        }
        if (args[0].equals("ColoradoM")) {
            System.out.println("Loading graph...");
            GraphTable graph = prepareColoradoGraph();
            //List<Long> puntos = getPoints(graph);
            int seed = 0;

            System.out.println("Start algorithm...");
            /*
            es.uma.lcc.neo.cintrano.RobustShortestPathMain.executeAlgorithm(
                    new MOShortestPathProblem(graph, 8L, 720L, 10),
                    seed,
                    0.9,
                    0.1,
                    1000,
                    10,
                    "MOCell"
            );
            */
        }
        if (args[0].equals("ColoradoS")) {
            System.out.println("Loading graph...");
            GraphTable graph = prepareColoradoGraph();
            //List<Long> puntos = getPoints(graph);
            int seed = 0;

            System.out.println("Start algorithm...");
            /*
            es.uma.lcc.neo.cintrano.RobustShortestPathMain.executeAlgorithm(
                    new MOShortestPathProblem(graph, 8L, 720L, 10),
                    seed,
                    0.9,
                    0.1,
                    1000,
                    10,
                    "SPEA2"
            );
            */
        }

        System.out.println("=== END EXPERIMENTS ===");
    }


    private static Long[] readPoints(String idString) {
        int id = Integer.parseInt(idString);
        Long[] points = null;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader("input.data"));
            String line = br.readLine();

            int count = 0;
            while (line != null && count <= id) {
                if (count == id) {
                    String[] array = line.split(" ");
                    points = new Long[]{Long.parseLong(array[0]), Long.parseLong(array[1])};
                }
                count++;
                line = br.readLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return points;
    }




    private static void rDijkstraExperiments(GraphTable graph, Long[] randomPoints) {
        Long[] objectives = new Long[]{0L, 1L, 2L, 3L};
        float[] weights = new float[]{0.0001f, 0.25f, 0.75f, 1.0f};

        for (int w0 = 0; w0 < weights.length; w0++) {
            for (int w1 = 0; w1 < weights.length; w1++) {
                for (int w2 = 0; w2 < weights.length; w2++) {
                    for (int w3 = 0; w3 < weights.length; w3++) {
                        if (weights[w0] + weights[w1] + weights[w2] + weights[w3] > 0) {
                            rDijkstraWD(graph, randomPoints, objectives, new float[]{w0, w1, w2, w3});
                        }
                    }
                }
            }
        }
    }

    private static void rDijkstraWD(GraphTable graph, Long[] randomPoints, Long[] objectives, float[] weights) {
        long timeInit = System.currentTimeMillis();
        DijkstraWNDim dj = new DijkstraWNDim(objectives);
        dj.setWeights(weights);
        dj.setGraph(graph);
        long timeStart = System.currentTimeMillis();
        List<Node> path = dj.getPath(randomPoints[0], randomPoints[1]);
        long timeEnd = System.currentTimeMillis();
        //sol.add(new NodePathSolution(graph.getFitness(path, ""), path));
        String wTag = weights[0] + "-" +weights[1] + "-" +weights[2] + "-" +weights[3];
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "Dijkstra", wTag);
    }




    private static void rAstarExperiments(GraphTable graph, Long[] randomPoints) {
        float[] weights = new float[]{0.0001f, 0.25f, 0.75f, 1.0f};

        for (int w0 = 0; w0 < weights.length; w0++) {
            for (int w1 = 0; w1 < weights.length; w1++) {
                for (int w2 = 0; w2 < weights.length; w2++) {
                    for (int w3 = 0; w3 < weights.length; w3++) {
                        if (weights[w0] + weights[w1] + weights[w2] + weights[w3] > 0) {
                            rAstarMO(graph, randomPoints, new float[]{w0, w1, w2, w3});
                        }
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
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "A*", wTag);

    }

    private static void rAstar(GraphTable graph, Long[] randomPoints) {
        long timeInit = System.currentTimeMillis();
        Astar astar = new Astar();
        astar.setGraph(graph);
        astar.setTarget(1F);
        long timeStart = System.currentTimeMillis();
        List<Node> path = astar.getPath(randomPoints[0], randomPoints[1]);
        long timeEnd = System.currentTimeMillis();
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "A*", "0");

        timeInit = System.currentTimeMillis();
        astar = new Astar();
        astar.setGraph(graph);
        astar.setTarget(0F);
        timeStart = System.currentTimeMillis();
        path = astar.getPath(randomPoints[0], randomPoints[1]);
        timeEnd = System.currentTimeMillis();
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "A*", "1");

        timeInit = System.currentTimeMillis();
        astar = new Astar();
        astar.setGraph(graph);
        astar.setTarget(0.5F);
        timeStart = System.currentTimeMillis();
        path = astar.getPath(randomPoints[0], randomPoints[1]);
        timeEnd = System.currentTimeMillis();
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "A*", "2");
    }

    private static Long[] selectRandomPoints(GraphTable graph, Random random) {
        List<Long> vertex = new ArrayList<>(graph.getIntersections().keySet());
        System.out.println("----- SIZE " + graph.getAdjacencyMatrix().rowKeySet().size() + " " + vertex.size());
        Long[] index = new Long[2];
        do {
            index[0] = vertex.get(random.nextInt(vertex.size()));
            index[1] = vertex.get(random.nextInt(vertex.size()));
        } while (index[0].equals(index[1]));
        return index;
    }

    private static void rPulse(GraphTable graph, Long[] randomPoints) {
        long timeInit = System.currentTimeMillis();
        PulseMO pulse = new PulseMO(4);
        pulse.setGraph(graph);

        long timeStart = System.currentTimeMillis();
        Set<NodePathSolution> solutions = pulse.pulseAlgorithm(randomPoints[0], randomPoints[1]);
        long timeEnd = System.currentTimeMillis();

        System.out.println("NUMBER OF SOLUTIONS: " + solutions.size());
        System.out.println("... method end");

        System.out.println("Initialization time: \t" + (timeStart - timeInit));
        System.out.println("Execution time: \t" + (timeEnd - timeStart));
        System.out.println("... method end");

        System.out.println("................");
        System.out.println("\nPrint solutions FUN:");

        printSolutions(solutions, (timeStart - timeInit), (timeEnd - timeStart), "PULSE4");
        /*for (NodePathSolution s : solutions) {
            for (float l : s.getObjectives()) {
                System.out.print(l + " ");
            }
            System.out.println();
        }
        System.out.println("\nPrint solutions VAR:");
        for (NodePathSolution s : solutions) {
            for (Long l : s.getVariables()) {
                System.out.print(l + " ");
            }
            System.out.println();
        }
        */
    }

    private static void rPulse2Experiments(GraphTable graph, Long[] randomPoints) {
        Long[] objectives = new Long[]{0L, 1L, 2L, 3L};
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

        System.out.println("NUMBER OF SOLUTIONS: " + solutions.size());
        System.out.println("... method end");

        System.out.println("Initialization time: \t" + (timeStart - timeInit));
        System.out.println("Execution time: \t" + (timeEnd - timeStart));
        System.out.println("... method end");

        System.out.println("................");
        System.out.println("\nPrint solutions FUN:");

        printSolutions(solutions, (timeStart - timeInit), (timeEnd - timeStart), "PULSE2", obj1 + "-" + obj2);
        /*for (NodePathSolution s : solutions) {
            for (float l : s.getObjectives()) {
                System.out.print(l + " ");
            }
            System.out.println();
        }
        System.out.println("\nPrint solutions VAR:");
        for (NodePathSolution s : solutions) {
            for (Long l : s.getVariables()) {
                System.out.print(l + " ");
            }
            System.out.println();
        }
        */
    }

    private static void printSolutions(GraphTable graph, NodePathSolution s, long t1, long t2, String algorithm, String i) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("resultsFUN.ssv", true));
                String line = graph.getMapping().get(s.getVariables()[0]) + " ";
                line += graph.getMapping().get(s.getVariables()[s.getVariables().length-1]) + " ";
                line += t1 + " ";
                line += t2 + " ";
                line += algorithm + " ";
                line += i + " ";
                line += s.getObjectives()[0] + " ";
                line += s.getObjectives()[1] + " ";
                line += s.getObjectives()[2] + " ";
                line += s.getObjectives()[3] + " ";
                line += "\n";
                out.write(line);

        } catch (IOException e) {
            // error processing code
        } finally {
            if (out != null) {
                try {
                    out.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }


    public static Map<Long, List<double[]>> fitnessAll = new HashMap<Long, List<double[]>>();
    public static Map<Long, List<Long[]>> solutions = new HashMap<Long, List<Long[]>>();
    private static void rNSGAII(GraphTable graph, Long[] points, long seed) {
        executeAlgorithm(
                new MOShortestPathProblem(graph, points[0], points[1], 10),
                seed,
                0.9, //crossoverProbability,
                0.1, //mutationProbability,
                1000,//numIterations,
                10,//populationSize,
                "NSGAII"
        );
    }

    private static void rNSGAIII(GraphTable graph, Long[] points, long seed) {
        executeAlgorithm(
                new MOShortestPathProblem(graph, points[0], points[1], 10),
                seed,
                0.9, //crossoverProbability,
                0.1, //mutationProbability,
                1000,//numIterations,
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

        double mutationDistributionIndex = 20.0 ;
        mutation = new Mutation1PChange((MOShortestPathProblem) problem, mutationProbability, mutationDistributionIndex);

        String label = "Colorado_" + seed + "_";
        Algorithm algorithm = null;
        if (metaheuristic.equals("NSGAII")) {
            label += "NSGAII";
            /*
            algorithm = new NSGAIIBuilder<NodePathSolution>(problem, crossover, mutation)
                    .setMaxIterations(numIterations)
                    .setPopulationSize(populationSize)
                    .build();
                    */
            final BinaryTournamentSelection<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution> selection = new BinaryTournamentSelection<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>(new RankingAndCrowdingDistanceComparator<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>());
            SolutionListEvaluator<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution> evaluator = new SequentialSolutionListEvaluator<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>();

            algorithm = new NSGAIIMeasures<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>(problem, numIterations, populationSize, crossover, mutation, selection, evaluator);
            final MeasureManager measures = ((NSGAIIMeasures) algorithm).getMeasureManager();
            PushMeasure<Object> pushMeasure = measures.getPushMeasure("currentIteration");
            //final Algorithm<List<NodePathSolution>> finalAlgorithm = algorithm;
            final Algorithm finalAlgorithm = algorithm;
            pushMeasure.register(new MeasureListener<Object>() {
                public void measureGenerated(Object o) {
                    fillLog((List<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>) ((NSGAIIMeasures) finalAlgorithm).getPopulation(), o);
                }
            });
        }
        if (metaheuristic.equals("MOCell")) {
            label += "MOCell";
            algorithm = new MOCellBuilder<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>(problem, crossover, mutation)
                    .setMaxEvaluations(numIterations)
                    .setPopulationSize(populationSize)
                    .build();
        }
        if (metaheuristic.equals("SPEA2")) {
            label += "SPEA2";

            algorithm = new SPEA2Builder<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>(problem, crossover, mutation)
                    .setMaxIterations(numIterations)
                    .setPopulationSize(populationSize)
                    .build();
        }

        if (algorithm != null) {
            AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm).execute();

            List<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution> population = (List<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution>) algorithm.getResult();
            long computingTime = algorithmRunner.getComputingTime();

            JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");

            fillLog(population, new Long(111));
            /*
            for (NodePathSolution solution : population) {
                System.out.println(solution);
            }
            */
            printStaticFinalLog(algorithm, computingTime, label);
//            printFinalLog(algorithm, computingTime, label);

            printFinalSolutionSet(population);
        }
    }
    private static void fillLog(List<es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution> population, Object o) {
        List<double[]> f = new ArrayList<double[]>();
        List<Long[]> v = new ArrayList<Long[]>();
        for (es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2.NodePathSolution solution : population) {
            //lista.add(solution.getObjective(0));
            f.add(solution.getObjectives());
            v.add(solution.getVariables());
        }
        fitnessAll.put((Long) o, f);
        solutions.put((Long) o, v);
    }


    public static void printStaticFinalLog(Algorithm<List<NodePathSolution>> algorithm, long time, String label) {
        Writer writer = null;
        //String filename = label + "_log_" + algorithm.getResult().get(0).getVariables()[0] + "_" + algorithm.getResult().get(0).getVariables()[algorithm.getResult().get(0).getVariables().length - 1] + "_" + JMetalRandom.getInstance().getSeed() + "_F" + ".csv";
        String filename = "result_F.csv";
        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filename), "utf-8"));
            int i;
            for (Long iter: fitnessAll.keySet()){
                i = 0;
                for (double[] lista : fitnessAll.get(iter)) {
                    writer.write(iter + " " + i + " ");
                    for (double e : lista) {
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


        //filename = label + "_log_" + algorithm.getResult().get(0).getVariables()[0] + "_" + algorithm.getResult().get(0).getVariables()[algorithm.getResult().get(0).getVariables().length - 1] + "_" + JMetalRandom.getInstance().getSeed() + "_V" + ".csv";
        filename = "result_V.csv";
        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filename), "utf-8"));
            int i;
            for (Long iter: solutions.keySet()){
                i = 0;
                for (Long[] lista : solutions.get(iter)) {
                    writer.write(iter + " " + i + " ");
                    for (Long e : lista) {
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



    private static void printSolutions(Set<NodePathSolution> solutions, long t1, long t2, String algorithm) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("resultsFUN.ssv", true));
            int i = 0;
            String line;
            for (NodePathSolution s : solutions) {
                line = s.getVariables()[0] + " ";
                line += s.getVariables()[s.getVariables().length-1] + " ";
                line += t1 + " ";
                line += t2 + " ";
                line += algorithm + " ";
                line += i + " ";
                line += s.getObjectives()[0] + " ";
                line += s.getObjectives()[1] + " ";
                line += s.getObjectives()[2] + " ";
                line += s.getObjectives()[3] + " ";
                line += "\n";
                out.write(line);
                i++;
            }
        } catch (IOException e) {
            // error processing code
        } finally {
            if (out != null) {
                try {
                    out.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private static void printSolutions(Set<NodePathSolution> solutions, long t1, long t2, String algorithm, String tag) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("resultsFUN.ssv", true));
            int i = 0;
            String line;
            for (NodePathSolution s : solutions) {
                line = s.getVariables()[0] + " ";
                line += s.getVariables()[s.getVariables().length-1] + " ";
                line += t1 + " ";
                line += t2 + " ";
                line += algorithm + " ";
                line += i + "-" + tag + " ";
                line += s.getObjectives()[0] + " ";
                line += s.getObjectives()[1] + " ";
                line += s.getObjectives()[2] + " ";
                line += s.getObjectives()[3] + " ";
                line += "\n";
                out.write(line);
                i++;
            }
        } catch (IOException e) {
            // error processing code
        } finally {
            if (out != null) {
                try {
                    out.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private static GraphTable prepareColoradoGraph() {
        GraphTable graph = ProcessGraph.parserFile("USA-road-d.COL.co");
        //System.out.println("Adding weight to the graph...");
        graph = ProcessGraph.applyArcs(graph, 4L, "USA-road-d.COL.gr");
        //System.out.println("Adding weight to the graph...");
        graph = ProcessGraph.applyArcs(graph, 0L, "USA-road-t.COL.gr");
        graph = ProcessGraph.computeNewWeight(graph, 1L, 5L, 4L, 0L);

        //ProcessGraph.printSCCs(graph);

        System.out.println("Normalizing...");
        graph = ProcessGraph.normalize(graph);
        //ProcessGraph.printGraph(graph, "gCOL.xml");
        //ProcessGraph.printWeights(graph, "wCOL.xml");

        System.out.println("Loading graph...");
        //graph = ProcessGraph.parserGraph("gCOL.xml");
        System.out.println("Adding weight to the graph...");
        //graph = ProcessGraph.applyWeights(graph, "wCOL.xml");
        //ProcessGraph.printRandomWeights(graph, "wCOL2.xml", 0, 10, 2);
        //ProcessGraph.printRandomWeights(graph, "wCOL4.xml", 0, 15, 2);
        graph = ProcessGraph.applyWeights(graph, "wCOL2.xml");
        graph = ProcessGraph.applyWeights(graph, "wCOL4.xml");
        return graph;
    }

    private static GraphTable prepareSimpleColoradoGraph() {
        GraphTable graph = ProcessGraph.parserFile("USA-road-d.COL.co");
        graph = ProcessGraph.applyArcs(graph, 4L, "USA-road-d.COL.gr");
        graph = ProcessGraph.applyArcs(graph, 0L, "USA-road-t.COL.gr");
        graph = ProcessGraph.computeNewWeight(graph, 1L, 5L, 4L, 0L);

    //    graph = ProcessGraph.normalizate(graph);

        graph = ProcessGraph.applyWeights(graph, "wCOL2.xml");
        graph = ProcessGraph.applyWeights(graph, "wCOL4.xml");

        // Remove auxiliary columns
        graph.getWeightsMatrix().column(4L).clear();
        graph.getWeightsMatrix().column(5L).clear();
        return graph;
    }

    private static GraphTable prepareGraph() {
        // Graph
        String graphFilePath = "new-malaga-graph.xml";//"graph_connected.xml";
        String weightFilePath0 = "weights_time-noise.xml";//"wNew.xml";
        GraphTable graph = ProcessGraph.parserFile(graphFilePath);
        graph = ProcessGraph.applyWeights(graph, weightFilePath0);
        //graph = ProcessGraph.applyWeights(graph, "wVar0.xml");
        //graph = ProcessGraph.applyWeights(graph, "wVar1.xml");
        graph.getWeightsMatrix().column(10L).clear();
        graph = ProcessGraph.applyMapping(graph, "mapping-malaga.txt");
        return graph;
    }


    private static void runCommandGoogle(GraphTable graph, int seed, Long start, Long end) {
        // Algorithm
        double crossoverProbability = 0.9d;
        double mutationProbability = 0.1d;
        int numIterations = 100;
        int populationSize = 10;

        System.out.println("Start algorithm...");
        /*
        executeAlgorithm(
                new SOShortestPathProblem(graph, start, end, 10),
                seed,
                crossoverProbability,
                mutationProbability,
                numIterations,
                populationSize
        );
        */
    }

    private static void runCommandGoogleMO(GraphTable graph, int seed, Long start, Long end) {
        // Algorithm
        double crossoverProbability = 0.9d;
        double mutationProbability = 0.1d;
        int numIterations = 10;
        int populationSize = 10;

        System.out.println("Start algorithm...");
        /*
        es.uma.lcc.neo.cintrano.RobustShortestPathMain.executeAlgorithm(
                new MOShortestPathProblem(graph, start, end, 10),
                seed,
                crossoverProbability,
                mutationProbability,
                numIterations,
                populationSize,
                "NSGAII"
        );
        */
    }

    private static void runCommandGoogle(int seed, Long p1, Long p2) {
        String cmd = "java -cp SixPack-1.0-SNAPSHOT-jar-with-dependencies.jar es.uma.lcc.neo.cintrano.GoogleShortestPathMain " +
            "graph_connected.xml wNew.xml ";
        cmd += p1 + " " + p2 + " " + seed + " ";
        cmd += "0.9 0.1 100 10";
        System.out.println("$$ running command:: " + cmd);
        final Process p;
        try {
            p = Runtime.getRuntime().exec(cmd);

            new Thread(new Runnable() {
                public void run() {
                    BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String line = null;

                    try {
                        while ((line = input.readLine()) != null) {
                            System.out.println(line);
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }).start();

            p.waitFor();
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }
}
