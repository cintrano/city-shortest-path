package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.astar.AstarMO;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDim;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.iterated.IteratedLS;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.iterated.MOIteratedLS;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.*;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.operators.Crossover1PLS;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.operators.Mutation1PChange;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.operators.Mutation1PChangeD;
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
import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIStoppingByTime;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;
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

    private final static int[] seed = new int[]{7, 3, 23, 29, 67, 79, 97, 101, 139, 199, 211, 269, 311, 347, 389,
            431, 461, 503, 569, 601, 607, 653, 691, 701, 733, 739, 761, 809, 811, 997};

    public final static Long[] objectives = new Long[]{0L, 1L, 2L, 3L};
    private final static float[] weights = new float[]{0.0001f, 0.25f, 0.5f, 0.75f, 1.0f};

    private final static int POPULATION_SIZE = 100;

    public static void main (String[] args) {
        System.out.println("=== START EXPERIMENTS ===");
        for (String s : args) {
            System.out.print(s + " ");
        }
        System.out.println();

        if (args[0].equals("Mapping")) {
            GraphTable graph = ProcessGraph.parserFile(args[1]);
            assert graph != null;
            ProcessGraph.printMapping(graph);
        } else if (args[0].equals("Reevaluate")) {
            String output_file = args[1];
            String path = args[2];

            GraphTable graph = ProcessGraph.prepareGraph("g_CC.xml", "w_filter.xml", "malaga.osm-mapping.txt", "MAL");
//            ProcessGraph.fixVertexIndex(graph);
            ProcessGraph.readTlLogics(graph, "malaga_tls.json");


            reevaluateSolutions(output_file, path, graph, args[3].equals("META"));
        } else if (args.length >= 3) {
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
                boolean robust = args[4].equals("Robust");
                boolean tlOption = args[5].equals("TL");
                int maxIterations = 10000;
                long maxTime = 0; // if default value want to be used must be 60*1000 = 1 min
                if (args.length > 6) {
                    maxIterations = Integer.parseInt(args[6]); // If it is 0, it will not use as stopping criteria
                    maxTime = Integer.parseInt(args[7]); // If it is 0, it will not use as stopping criteria
                }
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
                        executeAlgorithm(
                                new MOShortestPathProblem(graph, points[0], points[1], 10, robust, tlOption),
                                seed[Integer.parseInt(args[3])],
                                0.9, //crossoverProbability,
                                0.1, //mutationProbability,
                                maxIterations, //numIterations,
                                maxTime, //computationTime
                                POPULATION_SIZE, //populationSize,
                                "NSGAII"
                        );
                        break;
                    case "MOEAD":
                        System.out.println("Run MOEA-D...");
                        executeAlgorithm(
                                new MOShortestPathProblemDouble(graph, points[0], points[1], 10, robust, tlOption),
                                seed[Integer.parseInt(args[3])],
                                0.9, //crossoverProbability,
                                0.1, //mutationProbability,
                                maxIterations, //numIterations,
                                maxTime, //computationTime
                                POPULATION_SIZE, //populationSize,
                                "MOEAD"
                        );
                        break;
                    case "Iterated":
                        System.out.println("Run Iterated...");
                        rIterated(graph, points, seed[Integer.parseInt(args[3])]);
                        break;
                    case "MOIterated":
                        System.out.println("Run MOIterated...");
                        rMOIterated(graph, points, seed[Integer.parseInt(args[3])], robust, tlOption, maxIterations, maxTime);
                        break;
                }
            }
        } else {
            System.out.print("ERROR:: Incorrect parameters: " + Arrays.toString(args));
        }

        System.out.println("=== END EXPERIMENTS ===");
    }

    private static void reevaluateSolutions(String output_file, String path, GraphTable graph, boolean meta) {
        String FILE_NAME = "resultsVAR.ssv";
        if (meta) FILE_NAME = "result_V.csv";

        File file = new File(path);
        String[] directories = file.list((current, name) -> new File(current, name).isDirectory());
        System.out.println(Arrays.toString(directories));
        try {
            FileWriter fw = new FileWriter(output_file);
            assert directories != null;
            for (String directory : directories) {
                String[] tokens = directory.split("[.]");
                System.out.println(directory  + " " + tokens.length);
                BufferedReader br;
                float[] fitness;
                    br = new BufferedReader(new FileReader(path + directory + "/" + FILE_NAME));

                    String outLine;
                    String line = br.readLine();
                    int i = 0;
                    while (line != null) {
                        fitness = reevaluateSolution(line, graph, tokens[1].equals("Robust"), meta);
                        // TODO print line as output
                        outLine = "Reeval" + "\t";
                        outLine += tokens[0] + "\t" + (Integer.parseInt(tokens[3])%30) + "\t";
                        outLine += tokens[1] + "\t" + tokens[2] + "\t"; // Robust, TL
                        outLine += i + "\t" + tokens[0] + "\t"; // Dummy time

                        outLine += fitness[fitness.length - 1] + "\t";
                        outLine += fitness[0] + "\t" + fitness[1] + "\t" + fitness[2] + "\t" + fitness[3];
                        fw.write(outLine + "\n");
                        line = br.readLine();
                        i++;
                    }

                    br.close();
            }
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static float[] reevaluateSolution(String line, GraphTable graph, boolean robustFlag, boolean meta) {
        String[] tokens = line.split(" ");
        long[] sol;
        if (meta) {
            int start = 2;
            String[] solString = Arrays.copyOfRange(tokens, start, tokens.length);
            sol = new long[solString.length];
            for (int i = 0; i < sol.length; i++) {
                double d = Double.parseDouble(solString[i]);
                sol[i] = (new Double(d)).longValue();
            }
            Random rand = new Random();
            rand.setSeed(13);
            return MyUtility.fitnessTLMeta(graph, sol, 4, 30, robustFlag, rand);
        } else {
            int start = 7;
            String[] solString = Arrays.copyOfRange(tokens, start, tokens.length - 1);
            sol = Arrays.stream(solString).mapToLong(Long::parseLong).toArray();
            Random rand = new Random();
            rand.setSeed(13);
            return MyUtility.fitnessTL(graph, sol, 4, 30, robustFlag, rand);
        }
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

    private static void rMOIterated(GraphTable graph, Long[] points, int seed, boolean robust, boolean tlOption, int iterations, long time) {
        long timeInit = System.currentTimeMillis();
        MOIteratedLS algorithm = new MOIteratedLS(seed);
        algorithm.setGraph(graph);
        algorithm.setRobustFlag(robust);
        algorithm.setTLFlag(tlOption);
        algorithm.setMaxIterations(iterations);
        algorithm.setThresholdComputingTime((int) time);
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


//    private static Map<Long, List<double[]>> fitnessAll = new HashMap<>();
//    private static Map<Long, List<Double[]>> solutions = new HashMap<>();

    private static void executeAlgorithm(Problem problem, long seed, double crossoverProbability,
                                         double mutationProbability, int numIterations, long computationTime,
                                         int populationSize, String metaheuristic) {
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

            if (computationTime == 0) {
                algorithm = new MOEADSTM(problem, populationSize, populationSize, numIterations, mutation, new MyDifferentialEvolutionCrossover(0.5, 0.5, "rand/1/bin"), MOEAD.FunctionType.TCHE, "es/uma/lcc/neo/cintrano/robustness/mo/shortestpath", 0.9, 50, 50);//selection, evaluator);
                //algorithm = new MOEADSTM(problem, populationSize, populationSize, numIterations, mutation, new MyDifferentialEvolutionCrossover(0.5, 0.5, "rand/1/bin"), MOEAD.FunctionType.TCHE, "", 0.9, 50, 50);//selection, evaluator);
                InputStream in = algorithm.getClass().getResourceAsStream("/" + "" + "/" + "/W4D_100.dat");
                //System.out.println(in);
                //System.out.println(algorithm.getClass().getName());
            } else {
                algorithm = new MOEADSTMTime(problem, computationTime, populationSize, populationSize, numIterations, mutation, new MyDifferentialEvolutionCrossover(0.5, 0.5, "rand/1/bin"), MOEAD.FunctionType.TCHE, "es/uma/lcc/neo/cintrano/robustness/mo/shortestpath", 0.9, 50, 50);//selection, evaluator);
                InputStream in = algorithm.getClass().getResourceAsStream("/" + "" + "/" + "/W4D_100.dat");
                //System.out.println(in);
                //System.out.println(algorithm.getClass().getName());
            }


            AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm).execute();
            List<DoubleSolution> population = (List<DoubleSolution>) algorithm.getResult();
            long computingTime = algorithmRunner.getComputingTime();

            System.out.println("Final population size: " + population.size());
            JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");

//            fillLog(population, 111L);
            MyUtility.printStaticFinalLog(algorithm, computingTime);

            //printFinalSolutionSet(population);

        }

        if (metaheuristic.equals("NSGAII")) {
            double mutationDistributionIndex = 20.0;
            mutation = new Mutation1PChange((MOShortestPathProblem) problem, mutationProbability, mutationDistributionIndex);
            label += "NSGAII";
            final BinaryTournamentSelection<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.NodePathSolution> selection = new BinaryTournamentSelection<>(new RankingAndCrowdingDistanceComparator<>());
            SolutionListEvaluator<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.NodePathSolution> evaluator = new SequentialSolutionListEvaluator<>();

            if (computationTime == 0) {
                algorithm = new NSGAII<>(problem, numIterations, populationSize,
                        crossover, mutation, selection, new DominanceComparator<>(), evaluator);
            } else {
                algorithm = new MyNSGAIIStoppingByTime(problem, populationSize, computationTime,
                        crossover, mutation, selection, new DominanceComparator<>(), evaluator);
            }
            AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm).execute();
            List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.NodePathSolution> population =
                    (List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.NodePathSolution>) algorithm.getResult();
            long computingTime = algorithmRunner.getComputingTime();
            System.out.println("Final population size: " + population.size());

            JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");

//            fillLogN(population, 111L);
            MyUtility.printStaticFinalLog(algorithm, computingTime, label);

            //printFinalSolutionSet(population);
        }
    }
}
