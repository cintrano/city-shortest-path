package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.astar.AstarMO;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDim;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.iterated.IteratedLS;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MOShortestPathProblem;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MOShortestPathProblemDouble;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDifferentialEvolutionCrossover;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.operators.Crossover1PLS;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.operators.Mutation1PChangeD;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.pulse.Pulse;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.pulse.PulseMO;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.utilities.ProcessGraph;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEAD;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEADSTM;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalLogger;
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
                    graph = ProcessGraph.prepareGraph("hbefa-malaga-graph.xml", "weights_time-hbefa.xml", "mapping-malaga.txt", "MAL");
                    ProcessGraph.fixVertexIndex(graph);
                    ProcessGraph.readTlLogics(graph, "sumo_processed.json");
                    break;
                case "Colorado":
                    graph = ProcessGraph.prepareGraph("hbefa-col-graph.xml", "col_weights_time-hbefa-MOD.xml", "mapping-col.txt", "COL");
                    break;
                case "NY":
                    graph = ProcessGraph.prepareGraph("hbefa-ny-graph.xml", "ny_weights_time-hbefa.xml", "mapping-ny.txt", "NY");
                    break;
            }

            // Points
            System.out.print("Reading points...");
            Long[] points;
            if (!args[2].equals("default")) {
                points = readPoints(Integer.parseInt(args[2]), args[0]);
            } else {
                points = new Long[]{8L, 720L};
            }

            // Algorithm
            if (points != null) {
                System.out.println(points[0] + " " + points[1]);
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
                        rNSGAII(graph, points, seed[Integer.parseInt(args[3])]);
                        break;
                    case "MOEAD":
                        System.out.println("Run MOEA-D...");
                        rMOEAD(graph, points, seed[Integer.parseInt(args[3])]);
                        break;
                    case "NSGAIII":
                        System.out.println("Run NSGA-III...");
                        rNSGAIII(graph, points, seed[Integer.parseInt(args[3])]);
                        break;
                    case "Iterated":
                        System.out.println("Run Iterated...");
                        rIterated(graph, points, seed[Integer.parseInt(args[3])]);
                        break;
                }
            }
        }

        System.out.println("=== END EXPERIMENTS ===");
    }




    private static Long[] readPoints(int id, String city) {
        Long[] points = null;
        BufferedReader br = null;
        try {
            switch (city) {
                case "Malaga":
                    br = new BufferedReader(new FileReader("input-MAL.data"));
                    break;
                case "Colorado":
                    br = new BufferedReader(new FileReader("input-COL.data"));
                    break;
                case "NY":
                    br = new BufferedReader(new FileReader("input-NY.data"));
                    break;
            }
            assert br != null;
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



    private static void rIterated(GraphTable graph, Long[] points, int seed) {
        long timeInit = System.currentTimeMillis();
        IteratedLS algorithm = new IteratedLS(seed);
        algorithm.setGraph(graph);
        long timeStart = System.currentTimeMillis();
        List<Node> path = algorithm.getPath(points[0], points[1], new float[]{1, 1, 1, 1});
        long timeEnd = System.currentTimeMillis();
        String wTag = "1-1-1-1";
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "Iterated", wTag);

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
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "Dijkstra", wTag);
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
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "A*", wTag);

    }

    private static void rPulse(GraphTable graph, Long[] randomPoints) {
        long timeInit = System.currentTimeMillis();
        PulseMO pulse = new PulseMO(objectives.length);
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

        System.out.println("NUMBER OF SOLUTIONS: " + solutions.size());
        System.out.println("... method end");

        System.out.println("Initialization time: \t" + (timeStart - timeInit));
        System.out.println("Execution time: \t" + (timeEnd - timeStart));
        System.out.println("... method end");

        System.out.println("................");
        System.out.println("\nPrint solutions FUN:");

        printSolutions(solutions, (timeStart - timeInit), (timeEnd - timeStart), "PULSE2", obj1 + "-" + obj2);
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
            if (out != null) try {
                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }


    private static Map<Long, List<double[]>> fitnessAll = new HashMap<>();
    private static Map<Long, List<Double[]>> solutions = new HashMap<>();
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

    private static void rMOEAD(GraphTable graph, Long[] points, long seed) {
        executeAlgorithm(
                new MOShortestPathProblemDouble(graph, points[0], points[1], 10),
                seed,
                0.9, //crossoverProbability,
                0.1, //mutationProbability,
                1000,//numIterations,
                10,//populationSize,
                "MOEAD"
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

        }


        if (algorithm != null) {
            AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm).execute();

            List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDoubleSolution> population = (List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2.MyDoubleSolution>) algorithm.getResult();
            long computingTime = algorithmRunner.getComputingTime();

            JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");

            fillLog(population, 111L);
            printStaticFinalLog(algorithm, computingTime, label);

            printFinalSolutionSet(population);
        }
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

    private static void printStaticFinalLog(Algorithm<List<NodePathSolution>> algorithm, long time, String label) {
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
}
