package es.uma.lcc.neo.robustness.mo.shortestpath;

import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.astar.Astar;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWeighted;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.pulse.PulseMO;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.robustness.mo.shortestpath.utilities.ProcessGraph;

import java.io.*;
import java.lang.management.ManagementFactory;
import java.util.*;

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

        if (args[0].equals("PrepareGraph")) {
            GraphTable graph = prepareGraph();
            //ProcessGraph.printRandomWeights(graph, "wVar0.xml", 0L, 0.9f, 1.1f, 2);
            //ProcessGraph.printRandomWeights(graph, "wVar1.xml", 1L, 0.9f, 1.1f, 3);
//            graph = fixVertexIndex(graph);
            ProcessGraph.printMapping(graph);
        }
        if (args.length == 3) {
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
            }

            // Algorithm
            if (points != null) {
                System.out.println(points[0] + " " + points[1]);
                if (Objects.equals(args[1], "PULSE")) {
                    System.out.println("Run PULSE...");
                    rPulse(graph, points);
                }
                if (Objects.equals(args[1], "Dijkstra")) {
                    System.out.println("Run Dijkstra...");
                    rDijkstra(graph, points);
                }
                if (Objects.equals(args[1], "Astar")) {
                    System.out.println("Run A*...");
                    rAstar(graph, points);
                }
            }
        }

        if (args[0].equals("Normal")) {
            System.out.println("................");
            System.out.println("Loading graph...");
            GraphTable graph = prepareSimpleColoradoGraph();

            System.out.println("................");
            System.out.println("................");
            System.out.println("Running method...");
//            Set<NodePathSolution> solutions = NormalMethod.compute(graph,
//                    4, graph.getIntersections().get(8L), graph.getIntersections().get(720L));
            //Set<NodePathSolution> solutions = NormalMethod.compute(graph,
            //        4, graph.getIntersections().get(8L), graph.getIntersections().get(7324L));
            System.out.println("... method end");
            System.out.println("................");
            System.out.println("\nPrint solutions:");
//            for (NodePathSolution s : solutions) {
//                System.out.println(s);
//            }

        }
        if (args[0].equals("Google")) {
            GraphTable graph = prepareGraph();

            for (int s = 0; s < seed.length; s++) {
                for (int i = 0; i < points.length - 1; i++) {
                    for (int j = i + 1; j < points.length; j++) {
                        runCommandGoogle(graph, seed[s], points[i], points[j]);
                        runCommandGoogle(graph, seed[s], points[j], points[i]);
                    }
                }
            }
        }
        if (args[0].equals("GoogleMO")) {
            GraphTable graph = prepareGraph();

            for (int s = 0; s < seed.length; s++) {
                for (int i = 0; i < points.length - 1; i++) {
                    for (int j = i + 1; j < points.length; j++) {
                        runCommandGoogleMO(graph, seed[s], points[i], points[j]);
                        runCommandGoogleMO(graph, seed[s], points[j], points[i]);
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
/*
            System.out.println("A star");
            Astar astar = new Astar();
            astar.setGraph(graph);
            System.out.println("f01");
            astar.setTarget(0F);
            Float sum0 = 0F;
            Float sum1 = 0F;
            List<Node> path = astar.getPath(graph.getIntersections().get(8L), graph.getIntersections().get(720L));
            for (int i = 0; i < path.size() - 1; i++) {
                sum0 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 0L);
                sum1 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 1L);
            }
            System.out.println(sum0 + " " + sum1);

            /*

            Astar astar = new Astar();
            astar.setGraph(graph);
            System.out.println("f" + 0);
            astar.setTarget(0L);
            Float sum0 = 0F;
            Float sum1 = 0F;
            LinkedList<IntersectionInterface> path = astar.getPath(graph.getIntersections().get(8L), graph.getIntersections().get(720L));
            for (int i = 0; i < path.size() - 1; i++) {
                sum0 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 0L);
                sum1 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 1L);
            }
            System.out.println(sum0 + " " + sum1);
            System.out.println("f" + 1);
            astar.setTarget(1L);
            sum0 = 0F;
            sum1 = 0F;
            path = astar.getPath(graph.getIntersections().get(8L), graph.getIntersections().get(720L));
            for (int i = 0; i < path.size() - 1; i++) {
                sum0 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 0L);
                sum1 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 1L);
            }
            System.out.println(sum0 + " " + sum1);
*/
/*
            System.out.println("Dj");
            DijkstraWeighted dj = new DijkstraWeighted(0L, 1L, 0.5f);
            dj.setGraph(graph);
            System.out.println("f01");
            sum0 = 0F;
            sum1 = 0F;
            path = dj.getPath(graph.getIntersections().get(8L), graph.getIntersections().get(720L));
            for (int i = 0; i < path.size() - 1; i++) {
                sum0 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 0L);
                sum1 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 1L);
            }
            System.out.println(sum0 + " " + sum1);


            System.out.println("Dj");
            dj = new DijkstraWeighted(0L, 1L, 1);
            dj.setGraph(graph);
            System.out.println("f" + 0);
            sum0 = 0F;
            sum1 = 0F;
            path = dj.getPath(graph.getIntersections().get(8L), graph.getIntersections().get(720L));
            for (int i = 0; i < path.size() - 1; i++) {
                sum0 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 0L);
                sum1 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 1L);
            }
            System.out.println(sum0 + " " + sum1);
            System.out.println("f" + 1);
            dj = new DijkstraWeighted(0L, 1L, 0);
            dj.setGraph(graph);
            sum0 = 0F;
            sum1 = 0F;
            path = dj.getPath(graph.getIntersections().get(8L), graph.getIntersections().get(720L));
            for (int i = 0; i < path.size() - 1; i++) {
                sum0 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 0L);
                sum1 += graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(path.get(i).getId(), path.get(i+1).getId()), 1L);
            }
            System.out.println(sum0 + " " + sum1);

            System.out.println("-------");
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







    private static void rDijkstra(GraphTable graph, Long[] randomPoints) {
        long timeInit = System.currentTimeMillis();
        DijkstraWeighted dj = new DijkstraWeighted(0L, 1L, 1f);
        dj.setGraph(graph);
        long timeStart = System.currentTimeMillis();
        List<Node> path = dj.getPath(randomPoints[0], randomPoints[1]);
        long timeEnd = System.currentTimeMillis();
        //sol.add(new NodePathSolution(graph.getFitness(path, ""), path));
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "Dijkstra", 0);

        timeInit = System.currentTimeMillis();
        dj = new DijkstraWeighted(0L, 1L, 0f);
        dj.setGraph(graph);
        timeStart = System.currentTimeMillis();
        path = dj.getPath(randomPoints[0], randomPoints[1]);
        //sol.add(new NodePathSolution(graph.getFitness(path, ""), path));
        timeEnd = System.currentTimeMillis();
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "Dijkstra", 1);

        timeInit = System.currentTimeMillis();
        dj = new DijkstraWeighted(0L, 1L, 0.5f);
        dj.setGraph(graph);
        timeStart = System.currentTimeMillis();
        path = dj.getPath(randomPoints[0], randomPoints[1]);
        timeEnd = System.currentTimeMillis();
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "Dijkstra", 2);
    }

    private static void rAstar(GraphTable graph, Long[] randomPoints) {
        long timeInit = System.currentTimeMillis();
        Astar astar = new Astar();
        astar.setGraph(graph);
        astar.setTarget(1F);
        long timeStart = System.currentTimeMillis();
        List<Node> path = astar.getPath(randomPoints[0], randomPoints[1]);
        long timeEnd = System.currentTimeMillis();
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "A*", 0);

        timeInit = System.currentTimeMillis();
        astar = new Astar();
        astar.setGraph(graph);
        astar.setTarget(0F);
        timeStart = System.currentTimeMillis();
        path = astar.getPath(randomPoints[0], randomPoints[1]);
        timeEnd = System.currentTimeMillis();
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "A*", 1);

        timeInit = System.currentTimeMillis();
        astar = new Astar();
        astar.setGraph(graph);
        astar.setTarget(0.5F);
        timeStart = System.currentTimeMillis();
        path = astar.getPath(randomPoints[0], randomPoints[1]);
        timeEnd = System.currentTimeMillis();
        printSolutions(graph, new NodePathSolution(graph.getFitness(path, ""), path), timeStart - timeInit, timeEnd - timeStart, "A*", 2);
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

        printSolutions(solutions, (timeStart - timeInit), (timeEnd - timeStart), "PULSE");
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

    private static void printSolutions(GraphTable graph, NodePathSolution s, long t1, long t2, String algorithm, int i) {
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
        String graphFilePath = "graph_connected.xml";
        String weightFilePath0 = "wNew.xml";
        GraphTable graph = ProcessGraph.parserFile(graphFilePath);
        System.out.println("Adding weight to the graph...");
        graph = ProcessGraph.applyWeights(graph, weightFilePath0);
        graph = ProcessGraph.applyWeights(graph, "wVar0.xml");
        graph = ProcessGraph.applyWeights(graph, "wVar1.xml");
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
