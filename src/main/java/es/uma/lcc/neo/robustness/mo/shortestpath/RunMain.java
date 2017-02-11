package es.uma.lcc.neo.robustness.mo.shortestpath;

import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.Astar;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.DijkstraWeighted;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.Pulse;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.robustness.mo.shortestpath.utilities.ProcessGraph;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.util.List;
import java.util.Set;

/**
 * Created by cintrano on 22/01/17.
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

        if (args[0].equals("Pulse")) {
            System.out.println("................");
            System.out.println("Loading graph...");
            GraphTable graph = prepareSimpleColoradoGraph();

            System.out.println("................");
            System.out.println("................");
            System.out.println("Running method...");

            Pulse pulse = new Pulse();
            pulse.setGraph(graph);
            Set<List<Long>> solutions = pulse.pulseAlgorithm(8L, 720L);

            System.out.println("NUMBER OF SOLUTIONS: " + solutions.size());
            System.out.println("... method end");
            System.out.println("................");
            System.out.println("\nPrint solutions:");
            for (List<Long> s : solutions) {
                for (Long l : s) {
                    System.out.print(l + " ");
                }
                System.out.println();
            }

        }

        if (args[0].equals("Normal")) {
            System.out.println("................");
            System.out.println("Loading graph...");
            GraphTable graph = prepareSimpleColoradoGraph();

            System.out.println("................");
            System.out.println("................");
            System.out.println("Running method...");
            Set<NodePathSolution> solutions = NormalMethod.compute(graph,
                    4, graph.getIntersections().get(8L), graph.getIntersections().get(720L));
            //Set<NodePathSolution> solutions = NormalMethod.compute(graph,
            //        4, graph.getIntersections().get(8L), graph.getIntersections().get(7324L));
            System.out.println("... method end");
            System.out.println("................");
            System.out.println("\nPrint solutions:");
            for (NodePathSolution s : solutions) {
                System.out.println(s);
            }

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

            System.out.println("A star");
            Astar astar = new Astar();
            astar.setGraph(graph);
            System.out.println("f01");
            astar.setTarget(0L);
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

    private static GraphTable prepareColoradoGraph() {
        GraphTable graph = ProcessGraph.parserFile("USA-road-d.COL.co");
        //System.out.println("Adding weight to the graph...");
        graph = ProcessGraph.applyArcs(graph, 4L, "USA-road-d.COL.gr");
        //System.out.println("Adding weight to the graph...");
        graph = ProcessGraph.applyArcs(graph, 0L, "USA-road-t.COL.gr");
        graph = ProcessGraph.computeNewWeight(graph, 1L, 5L, 4L, 0L);

        //ProcessGraph.printSCCs(graph);

        System.out.println("Normalizing...");
        graph = ProcessGraph.normalizate(graph);
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

        graph = ProcessGraph.normalizate(graph);

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
        System.out.println("Loading graph...");
        GraphTable graph = ProcessGraph.parserFile(graphFilePath);
        System.out.println("Adding weight to the graph...");
        graph = ProcessGraph.applyWeights(graph, weightFilePath0);
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
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
}
