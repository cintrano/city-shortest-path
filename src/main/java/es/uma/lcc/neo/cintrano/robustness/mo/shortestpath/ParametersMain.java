package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDim;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.utilities.ProcessGraph;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Created by Christian Cintrano on 23/05/17.
 * Parameters Main
 */
public class ParametersMain {

    private static final int[] SEED = new int[]{0, 3, 23, 29, 67, 79, 97, 101, 139, 199, 211, 269, 311, 347, 389,
            431, 461, 503, 569, 601, 607, 653, 691, 701, 733, 739, 761, 809, 811, 997};

    private static final int WEIGHTS_SIZE = 5;
    private static final int OBJECTIVE_SIZE = 4;

    private static int weightSeed, percentString, ordinationSeed;

    public static void main (String[] args) {
        System.out.println("=== START EXPERIMENTS ===");
        for (String s : args) {
            System.out.print(s + " ");
        }
        System.out.println();

        String map;
        int pointId;

        if (args.length == 4) {
            map = args[0];
            GraphTable graph = getGraphMap(map);
            pointId = Integer.parseInt(args[1]);
            weightSeed = Integer.parseInt(args[2]);
            percentString = Integer.parseInt(args[3]);
            ordinationSeed = Integer.parseInt(args[4]);
            executeProcess(graph, map, pointId, weightSeed, percentString, ordinationSeed);
        }
        if (args.length == 2) {// Only Random options
            map = args[0];
            GraphTable graph = getGraphMap(map);

            int number = Integer.parseInt(args[1]);
            int tens = number / 10;
            int units = number % 10; // point
            int weightLB = (tens / 3) * 10;
            int orderLB = (tens % 3) * 10;
            for (int i = weightLB; i < (weightLB + 10); i++) {
                for (int j = orderLB; j < (orderLB + 10); j++) {
                    executeProcess(graph, map, units, i, 100, j);
                }
            }
        }

        System.out.println("=== END EXPERIMENTS ===");
    }

    private static void executeProcess(GraphTable graph, String map, int pointId, int weightSeed, int percentage, int ordinationSeed) {

        // Points
        System.out.print("Reading points...");
        Long[] points = RunMain.readPoints(pointId, map);
        System.out.println("end");

        // Algorithm
        if (points != null) {
            System.out.println(points[0] + " " + points[1]);

            System.out.print("Calculate the " + WEIGHTS_SIZE + " weights...");
            float[] weights = getWeights(weightSeed);
            System.out.println("end");

            System.out.print("Get ordination...");
            List<float[]> order = getOrdination(ordinationSeed, weights, percentage);
            System.out.println("...end");

            System.out.println("Run Dijkstra...");
            rDijkstraExperiments(graph, points, order);
        }
    }

    private static GraphTable getGraphMap(String map) {
        GraphTable graph = null;
        // Map
        switch (map) {
            case "Malaga":
                System.out.println("Loading graph...");
                graph = prepareGraph("hbefa-malaga-graph.xml", "weights_time-hbefa.xml", "mapping-malaga.txt", "MAL");
                graph = ProcessGraph.fixVertexIndex(graph);
                break;
            case "Colorado":
                graph = prepareGraph("hbefa-col-graph.xml", "col_weights_time-hbefa-MOD.xml", "mapping-col.txt", "COL");
                //graph = prepareSimpleColoradoGraph();
                break;
            case "NY":
                graph = prepareGraph("hbefa-ny-graph.xml", "ny_weights_time-hbefa.xml", "mapping-ny.txt", "NY");
                //graph = prepareSimpleColoradoGraph();
                break;
        }
        return graph;
    }

    private static List<float[]> getOrdination(int arg, float[] weights, int percentage) {
        int sel = 1;
        for (int i = 0; i < OBJECTIVE_SIZE; i++) {
            sel = sel * WEIGHTS_SIZE;
        }
        sel = (sel * percentage) / 100;

        List<float[]> order = new ArrayList<>();
        if (arg == -1) {
            int iter = 0;
            for (float w0 : weights) {
                for (float w1 : weights) {
                    for (float w2 : weights) {
                        for (float w3 : weights) {
                            if (w0 + w1 + w2 + w3 > 0) {
                                order.add(new float[]{w0, w1, w2, w3});
                                iter++;
                                if (iter == sel) {
                                    return order;
                                }
                            }
                        }
                    }
                }
            }
            return order;
        } else {
            List<float[]> combination = new ArrayList<>();
            for (float w0 : weights) {
                for (float w1 : weights) {
                    for (float w2 : weights) {
                        for (float w3 : weights) {
                            if (w0 + w1 + w2 + w3 > 0) {
                                combination.add(new float[]{w0, w1, w2, w3});
                            }
                        }
                    }
                }
            }

            Random random = new Random(SEED[arg]);
            int index;
            for (int i = 0; i < sel; i++) {
                index = combination.size() == 0 ? 0 : random.nextInt(combination.size());
                order.add(combination.get(index));
                combination.remove(index);
            }
            return order;
        }
    }

    private static float[] getWeights(int arg) {
        if (arg == -1) {
            return new float[]{0.0001f, 0.25f, 0.5f, 0.75f, 1.0f};
        } else {
            Random random = new Random(SEED[arg]);
            float[] weights = new float[WEIGHTS_SIZE];
            for (int i = 0; i < WEIGHTS_SIZE; i++) {
                weights[i] = random.nextFloat();
            }
            return weights;
        }
    }

    private static void rDijkstraExperiments(GraphTable graph, Long[] randomPoints, List<float[]> order) {
        Long[] objectives = new Long[]{0L, 1L, 2L, 3L};
        System.out.println("Number of combinations " + order.size());

        for (float[] weights : order) {
            System.out.println("weights: " + weights[0] + " " + weights[1] + " " + weights[2] + " " + weights[3]);
            rDijkstraWD(graph, randomPoints, objectives, weights);
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

    private static void printSolutions(GraphTable graph, NodePathSolution s, long t1, long t2, String algorithm, String i) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("resultsFUN.ssv", true));
            String line = graph.getMapping().get(s.getVariables()[0]) + " ";
            line += graph.getMapping().get(s.getVariables()[s.getVariables().length-1]) + " ";
            line += t1 + " ";
            line += t2 + " ";
            line += algorithm + " ";
            line += weightSeed + " ";
            line += percentString + " ";
            line += ordinationSeed + " ";
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

    private static GraphTable prepareGraph(String graphFilePath, String weightFilePath0, String mapping, String tag) {
        // Graph
        //String graphFilePath = "hbefa-malaga-graph.xml";//"new-malaga-graph.xml";//"graph_connected.xml";
        //String weightFilePath0 = "weights_time-hbefa.xml";//"weights_time-noise.xml";//"wNew.xml";
        GraphTable graph = ProcessGraph.parserFile(graphFilePath);
        if (tag.equals("COL")) {
            graph = ProcessGraph.readWeights(graph, weightFilePath0);
        } else {
            graph = ProcessGraph.applyWeights(graph, weightFilePath0);
        }
        //graph = ProcessGraph.applyWeights(graph, "wVar0.xml");
        //graph = ProcessGraph.applyWeights(graph, "wVar1.xml");
        graph.getWeightsMatrix().column(10L).clear();
        graph = ProcessGraph.applyMapping(graph, mapping);
        return graph;
    }

}
