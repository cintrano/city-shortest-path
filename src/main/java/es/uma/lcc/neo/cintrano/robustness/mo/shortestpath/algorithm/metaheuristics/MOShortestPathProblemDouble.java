package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.astar.Astar;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.TlLogic;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by cintrano on 22/12/16.
 * Multi Objective Shortest Path Problem
 */
public class MOShortestPathProblemDouble implements Problem<MyDoubleSolution> {

    private static int MAX_SAMPLES = 30;

    private GraphTable graph;
    private Long start;
    private Long end;
    private int numVariables;
    private int numObjectives = 4;
    private final int numConstraints = 0;
    private boolean robustFlag = true;
    private boolean TLFlag = true;

    /**
     * Constructor.
     * Creates a default instance of the Bi-Objective Shortest Path problem.
     */
    private MOShortestPathProblemDouble(GraphTable graph, Integer numVariables) {
        this.graph = graph;
        this.numVariables = numVariables;
    }

    public MOShortestPathProblemDouble(GraphTable graph, Long start, Long end, int numVariables) {
        this(graph, numVariables);
        this.start = start;
        this.end = end;
    }

    public MOShortestPathProblemDouble(GraphTable graph, Long start, Long end, int numVariables, boolean robustFlag, boolean TLFlag) {
        this(graph, start, end, numVariables);
        this.robustFlag = robustFlag;
        this.TLFlag = TLFlag;
        if (!robustFlag) MAX_SAMPLES = 1;
    }

    public int getNumberOfVariables() {
        return numVariables;
    }

    public int getNumberOfObjectives() {
        return numObjectives;
    }

    public int getNumberOfConstraints() {
        return numConstraints;
    }

    public String getName() {
        return "MOShortestPathProblem";
    }

    /**
     * Traffic light version
     * @param pathSolution path to be evaluate
     */
    public void evaluate(MyDoubleSolution pathSolution) {
        // System.out.println("EVALUATION:: " + pathSolution);
        float[] fit = new float[numObjectives + 1]; // To return the num of TL as the last one
        int tlCount = 0;
        long n1, n2;
        float[] mu, sigma;
        float[] sum = new float[2];
        float[] sumSq = new float[2];

        for (int sample = 0; sample < MAX_SAMPLES; sample++) {
            float[] cost = new float[numObjectives/2];
            boolean continuar = true;
            for (int i = 0; i < pathSolution.getNumberOfVariables() - 1; i++) {
                if ((pathSolution.getVariableValue(i)).longValue() == end) {
                    continuar = false;
                }
                if (continuar) {
                    n1 = pathSolution.getVariableValue(i).longValue();
                    n2 = (pathSolution.getVariableValue(i + 1)).longValue();
                    long arc = graph.getAdjacencyMatrix().get(n1, n2);
                    mu = new float[numObjectives / 2];
                    sigma = new float[numObjectives / 2];

                    for (int objective = 0; objective < numObjectives/2; objective ++) {
                        mu[objective] = graph.getWeightsMatrix().get(arc, (long) objective);
                        sigma[objective] = graph.getWeightsMatrix().get(arc, (long) objective + 1);
                    }
                    float tentativeTime;
                    float tentativePollution;
                    if (robustFlag) {
                        tentativeTime = (float) randNormal(mu[0], sigma[0]);
                        tentativePollution = (float) randNormal(mu[1], sigma[1]);//mu[0] * (tentativeTime/mu[0]); // CO2
                    } else {
                        tentativeTime = mu[0];
                        tentativePollution = mu[1];
                    }
                    cost[0] += tentativeTime; // TODO Add drive profile
                    cost[1] += tentativePollution;
                    // TrafficLight phase
                    if (TLFlag) {
                        TlLogic tl = graph.getTlMatrix().get(n1, n2);
                        if (tl != null) {
                            tlCount++;
                            int time = tl.calculateTimeStop(Math.round(cost[0] + 0.5d));
                            cost[0] += time;
                            float pollution = 2166.094f * (float) time; // c0 * time
                            cost[1] += pollution;
                        }
                    }
                }
            }
            // MEANS
            if (robustFlag) {
                fit[0] = (fit[0] * sample + cost[0]) / (sample + 1);
                fit[2] = (fit[2] * sample + cost[1]) / (sample + 1);
                // Variances
                sum[0] += cost[0];
                sumSq[0] += cost[0] * cost[0];
                sum[1] += cost[1];
                sumSq[1] += cost[1] * cost[1];
            } else {
                fit[0] = cost[0];
                fit[2] = cost[1];
            }
        }
        if (robustFlag) {
            fit[1] = (sumSq[0] - (sum[0] * sum[0]) / (float) MAX_SAMPLES) / (float) MAX_SAMPLES;
            fit[3] = (sumSq[1] - (sum[1] * sum[1]) / (float) MAX_SAMPLES) / (float) MAX_SAMPLES;
        } else {
            fit[1] = 0;
            fit[3] = 0;
        }

        fit[fit.length - 1] = tlCount;
        for (int i = 0; i < pathSolution.getNumberOfObjectives(); i++) {
            pathSolution.setObjective(i, fit[i]);
        }
        pathSolution.setAttribute("tl", tlCount);
    }

    /**
     *  Use the polar form of the Box-Muller transformation to obtain
     * a pseudo random number from a Gaussian distribution
     * Code taken from Maurice Clerc's implementation
     *
     * @param mean Mean
     * @param standardDeviation Sigma
     * @return A pseudo random number
     */
    private double randNormal(double mean, double standardDeviation) {
        double x1, x2, w, y1;
        do {
            x1 = 2.0 * JMetalRandom.getInstance().nextDouble() - 1.0;
            x2 = 2.0 * JMetalRandom.getInstance().nextDouble() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);
        w = Math.sqrt((-2.0 * Math.log(w)) / w);
        y1 = x1 * w * (standardDeviation + mean);
        return Math.abs(y1);
    }


    public MyDoubleSolution createSolution() {
        Long[] aux = getRandomAstarPath(this.graph, this.start, this.end);
        Double[] auxD = new Double[aux.length];
        for (int i = 0; i < aux.length; i++) {
            auxD[i] = aux[i].doubleValue();
        }
        MyDoubleSolution d = new MyDoubleSolution(new double[getNumberOfObjectives()], auxD);
        evaluate(d);
        return d;
    }

    private Long[] getRandomAstarPath(GraphTable graph, Long start, Long end) {
        JMetalRandom randomGenerator = JMetalRandom.getInstance();
        int index = randomGenerator.nextInt(0, graph.getAdjacencyMatrix().rowKeySet().size()-1);
        Long nextNode = new ArrayList<>(graph.getAdjacencyMatrix().rowKeySet()).get(index);

        Long[] head = getRandomAstarSearchPath(graph, start, nextNode); //getShortestPathBetween(start, nextNode);
        Long[] tail = getRandomAstarSearchPath(graph, nextNode, end); //getShortestPathBetween(nextNode, end);

        Long[] out = new Long[head.length + tail.length - 1];
        System.arraycopy(head, 0, out, 0, head.length);
        if (tail.length - 1 >= 0) System.arraycopy(tail, 1, out, 1 + head.length - 1, tail.length - 1);
        return out;
    }

    private static Long[] getRandomAstarSearchPath(GraphTable graph, Long start, Long end) {
        List<Node> path;
        Astar algorithm = new Astar();
        algorithm.setGraph(graph);
        algorithm.setTarget(0.5f);

        path = algorithm.getPath(start, end);

        Long[] pathArray = new Long[path.size()];
        for (int i = 0; i < path.size(); i++) {
            pathArray[i] = graph.getMapping().get(path.get(i).getId());
        }
        return pathArray;
    }


    @Override
    public String toString() {
        return "BiObjectiveShortestPathProblem{" +
                "graph=" + graph +
                ", start=" + start +
                ", end=" + end +
                ", numVariables=" + numVariables +
                ", numObjectives=" + numObjectives +
                ", numConstraints=" + numConstraints +
                '}';
    }

    public GraphTable getGraph() {
        return graph;
    }

    public Long[] getShortestPathBetween(Long from, Long to) {
        return getRandomAstarSearchPath(this.graph, from, to);
    }
}
