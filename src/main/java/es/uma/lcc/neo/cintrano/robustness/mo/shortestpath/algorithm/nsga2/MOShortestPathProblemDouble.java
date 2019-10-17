package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.MyUtility;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.astar.Astar;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWeighted;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.TlLogic;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.utilities.ProcessGraph;
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Stack;

/**
 * Created by cintrano on 22/12/16.
 * Multi Objective Shortest Path Problem
 */
public class MOShortestPathProblemDouble implements Problem<DoubleSolution>, DoubleProblem {

    private static int MAX_SAMPLES = 30;

    private GraphTable graph;
    private Long start;
    private Long end;
    private int numVariables;
    private int numObjectives = 4;
    private final int numConstraints = 0;
    private Map<String, Integer> objectivesTags; // <tag, fitness array id>
    private boolean robustFlag = true;
    private boolean TLFlag = true;

    private Algorithm<DoubleSolution> algorithm;
    public void setAlgorithm(Algorithm<DoubleSolution> algorithm) {
        this.algorithm = algorithm;
    }

    /**
     * Constructor.
     * Creates a default instance of the Bi-Objective Shortest Path problem.
     */
    public MOShortestPathProblemDouble() {

    }

    private MOShortestPathProblemDouble(GraphTable graph, Integer numVariables) {
        this.graph = graph;
        setNumVariables(numVariables);
    }

    public MOShortestPathProblemDouble(GraphTable graph) {
        this(graph, 10);
    }

    public MOShortestPathProblemDouble(GraphTable graph, Long start, Long end, int numVariables, Map<String, Integer> objectivesTags) {
        this(graph, numVariables);
        this.start = start;
        this.end = end;
        this.objectivesTags = objectivesTags;
    }
    public MOShortestPathProblemDouble(GraphTable graph, Long start, Long end, int numVariables, boolean robustFlag, boolean TLFlag) {
        this(graph, start, end, numVariables);
        this.robustFlag = robustFlag;
        this.TLFlag = TLFlag;
        if (!robustFlag) MAX_SAMPLES = 1;
    }

    public MOShortestPathProblemDouble(GraphTable graph, Long start, Long end, int numVariables) {
        this(graph, numVariables);
        this.start = start;
        this.end = end;
    }

    private void setNumVariables(int numVariables) {
        this.numVariables = numVariables;
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
    public void evaluateAnterior(DoubleSolution pathSolution) {
        for (int objective = 0; objective < numObjectives; objective += 2) {
            int tlCount = 0;
            long n1, n2;
            float amount = 0;
            List<Double> samples = new ArrayList<>();
            for (int sample = 0; sample < MAX_SAMPLES; sample++) {
                double cost = 0;
                boolean continuar = true;
                for (int i = 0; i < pathSolution.getNumberOfVariables() - 1; i++) {
                    if ((pathSolution.getVariableValue(i)).longValue() == end) {
                        continuar = false;
                    }
                    if (continuar) {
                        n1 = pathSolution.getVariableValue(i).longValue();
                        n2 = (pathSolution.getVariableValue(i + 1)).longValue();
                        Long arc = graph.getAdjacencyMatrix().get(n1, n2);
                        float mu = graph.getWeightsMatrix().get(arc, (long) objective);
                        float sigma = graph.getWeightsMatrix().get(arc, (long) objective + 1L);
                        float tentativeTime = (float) randNormal(mu, sigma);
//                        System.out.println("------- JMetal Gaussian Random: " + tentativeTime);
                        cost += tentativeTime; // TODO Add drive profile
                        if (TLFlag) {
                            // TrafficLight phase
                            TlLogic tl = graph.getTlMatrix().get(n1, n2);
                            if (tl != null) {
                                tlCount++;
                                int time = tl.calculateTimeStop(Math.round(cost + 0.5d));
                                cost += time;
                            }
                        }
                    }
                }
                samples.add(cost);
                amount += cost;
            }
            float mean = amount / (float) MAX_SAMPLES;
            float variance = MyUtility.variance(samples, mean);
            System.out.println("F " + mean + " " + variance + " " + tlCount);
            pathSolution.setObjective(objective, mean);
            pathSolution.setObjective(objective + 1, variance);
            pathSolution.setAttribute("tl", tlCount);
        }
    }

    public void evaluate(DoubleSolution pathSolution) {
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

                    for (int objective = 0; objective < numObjectives; objective += 2) {
                        mu[objective] = graph.getWeightsMatrix().get(arc, (long) objective);
                        sigma[objective] = graph.getWeightsMatrix().get(arc, (long) objective + 1);
                    }
                    float tentativeTime = (float) randNormal(mu[0], sigma[0]);
                    float tentativePollution = (float) randNormal(mu[1], sigma[1]);//mu[0] * (tentativeTime/mu[0]); // CO2
                    cost[0] = tentativeTime; // TODO Add drive profile
                    cost[1] = tentativePollution;
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
            fit[0] = (fit[0] * sample + cost[0]) / (sample + 1);
            fit[2] = (fit[2] * sample + cost[1]) / (sample + 1);
            // Variances
            sum[0] += cost[0];
            sumSq[0] += cost[0] * cost[0];
            sum[1] += cost[1];
            sumSq[1] += cost[1] * cost[1];
        }
        fit[1] = (sumSq[0] - (sum[0] * sum[0]) / (float) MAX_SAMPLES) / (float) MAX_SAMPLES;
        fit[3] = (sumSq[1] - (sum[1] * sum[1]) / (float) MAX_SAMPLES) / (float) MAX_SAMPLES;

        fit[fit.length - 1] = tlCount;
        pathSolution.setObjective(0, fit[0]);
        pathSolution.setObjective(1, fit[1]);
        pathSolution.setObjective(2, fit[2]);
        pathSolution.setObjective(3, fit[3]);
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


    public DoubleSolution createSolution() {
        JMetalRandom randomGenerator = JMetalRandom.getInstance();

        float peso = 0.5f;
        peso = new Double(randomGenerator.nextDouble()).floatValue();
//        System.out.println(peso);

        Long[] aux = getRandomAstarPath(this.graph, this.start, this.end);
        Double[] auxD = new Double[aux.length];

        for (int i = 0; i < aux.length; i++) {
            auxD[i] = aux[i].doubleValue();
        }
/*
        for (int i = 0; i < aux.length; i++) {
            d.setVariableValue(i, auxD[i]);
        }
        */
        MyDoubleSolution d = new MyDoubleSolution(new double[getNumberOfObjectives()], auxD)  ;
        /*
        DoubleSolution d = new DoubleSolution(new double[getNumberOfObjectives()], auxD) {
            double[] objectives;
            Double[] variables;

            public DoubleSolution(double[] objectives, Double[] variables) {
                this.objectives = objectives;
                this.variables = variables;
            }

            @Override
            public Double getLowerBound(int i) {
                double d = variables[0];
                for (int j = 0; j < variables.length; j++) {
                    if (d > variables[i]) d = variables[i];
                }
                return d;
            }

            @Override
            public Double getUpperBound(int i) {
                double d = variables[0];
                for (int j = 0; j < variables.length; j++) {
                    if (d < variables[i]) d = variables[i];
                }
                return d;
            }

            @Override
            public void setObjective(int i, double v) {
                this.objectives[i] = v;
            }

            @Override
            public double getObjective(int i) {
                return this.objectives[i];
            }

            @Override
            public double[] getObjectives() {
                return this.objectives;
            }

            @Override
            public Double getVariableValue(int i) {
                return this.variables[i];
            }

            @Override
            public void setVariableValue(int i, Double aDouble) {
                this.variables[i] = aDouble;
            }

            @Override
            public String getVariableValueString(int i) {
                return this.variables[i].toString();
            }

            @Override
            public int getNumberOfVariables() {
                return variables.length;
            }

            @Override
            public int getNumberOfObjectives() {
                return objectives.length;
            }

            @Override
            public Solution<Double> copy() {
                return null;
            }

            @Override
            public void setAttribute(Object o, Object o1) {

            }

            @Override
            public Object getAttribute(Object o) {
                return null;
            }
        };
        */

        return d;

    }

    private Long[] getRandomAstarPath(GraphTable graph, Long start, Long end) {
        JMetalRandom randomGenerator = JMetalRandom.getInstance();
        int index = randomGenerator.nextInt(0, graph.getAdjacencyMatrix().rowKeySet().size()-1);
        Long nextNode = new ArrayList<Long>(graph.getAdjacencyMatrix().rowKeySet()).get(index);

        Long[] head = getRandomAstarSearchPath(graph, start, nextNode);//getShortestPathBetween(start, nextNode);
        Long[] tail = getRandomAstarSearchPath(graph, nextNode, end);//getShortestPathBetween(nextNode, end);

        Long[] out = new Long[head.length + tail.length - 1];
            for (int i = 0; i < head.length; i++) {
                out[i] = head[i];
            }
            for (int i = 1; i < tail.length; i++) {
                out[i + head.length - 1] = tail[i];
            }
        return out;
    }

    private static Long[] getRandomPath(GraphTable graph, Long start, Long end) {
        List<Long> aux;
        Long current = start;
        List<Long> path = new ArrayList<Long>();
        path.add(start);
        JMetalRandom randomGenerator = JMetalRandom.getInstance();
        boolean valid;
        do {
            aux = new ArrayList<Long>(graph.getAdjacencyMatrix().row(current).keySet());
//            System.out.println(aux);
            valid = false;
            while (aux.size() > 0 && !valid) {
                int index = randomGenerator.nextInt(0, aux.size()-1);
                //System.out.println(aux.size() + " " + index);
                if (path.contains(aux.get(index))) {
                    // remove next node
                    aux.remove(index);
                } else {
                    path.add(aux.get(index));
                    current = aux.get(index);
                    valid = true;
                }
            }
        } while (!current.equals(end) && aux.size() > 0);
        Long[] pathArray = new Long[path.size()];
//        System.out.println("=== NEW INDIVIDUAL ===");
        for (int i = 0; i < path.size(); i++) {
            pathArray[i] = path.get(i);
//            System.out.println(pathArray[i]);
        }
//        System.out.println("======");
        return pathArray;
    }

    private static Long[] getRandomDijkstraSearchPath(GraphTable graph, Long start, Long end, float peso) {
        List<Long> path;
        DijkstraWeighted algorithm = new DijkstraWeighted(0L, 1L, peso);
        algorithm.setGraph(graph);

        /*
        Long middleNode = new ArrayList<Long>(graph.getIntersections().keySet()).get(randomGenerator.nextInt(0, graph.getIntersections().keySet().size() -1));
        path = algorithm.getPath(graph.getIntersections().get(start), graph.getIntersections().get(middleNode));
        path.remove(path.size()-1);
        path.addAll(algorithm.getPath(graph.getIntersections().get(middleNode), graph.getIntersections().get(end)));
        */
        List<Node> pathN = algorithm.getPath(start, end);
        path = ProcessGraph.nodeToLong(graph, pathN);

        Long[] pathArray = new Long[path.size()];
        //System.out.println("=== NEW INDIVIDUAL ===");
        for (int i = 0; i < path.size(); i++) {
            pathArray[i] = path.get(i);
//            System.out.print(pathArray[i] + " ");
        //    System.out.println(pathArray[i]);
        }
//        System.out.println();
        //System.out.println("======");
/*
        System.out.print("N " + pathArray.length + " --> ");
        for (Long l : pathArray) {
            System.out.print(" " + l);
        }
        System.out.println();
*/
        return pathArray;
    }
    private static Long[] getRandomAstarSearchPath(GraphTable graph, Long start, Long end) {
        List<Node> path;
        Astar algorithm = new Astar();
        algorithm.setGraph(graph);
        algorithm.setTarget(0.5f);

        /*
        Long middleNode = new ArrayList<Long>(graph.getIntersections().keySet()).get(randomGenerator.nextInt(0, graph.getIntersections().keySet().size() -1));
        path = algorithm.getPath(graph.getIntersections().get(start), graph.getIntersections().get(middleNode));
        path.remove(path.size()-1);
        path.addAll(algorithm.getPath(graph.getIntersections().get(middleNode), graph.getIntersections().get(end)));
        */
        path = algorithm.getPath(start, end);

        //System.out.print("PATH " + path.size() + " " + start + " " + end+ " -> ");
        //for (IntersectionInterface l : path) {
        //    System.out.print(l + " ");
        //}
        //System.out.println();


        Long[] pathArray = new Long[path.size()];
        //System.out.println("=== NEW INDIVIDUAL ===");
        for (int i = 0; i < path.size(); i++) {
            pathArray[i] = graph.getMapping().get(path.get(i).getId());
            //    System.out.println(pathArray[i]);
        }
        //System.out.println("======");
/*
        System.out.print("N " + pathArray.length + " --> ");
        for (Long l : pathArray) {
            System.out.print(" " + l);
        }
        System.out.println();
*/
        return pathArray;
    }


    private static Long[] getDepthFirstSearchPath(GraphTable graph, Long start, Long end) {
        List<Long> path = new ArrayList<Long>();
        List<Long> aux;

        Stack<Long> stack = new Stack<Long>();
        stack.add(start);
        boolean foundGoal = false;

        while (!stack.isEmpty() && !foundGoal) {
            Long node = stack.pop();

            if (path.size() > 1 && isCorrelative(graph, path.get(path.size()-1), node)) {
                if (node.equals(end)) {
                    path.add(node);
                    foundGoal = true;
                } else {
                    path.add(node);
                    aux = new ArrayList<Long>(graph.getAdjacencyMatrix().row(node).keySet());
                    for (Long n : aux) {
                        stack.push(n);
                    }
                }
            } else if (path.isEmpty()) {
                path.add(node);
                aux = new ArrayList<Long>(graph.getAdjacencyMatrix().row(node).keySet());
                for (Long n : aux) {
                    stack.push(n);
                }
            } else {
                path.remove(path.size()-1);
                stack.push(node);
            }

        }

        Long[] pathArray = new Long[path.size()];
//        System.out.println("=== NEW INDIVIDUAL ===");
        for (int i = 0; i < path.size(); i++) {
            pathArray[i] = path.get(i);
//            System.out.println(pathArray[i]);
        }
//        System.out.println("======");
        return pathArray;
    }

    private static boolean isCorrelative(GraphTable graph, Long start, Long end) {
        return graph.getAdjacencyMatrix().row(start).keySet().contains(end);
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

    /**
     * Dijkstra method to compute a path between the node and the end-node
     * @param node start node
     * @return path
     */
    public Long[] getShortestPartialPath(Long node) {
        // TODO: change to dijkstra algorithm
        return getRandomPath(graph, node, end);
    }

    public Long[] getShortestPathWithMiddlePoint() {
        JMetalRandom randomGenerator = JMetalRandom.getInstance();
        Long middle = new ArrayList<Long>(this.graph.getIntersections().keySet()).get(randomGenerator.nextInt(0, this.graph.getIntersections().keySet().size() - 1));
        Long[] head = getRandomDijkstraSearchPath(this.graph, this.start, middle, 0.5f);
        Long[] tail = getRandomDijkstraSearchPath(this.graph, middle, this.end, 0.5f);
        Long[] out = new Long[head.length + tail.length -1];
        for (int i = 0; i < head.length -1; i++) {
            out[i] = head[i];
        }
        for (int i = 0; i < tail.length; i++) {
            out[i + head.length - 1] = tail[i];
        }
        return out;
    }

    public Long[] getShortestPathWithMiddlePoint(Long variableValue) {
        //return getRandomDijkstraSearchPath(this.graph, variableValue, this.end, 0.5f);
        return getRandomAstarSearchPath(this.graph, variableValue, this.end);
    }

    public Long findRandomNextNode(Long variableValue) {
        JMetalRandom randomGenerator = JMetalRandom.getInstance();
        return new ArrayList<Long>(graph.getAdjacencyMatrix().row(variableValue).keySet())
                .get(randomGenerator.nextInt(0,graph.getAdjacencyMatrix().row(variableValue).keySet().size()-1));
    }

    public GraphTable getGraph() {
        return graph;
    }

    public Long[] getShortestPathBetween(Long from, Long to) {
        return getRandomAstarSearchPath(this.graph, from, to);
    }

    private List<Double> lowerLimit ;
    private List<Double> upperLimit ;

    /* Getters */
    @Override
    public Double getUpperBound(int index) {
        return upperLimit.get(index);
    }

    @Override
    public Double getLowerBound(int index) {
        return lowerLimit.get(index);
    }

    /* Setters */
    protected void setLowerLimit(List<Double> lowerLimit) {
        this.lowerLimit = lowerLimit;
    }

    protected void setUpperLimit(List<Double> upperLimit) {
        this.upperLimit = upperLimit;
    }
}
