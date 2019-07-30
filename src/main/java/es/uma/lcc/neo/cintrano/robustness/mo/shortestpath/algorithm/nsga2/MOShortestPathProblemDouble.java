package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2;

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
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.solution.impl.DefaultDoubleSolution;
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

    private static final int MAX_SAMPLES = 30;

    private GraphTable graph;
    private Long start;
    private Long end;
    private int numVariables;
    private int numObjectives = 4;
    private final int numConstraints = 0;
    private Map<String, Integer> objectivesTags; // <tag, fitness array id>

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

    /*
      No traffic light version
     */
//    public void evaluate(DoubleSolution pathSolution) {
//        float[] fitness = new float[4];//pathSolution.getNumberOfObjectives()];
//        for (int i = 0; i < fitness.length; i++) {
//            fitness[i] = 0f;
//        }
//        boolean continuar = true;
//        for (int i = 0; i < pathSolution.getNumberOfVariables() - 1; i++) {
//            if ((pathSolution.getVariableValue(i)).longValue() == end) {
//                continuar = false;
//            }
//            if (continuar){
//                Long arco = graph.getAdjacencyMatrix().get((pathSolution.getVariableValue(i)).longValue(), (pathSolution.getVariableValue(i + 1)).longValue());
//                for (int j = 0; j < fitness.length; j++) {
//                    fitness[j] += graph.getWeightsMatrix().get(arco, (long) j);
//                }
//            }
//        }
//        for (int i = 0; i < fitness.length; i++) {//fitness.length; i++) {
//            pathSolution.setObjective(i, fitness[i]);
//        }
//    }

    /**
     * Traffic light version
     * @param pathSolution path to be evaluate
     */
    public void evaluate(DoubleSolution pathSolution) {
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
            if (continuar){
                n1 = pathSolution.getVariableValue(i).longValue();
                n2 = (pathSolution.getVariableValue(i + 1)).longValue();
                Long arc = graph.getAdjacencyMatrix().get(n1, n2);
                float mu = graph.getWeightsMatrix().get(arc, 0L); // TODO Change to the final weight
                float sigma = graph.getWeightsMatrix().get(arc, 1L); // TODO Change to the final weight
                float tentativeTime = (float) randNormal(mu, sigma);
                System.out.println("------- JMetal Gaussian Random: " + tentativeTime);
                // TrafficLight phase
                TlLogic tl = graph.getTlMatrix().get(n1, n2);
                cost += tentativeTime; // TODO Add drive profile
                if (tl != null) {
                    int time = tl.calculateTimeStop(Math.round(cost + 0.5d));
                    cost += time;
                }
            }
        }
            samples.add(cost);
            amount += cost;
        }
        float mean = amount / (float) MAX_SAMPLES;
        pathSolution.setObjective(0, mean);
        pathSolution.setObjective(1, variance(samples, mean));
    }

    private float variance(List<Double> samples, double mean) {
        float var = 0;
        for (Double sample : samples) {
            var += Math.pow(sample - mean, 2f);
        }
        var = var / (samples.size() -1);
        return var;
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
        return y1;
    }


    public DoubleSolution createSolution() {
        JMetalRandom randomGenerator = JMetalRandom.getInstance();

        float peso = 0.5f;
        peso = new Double(randomGenerator.nextDouble()).floatValue();
        System.out.println(peso);

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
            System.out.println(aux);
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
        System.out.println("=== NEW INDIVIDUAL ===");
        for (int i = 0; i < path.size(); i++) {
            pathArray[i] = path.get(i);
            System.out.println(pathArray[i]);
        }
        System.out.println("======");
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
            System.out.print(pathArray[i] + " ");
        //    System.out.println(pathArray[i]);
        }
        System.out.println();
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
        System.out.println("=== NEW INDIVIDUAL ===");
        for (int i = 0; i < path.size(); i++) {
            pathArray[i] = path.get(i);
            System.out.println(pathArray[i]);
        }
        System.out.println("======");
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
