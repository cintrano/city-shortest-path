package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.astar.Astar;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWeighted;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.utilities.ProcessGraph;
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.util.*;

/**
 * Created by cintrano on 22/12/16.
 * Multi Objective Shortest Path Problem
 */
public class MOShortestPathProblem implements Problem<NodePathSolution> {
    private GraphTable graph;
    private Long start;
    private Long end;
    private int numVariables;
    private int numObjectives = 4;
    private final int numConstraints = 0;
    private Map<String, Integer> objectivesTags; // <tag, fitness array id>

    private Algorithm<NodePathSolution> algorithm;
    public void setAlgorithm(Algorithm<NodePathSolution> algorithm) {
        this.algorithm = algorithm;
    }

    /**
     * Constructor.
     * Creates a default instance of the Bi-Objective Shortest Path problem.
     */
    public MOShortestPathProblem() {

    }

    private MOShortestPathProblem(GraphTable graph, Integer numVariables) {
        this.graph = graph;
        setNumVariables(numVariables);
    }

    public MOShortestPathProblem(GraphTable graph) {
        this(graph, 10);
    }

    public MOShortestPathProblem(GraphTable graph, Long start, Long end, int numVariables, Map<String, Integer> objectivesTags) {
        this(graph, numVariables);
        this.start = start;
        this.end = end;
        this.objectivesTags = objectivesTags;
    }
    public MOShortestPathProblem(GraphTable graph, Long start, Long end, int numVariables) {
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

    public void evaluate(NodePathSolution pathSolution) {
        float[] fitness = new float[4];//pathSolution.getNumberOfObjectives()];
        for (int i = 0; i < fitness.length; i++) {
            fitness[i] = 0f;
        }
        //System.out.println("evaluating " + pathSolution);
        for (int i = 0; i < pathSolution.getNumberOfVariables() - 1; i++) {
            Long arco = graph.getAdjacencyMatrix().get(pathSolution.getVariableValue(i), pathSolution.getVariableValue(i+1));
            //System.out.println(arco);
            for (int j = 0; j < fitness.length; j++) {
                fitness[j] += graph.getWeightsMatrix().get(arco, (long) j);
            }
            /*
            fitness[1] = 0f;
            fitness[2] = 0f;
            fitness[3] = 0f;
            /*
            fitness[0] += graph.getWeightsMatrix().row(arco).get(0L);
            fitness[1] += graph.getWeightsMatrix().row(arco).get(1L);
            fitness[2] += graph.getWeightsMatrix().row(arco).get(2L);
            fitness[3] += graph.getWeightsMatrix().row(arco).get(3L);
            */
            /*
            weights = graph.getWeights(pathSolution.getVariableValue(i), pathSolution.getVariableValue(i+1));
            for (Weight weight: weights) {
                fitness[objectivesTags.get(weight.getType())] += weight.getValue();
            }
            */
        }
        /*
        fitness[0] = (fitness[0] * 0.5f) + (fitness[1] * 0.5f);
        List<Double> lista = new ArrayList<Double>();
        for (Object n : ((SteadyStateGeneticAlgorithm) algorithm).getPopulation()) {
            lista.add(((NodePathSolution) n).getObjective(0));
        }
        Collections.sort(lista);
        System.out.println(lista.get(4));
        */
        // Penalty
        // TODO: Dijkstra penalty


        // Set fitness in the individual
        //System.out.print("i--> ");
        for (int i = 0; i < fitness.length; i++) {//fitness.length; i++) {
            pathSolution.setObjective(i, fitness[i]);
        //    System.out.print(fitness[i] + " ");
        }

        //((SteadyStateGeneticAlgorithm) algorithm).getPopulation();
        //System.out.println();
    }

    public NodePathSolution createSolution() {
        JMetalRandom randomGenerator = JMetalRandom.getInstance();

        float peso = 0.5f;
        peso = new Double(randomGenerator.nextDouble()).floatValue();
        System.out.println(peso);
        return new NodePathSolution(new double[getNumberOfObjectives()],
                //getRandomPath(this.graph, this.start, this.end));
                //getDepthFirstSearchPath(this.graph, this.start, this.end));
                //getRandomDijkstraSearchPath(this.graph, this.start, this.end, peso));
                getRandomAstarPath(this.graph, this.start, this.end));
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
}
