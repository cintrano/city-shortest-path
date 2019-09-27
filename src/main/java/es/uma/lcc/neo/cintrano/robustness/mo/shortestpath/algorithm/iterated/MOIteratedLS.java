package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.iterated;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.MyUtility;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.RunTLMain;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDim;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDimForbidden;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.TlLogic;

import java.util.*;

public class MOIteratedLS {

    private static int MAX_SAMPLES = 30;
    private GraphTable graph;
    private Random rand;
    private int maxIterations = 100;
    private boolean robustFlag = true;
    private boolean TLFlag = true;
    private int numObjectives = 4;

    private MOIteratedLS() {
        rand = new Random();
    }

    public MOIteratedLS(int seed) {
        this();
        rand.setSeed(seed);
    }

    public void setMaxIterations(int maxIterations) {
        this.maxIterations = maxIterations;
    }

    public void setGraph(GraphTable graph) {
        this.graph = graph;
    }

    public Set<NodePathSolution> getPath(Long start, Long end, float[] weights) {
        // Dijkstra algorithm helper
        DijkstraWNDim dj = new DijkstraWNDim(RunTLMain.objectives);
        dj.setWeights(weights);
        dj.setGraph(graph);

        List<Node> s0 = dj.getPath(start, end); // initial solution
        List<Node> current, best = new ArrayList<>();
        float[] fCurrent;
        Set<NodePathSolution> nonDominated = new HashSet<>();
        current = s0;
        best = s0;
        fCurrent = fitness(current);
        nonDominated.add(new NodePathSolution(fCurrent, node2arrayLong(current)));
        int iteration = 0;
        while (iteration < maxIterations) {
            current = mutate(best, weights); // weights is a Dijkstra param
            fCurrent = fitness(current);
            if (updateNonDominated(nonDominated, fCurrent, current)) {
                best = current;
            }
            iteration++;
        }
        return nonDominated;
    }

    public void setRobustFlag(boolean robustFlag) {
        this.robustFlag = robustFlag;
        if (!robustFlag) MAX_SAMPLES = 1;
    }

    public void setTLFlag(boolean TLFlag) {
        this.TLFlag = TLFlag;
    }

    private boolean updateNonDominated(Set<NodePathSolution> nonDominated, float[] fitness, List<Node> values) {
        for (NodePathSolution n : nonDominated) {
            boolean dom = true;
            for (int i = 0; i < fitness.length; i++) {
                dom = dom && (fitness[i] <= n.getObjectives()[i]);
            }
            if (dom) {
                nonDominated.add(new NodePathSolution(fitness, node2arrayLong(values)));
                return true;
            }
        }
        return false;
    }

    private Long[] node2arrayLong(List<Node> nodes) {
        Long[] indexes = new Long[nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            indexes[i] = nodes.get(i).getId();
        }
        return indexes;
    }

    private List<Node> mutate(List<Node> s, float[] weights) {
        Set<Long> neighbours;
        int point;
        do {
            point = rand.nextInt(s.size() - 1);
            neighbours = graph.getAdjacencyMatrix().row(graph.getMapping().get(s.get(point).getId())).keySet();
        } while (neighbours.size() < 2);
        List<Node> newSolution = new ArrayList<>(s.subList(0, point+1));

        // Dijkstra algorithm helper
        // TODO: Change the implementation of Dijkstra to take into account a set of forbidden nodes
        DijkstraWNDimForbidden dj = new DijkstraWNDimForbidden(RunTLMain.objectives);
        dj.setWeights(weights);
        dj.setGraph(graph);
        dj.setForbidden(newSolution);

        // Select the new node
        List<Long> neighboursList = new ArrayList<>(neighbours);
        Collections.shuffle(neighboursList);
        int j = 0;
        while (j < neighboursList.size() && neighboursList.get(j).equals(graph.getMapping().get(s.get(point + 1).getId()))) {
            j++;
        }
        List<Node> subPath = dj.getPath(neighboursList.get(j), graph.getMapping().get(s.get(s.size()-1).getId()));
        newSolution.addAll(subPath);
        return newSolution;
    }

    // Fitness robust mono-objective
    private float[] fitness(List<Node> s) {
        float[] fit = new float[numObjectives];
        for (int objective = 0; objective < numObjectives; objective += 2) {
            int tlCount = 0;
            long n1, n2;
            float amount = 0;
            List<Double> samples = new ArrayList<>();
            for (int sample = 0; sample < MAX_SAMPLES; sample++) {
                double cost = 0;
                for (int i = 0; i < s.size() - 1; i++) {
                    n1 = graph.getMapping().get(s.get(i).getId());
                    n2 = graph.getMapping().get(s.get(i + 1).getId());
                    long arc = graph.getAdjacencyMatrix().get(n1, n2);
                    float mu = graph.getWeightsMatrix().get(arc, (long) objective);
                    float sigma = graph.getWeightsMatrix().get(arc, (long) objective + 1);
                    float tentativeTime = (float) rand.nextGaussian() * sigma + mu;
                    cost += tentativeTime; // TODO Add drive profile
                    // TrafficLight phase
                    if (TLFlag) {
                        TlLogic tl = graph.getTlMatrix().get(n1, n2);
                        if (tl != null) {
                            tlCount++;
                            int time = tl.calculateTimeStop(Math.round(cost + 0.5d));
                            cost += time;
                        }
                    }
                }
                samples.add(cost);
                amount += cost;
            }
            float mean = amount / (float) MAX_SAMPLES;
            float variance = MyUtility.variance(samples, mean);
            System.out.println("F " + mean + " " + variance + " " + tlCount);
            fit[objective] = mean;
            fit[objective+1] = variance;
        }
        return fit;
    }
}
