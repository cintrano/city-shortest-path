package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.iterated;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.RunTLMain;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDim;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWNDimForbidden;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.TlLogic;

import java.util.*;

public class IteratedLS {

    private static final int MAX_SAMPLES = 30;
    private GraphTable graph;
    private Random rand;
    private int maxIterations = 100;

    private IteratedLS() {
        rand = new Random();
    }

    public IteratedLS(int seed) {
        this();
        rand.setSeed(seed);
    }

    public void setMaxIterations(int maxIterations) {
        this.maxIterations = maxIterations;
    }

    public void setGraph(GraphTable graph) {
        this.graph = graph;
    }

    public List<Node> getPath(Long start, Long end, float[] weights) {
//        System.out.println(start + " " + end);
        // Dijkstra algorithm helper
        DijkstraWNDim dj = new DijkstraWNDim(RunTLMain.objectives);
        dj.setWeights(weights);
        dj.setGraph(graph);

        List<Node> s0 = dj.getPath(start, end); // initial solution
        List<Node> current, best = new ArrayList<>();
        double fCurrent, fBest;
        current = s0;
//        System.out.println(Arrays.toString(current.toArray()));
        fBest = fitness(current);
        int iter = 0;
        while (iter < maxIterations) {
            current = mutate(current, weights); // weights is a Dijkstra param
            fCurrent = fitness(current);
            if (fCurrent < fBest) {
                fBest = fCurrent;
                best = current;
            }
            iter++;
        }
        return best;
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
        while (j < neighboursList.size() && neighboursList.get(j) == graph.getMapping().get(s.get(point + 1).getId())) {
            j++;
        }
        List<Node> subPath = dj.getPath(neighboursList.get(j), graph.getMapping().get(s.get(s.size()-1).getId()));
        newSolution.addAll(subPath);
        return newSolution;
    }

    private double fitness(List<Node> s) {
        int tlcount = 0;
        long n1, n2;
        double amount = 0;
        List<Double> samples = new ArrayList<>();
        for (int sample = 0; sample < MAX_SAMPLES; sample++) {
            double cost = 0;
            for (int i = 0; i < s.size() - 1; i++) {
                n1 = graph.getMapping().get(s.get(i).getId());
                n2 = graph.getMapping().get(s.get(i+1).getId());
                long arc = graph.getAdjacencyMatrix().get(n1, n2);
                float mu = graph.getWeightsMatrix().get(arc, 0L); // TODO Change to the final weight
                float sigma = graph.getWeightsMatrix().get(arc, 1L); // TODO Change to the final weight
                double tentativeTime = rand.nextGaussian() * sigma + mu;
                // TrafficLight phase
                TlLogic tl = graph.getTlMatrix().get(n1, n2);
                cost += tentativeTime;
                if (tl != null) {
                    tlcount++;
                    int time = tl.calculateTimeStop(Math.round(cost + 0.5d));
                    cost += time;
                }
            }
            amount += cost;
        }
        System.out.println("F " + amount / (double) MAX_SAMPLES + " " + tlcount);
        return amount / (double) MAX_SAMPLES;
    }
}
