package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.pulse;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraSimpleWNDim;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;

import java.util.*;

/**
 * Created by Cintrano on 8/02/17.
 * Duque, Lozano & Medaglia, An exact method for the biobjective shortest path problem for large-scale road networks
 */
public class PulseMO {

    private GraphTable graph;
    private Long end;
    private Set<NodePathSolution> Xe;

    // Vstart -> V; Map<V, Float[cSup, cInf, tSup, tInf]>
    private int objectivesNumber;
    private Map<Long, Float[]> cInf; // array.size = objectivesNumber
    private float[] zn; // |zn| = objectivesNumber
    //private Map<Long, Float[]> extremesPaths;

    // TODO No se si es correcto o no
    private Map<Long, List<Float[]>> L;
    private Long start;


    final static private Long[] LABS = new Long[]{0L, 1L, 2L, 3L};

    public PulseMO(int objectivesNumber) {
        this.objectivesNumber = objectivesNumber;
        zn = new float[objectivesNumber];
        cInf = new HashMap<Long, Float[]>();
        //extremesPaths = new HashMap<Long, Float[]>();
        L = new HashMap<Long, List<Float[]>>();
        Xe = new HashSet<NodePathSolution>();
    }

    public void setGraph(GraphTable graph) {
        this.graph = graph;
    }

    public Set<NodePathSolution> pulseAlgorithm(Long start, Long end) {
        this.start = start;
        this.end = end;
        List<Long> P = new ArrayList<Long>();
        System.out.println("___s0___");
        initialization(graph);
        System.out.println("___s1___");
        float[] c = new float[objectivesNumber];
        pulse(start, c, P);
        System.out.println("___s2___");
        return Xe;
    }

    private void initialization(GraphTable graph) {
        // fill c_, t_
        GraphTable graphInv = inverseGraph(graph);
        for (Long lab : LABS) {
            System.out.println("LAB " + lab);
            initialization(graphInv, lab);
        }
        System.out.println("___i2___");
        pseudoNadirPoint(graph);
        System.out.println("___i5___");
    }

    private GraphTable inverseGraph(GraphTable graph) {
        GraphTable graphInv = new GraphTable();
        graphInv.setIntersections(graph.getIntersections());
        graphInv.setWeightsMatrix(graph.getWeightsMatrix());
        Table<Long, Long, Long> adjacencyMatrix = HashBasedTable.create();
        for (Long r : graph.getAdjacencyMatrix().rowKeySet()) {
            for (Long c : graph.getAdjacencyMatrix().row(r).keySet()) {
                adjacencyMatrix.put(c, r, graph.getAdjacencyMatrix().get(r, c));
            }
        }
        graphInv.setAdjacencyMatrix(adjacencyMatrix);
        return graphInv;
    }

    private void initialization(GraphTable graph, Long type) {
        boolean[] visited = new boolean[graph.getIntersections().keySet().size()];
        visited[(int) (end - 1)] = true;
        int i = visited[0] ? 1 : 0;

        DijkstraSimpleWNDim algorithm = new DijkstraSimpleWNDim(LABS);
        float[] weight = new float[objectivesNumber];
        weight[type.intValue()] = 1;
        algorithm.setWeights(weight);
        algorithm.setGraph(graph);

        while (i < visited.length) {
            List<Long> path = algorithm.getPath(end, i + 1L);
            markSubPaths(path, visited);
            addFitness(graph, path, type);

            i = move(i, visited);
        }
        cInf.put(end, new Float[]{0F, 0F, 0F, 0F});
    }

    private void pseudoNadirPoint(GraphTable graph) {
        for (Long lab : LABS) {
            DijkstraSimpleWNDim algorithm = new DijkstraSimpleWNDim(LABS);
            float[] weight = new float[objectivesNumber];
            weight[lab.intValue()] = 1;
            algorithm.setGraph(graph);
            algorithm.setWeights(weight);
            List<Long> path = algorithm.getPath(start, end);
            float[] fitness = graph.getFitness(path);
            for (int i = 0; i < fitness.length; i++) {
                if (zn[i] < fitness[i]) {
                    zn[i] = fitness[i];
                }
            }
        }
    }

    private void markSubPaths(List<Long> path, boolean[] visited) {
        for (Long v : path) {
            visited[(int) (v - 1)] = true;
        }
    }

    private void addFitness(GraphTable graph, List<Long> path, Long type) {
        Float[] fitness;
        float sum = 0f;
        //float[] sumInv = new float[2];
        for (int i = 1; i < path.size(); i++) {
            if (cInf.containsKey(path.get(i))) {
                fitness = cInf.get(path.get(i));
            } else {
                fitness = new Float[objectivesNumber];
                for (int j = 0; j < fitness.length; j++) {
                    fitness[j] = 0F;
                }
            }
            Long arc = graph.getAdjacencyMatrix().get(path.get(i-1), path.get(i));
            sum += graph.getWeightsMatrix().get(arc, type);
            fitness[type.intValue()] = sum;
            cInf.put(path.get(i), fitness);
        }
    }

    private int move(int i, boolean[] visited) {
        while (i < visited.length && visited[i]) {
            i++;
        }
        return i;
    }

    private void pulse(Long v, float[] c, List<Long> P) {
        if (isAcyclic(v, P)) { //adapted
            if (!checkNadirPoint(v, c)) {
                if (!checkEfficientSet(v, c)) {
                    if (!checkLabels(v, c)) {
                        //store(c, t);
                        store(v, c);
                        List<Long> Pnew = union(P, v);
                        float[] cP = new float[c.length];
                        for (Long vj : outgoingNeighbors(v)) {
                            for (int i = 0; i < c.length; i++) {
                                cP[i] = c[i] + c(v, vj, (long) i);
                            }
                            if (vj.equals(end)) {
                                pulseEnd(vj, cP, Pnew);
                            } else {
                                pulse(vj, cP, Pnew);
                            }
                        }
                    }
                }
            }
        }
    }

    private void store(Long v, float c[]) {
        List<Float[]> value = L.get(v);
        if (value == null) {
            value = new ArrayList<Float[]>();
        }
        Float[] aux = new Float[c.length];
        for (int i = 0; i < c.length; i++) {
            aux[i] = c[i];
        }
        value.add(aux);
        L.put(v, value);
    }

    private void pulseEnd(Long v, float[] c, List<Long> P) {
        if (!checkEfficientSet(v, c)) {
            List<Long> Pnew = union(P, v);
            //x = mapPathToSolution(P);
            //updateEfficientSet(x);
            // TODO cambiar por NodePathSolutions el x del mapPathToSolution(P)
            updateEfficientSet(Pnew);
        }
    }

    private List<Long> union(List<Long> p, Long v) {
        List<Long> path = new ArrayList<Long>(p);
        path.add(v);
        return path;
    }

    private void updateEfficientSet(List<Long> x) {
        Iterator<NodePathSolution> i = Xe.iterator();
        NodePathSolution solutionNew = new NodePathSolution(graph, x);
        while (i.hasNext()) {
            NodePathSolution xP = i.next();
            int sum = 0;
            for (int j = 0; j < xP.getObjectives().length; j++) {
                if (solutionNew.getObjectives()[j] <= xP.getObjectives()[j]) {
                    sum++;
                }
            }
            if (sum == xP.getObjectives().length) {
                i.remove();
            }
        }
        Xe.add(solutionNew);
    }

    private float c(Long vi, Long vj, Long weight) {
        return graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(vi, vj), weight);
    }

    private Set<Long> outgoingNeighbors(Long v) {
        return graph.getAdjacencyMatrix().row(v).keySet();
    }

    private boolean isAcyclic(Long v, List<Long> p) {
        /*
        for (Long node : p) {
            if (v.equals(node)) {
                return false;
            }
        }
        return true;
        */
        return !p.contains(v);
    }

    private boolean checkNadirPoint(Long v, float[] c) {
        Float[] cost = cInf.get(v);
        for (int i = 0; i < c.length; i++) {
            //if (c[i] > cSup(v, i)) {
            if (c[i] + cost[i] > zn[i]) {
                return true;
            }
        }
        return false;
    }

    private float cInf(Long v, int i) {
        return cInf.get(v)[i];
    }

    private boolean checkEfficientSet(Long v, float[] c) {
        for (NodePathSolution x : Xe) {
            int sum = 0;
            for (int i = 0; i < c.length; i++) {
                if (c(x, i) <= c[i] + cInf(v, i)) {
                    sum++;
                }
            }
            if (sum == c.length) {
                return true;
            }
        }
        return false;
    }

    private float c(NodePathSolution x, int i) {
        return x.getObjectives()[i];
    }

    private boolean checkLabels(Long v, float[] c) {
        /*
        TODO terminal este check
         */
        if (L.get(v) != null) {
            for (Float[] costs : L.get(v)) {
                int sum = 0;
                for (int i = 0; i < c.length; i++) {
                    if (costs[i] <= c[i]) {
                        sum++;
                    }
                }
                if (sum == c.length) {
                    return true;
                }
            }
        }
        return false;
    }
}
