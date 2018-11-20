package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.pulse;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.dijkstra.DijkstraWeighted;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.utilities.ProcessGraph;

import java.util.*;

/**
 * Created by Cintrano on 8/02/17.
 * Duque, Lozano & Medaglia, An exact method for the biobjective shortest path problem for large-scale road networks
 */
public class Pulse {


    private GraphTable graph;
    private Long end;
    private Set<NodePathSolution> Xe;

    // Vstart -> V; Map<V, Float[cSup, cInf, tSup, tInf]>
    private Map<Long, float[]> extremesPaths;

    // TODO No se si es correcto o no
    private Map<Long, List<float[]>> L;
    private Long start;

    private Long obj1, obj2;

    public Pulse(Long obj1, Long obj2) {
        extremesPaths = new HashMap<Long, float[]>();
        L = new HashMap<Long, List<float[]>>();
        Xe = new HashSet<NodePathSolution>();
        this.obj1 = obj1;
        this.obj2 = obj2;
    }

    public void setGraph(GraphTable graph) {
        this.graph = graph;
    }

    public List<Node> getPath(Node start, Node end) {
        List<Node> P = new ArrayList<Node>();

        return null;
    }

    public Set<NodePathSolution> pulseAlgorithm(Long start, Long end) {
        this.start = start;
        this.end = end;
        List<Long> P = new ArrayList<Long>();
        float c = 0.0f;
        float t = 0.0f;
        System.out.println("___s0___");
        initialization(graph);
        System.out.println("___s1___");
        pulse(start, c, t, P);
        System.out.println("___s2___" + Xe.size());
        return Xe;
    }

    private void initialization(GraphTable graph) {
        // fill c_, t_
        GraphTable graphInv = inverseGraph(graph);
        System.out.println("___i0___");
        initialization(graphInv, obj1);
        System.out.println("___i1___");
        initialization(graphInv, obj2);
        System.out.println("___i2___");
        // Vstart -> Vend
        float C = getCLimit(graph);
        System.out.println("___i3___" + C);
        float T = getTLimit(graph);
        System.out.println("___i4___" + T);
        // fill c~, t~
        extremesPaths.put(end, new float[]{C, 0.0f, T, 0.0f});
        Map<Long, float[]> newExtremesPaths = new HashMap<Long, float[]>();
        for (Long v : extremesPaths.keySet()) {
            float[] value = extremesPaths.get(v);
            value[0] = C - value[1];
            value[2] = T - value[3];
            newExtremesPaths.put(v, value);
        }
        System.out.println("___i5___");
        extremesPaths = newExtremesPaths;
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
        graphInv.setMapping(graph.getMapping());
        return graphInv;
    }

    private void initialization(GraphTable graph, Long type) {
        boolean[] visited = new boolean[graph.getIntersections().keySet().size()];
        visited[(int) (end - 1)] = true;
        int i = visited[0] ? 1 : 0;

        DijkstraWeighted algorithm = null;
        if (Objects.equals(type, obj1)) {
            algorithm = new DijkstraWeighted(obj1, obj2, 1f);
        }
        if (Objects.equals(type, obj2)) {
            algorithm = new DijkstraWeighted(obj1, obj2, 0f);
        }
        if (algorithm != null) {
            algorithm.setGraph(graph);

//            System.out.print(i + " ");
            while (i < visited.length) {
                List<Node> pathN = algorithm.getPath(end, i + 1L);
                List<Long> path = ProcessGraph.nodeToLong(graph, pathN);
                markSubPaths(path, visited);
                addFitness(graph, path, type);

                i = move(i, visited);
            }
        }
    }


    private float getCLimit(GraphTable graph) {
        DijkstraWeighted algorithm = new DijkstraWeighted(obj1, obj2, (float) 0);
        algorithm.setGraph(graph);
        List<Node> pathN = algorithm.getPath(start, end);
        List<Long> path = ProcessGraph.nodeToLong(graph, pathN);
        float[] fitness = graph.getFitness(path);
        return fitness[obj1.intValue()];
    }
    private float getTLimit(GraphTable graph) {
        DijkstraWeighted algorithm = new DijkstraWeighted(obj1, obj2, (float) 1);
        algorithm.setGraph(graph);
        List<Node> pathN = algorithm.getPath(start, end);
        List<Long> path = ProcessGraph.nodeToLong(graph, pathN);
        float[] fitness = graph.getFitness(path);
        return fitness[obj2.intValue()];
    }
    private float getLimit(GraphTable graph, int i) {
        DijkstraWeighted algorithm = new DijkstraWeighted(obj1, obj2, (float) i);
        algorithm.setGraph(graph);
        List<Node> pathN = algorithm.getPath(start, end);
        List<Long> path = ProcessGraph.nodeToLong(graph, pathN);
        float[] fitness = graph.getFitness(path);
        return fitness[i];
    }
    /*
    private void initialization(GraphTable graph) {
        System.out.println();
        initialization(graph, obj1);
        System.out.println();
        System.out.println("___obj1___");
        initialization(graph, obj2);
        System.out.println();
        System.out.println("___obj2___");
        extremesPaths.put(start, new Float[]{0F,0F,0F,0F});
    }
    private void initialization(GraphTable graph, Long type) {
        boolean[] visited = new boolean[graph.getIntersections().keySet().size()];
        visited[(int) (start - 1)] = true;
        int i = visited[0] ? 1 : 0;

        DijkstraWeighted algorithm = null;
        if (type == obj1) {
            algorithm = new DijkstraWeighted(obj1, obj2, 1f);
        }
        if (type == obj2) {
            algorithm = new DijkstraWeighted(obj1, obj2, 0f);
        }
        if (algorithm != null) {
            algorithm.setGraph(graph);

            while (i < visited.length) {
                List<Node> path = algorithm.getPath(graph.getIntersections().get(start), graph.getIntersections().get(i + obj2));
                markSubPaths(path, visited);
                addFitness(path, type);

                i = move(i, visited);
            }
        }
    }
    */

    private void markSubPaths(List<Long> path, boolean[] visited) {
        for (Long v : path) {
            visited[(int) (v - 1)] = true;
        }
    }

    private void addFitness(GraphTable graph, List<Long> path, long type) {
        float[] fitness;
        float[] sum = new float[]{0.0f, 0.0f};
        //float[] sumInv = new float[2];
        for (int i = 1; i < path.size(); i++) {
            if (extremesPaths.containsKey(path.get(i))) {
                fitness = extremesPaths.get(path.get(i));
            } else {
                fitness = new float[]{0.0f,0.0f,0.0f,0.0f};
            }
            Long arc = graph.getAdjacencyMatrix().get(path.get(i-1), path.get(i));
            //Long arcInv = graph.getAdjacencyMatrix().get(path.get(path.size() - i - 1).getId(), path.get(path.size() - i).getId());
            sum[0] += graph.getWeightsMatrix().get(arc, obj1);
            sum[1] += graph.getWeightsMatrix().get(arc, obj2);
            //sumInv[0] += graph.getWeightsMatrix().get(arcInv, obj1);
            //sumInv[1] += graph.getWeightsMatrix().get(arcInv, obj2);
            if (type == obj1) {
                fitness[1] = sum[0];
                //fitness[1] = sumInv[0];
            }
            if (type == obj2) {
                //fitness[2] = sum[1];
                fitness[3] = sum[1];
            }
            extremesPaths.put(path.get(i), fitness);
        }
    }

    private int move(int i, boolean[] visited) {
        while (i < visited.length && visited[i]) {
            i++;
        }
        return i;
    }

    private void pulse(Long v, float c, float t, List<Long> P) {
        if (isAcyclic(v, P)) {
//            System.out.println("isAcyclic " + v );
//            System.out.println(c + " " + t );
//            System.out.println(extremesPaths.get(v)[0] + " " +  + extremesPaths.get(v)[1] + " " +  + extremesPaths.get(v)[2] + " " +  + extremesPaths.get(v)[3]);
            if (!checkNadirPoint(v, c, t)) {
//                System.out.println("\tcheckNadirPoint " + v);
                if (!checkEfficientSet(v, c, t)) {
//                    System.out.println("\t\tcheckEfficientSet " + v);
                    if (!checkLabels(v, c, t)) {
//                        System.out.println("\t\t\tcheckLabels " + v);
                        //store(c, t);
                        store(v, c, t);
                        List<Long> Pnew = union(P, v);
                        float cP, tP;
                        for (Long vj : outgoingNeighbors(v)) {
                            cP = c + c(v, vj);
                            tP = t + t(v, vj);
                            if (vj.equals(end)) {
//                                System.out.println("pulseend " + vj);
                                pulseEnd(vj, cP, tP, Pnew);
                            } else {
                                pulse(vj, cP, tP, Pnew);
                            }
                        }
                    }
                }
            }
        }
    }

    private void store(Long v, float c, float t) {
        List<float[]> value = L.get(v);
        if (value == null) {
            value = new ArrayList<float[]>();
        }
        value.add(new float[]{c, t});
        L.put(v, value);
    }

    private void pulseEnd(Long v, float c, float t, List<Long> P) {
        if (!checkEfficientSet(v, c, t)) {
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
        while (i.hasNext()) {
            NodePathSolution xP = i.next();
            if ((c(x) <= c(xP)) && (t(x) <= t(xP))) {
                i.remove();
            }
        }
        Xe.add(new NodePathSolution(graph, x));
        System.out.println("Update efficient set. " + Xe.size());
    }

    private float c(Long vi, Long vj) {
        return graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(vi, vj), obj1);
    }

    private float t(Long vi, Long vj) {
        return graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(vi, vj), obj2);
    }

    private Set<Long> outgoingNeighbors(Long v) {
        return graph.getAdjacencyMatrix().row(v).keySet();
    }

    private boolean isAcyclic(Long v, List<Long> p) {
        for (Long node : p) {
            if (v.equals(node)) {
                return false;
            }
        }
        return true;
    }

    private boolean checkNadirPoint(Long v, float c, float t) {
        return (c > (cSup(v)+0.01) || t > (tSup(v)+0.01));
    }

    private float cSup(Long v) {
        return extremesPaths.get(v)[0];
    }
    private float cInf(Long v) {
        return extremesPaths.get(v)[1];
    }
    private float tSup(Long v) {
        return extremesPaths.get(v)[2];
    }
    private float tInf(Long v) {
        return extremesPaths.get(v)[3];
    }


    private boolean checkEfficientSet(Long v, float c, float t) {
        for (NodePathSolution x : Xe) {
            if ((c(x) <= c + cInf(v)) && (t(x) <= t + tInf(v))) {
                return true;
            }
        }
        return false;
    }

    private float c(List<Long> x) {
        return graph.getFitness(x, obj1);
    }

    private float t(List<Long> x) {
        return graph.getFitness(x, obj2);
    }

    private float c(NodePathSolution x) {
        return x.getObjectives()[obj1.intValue()];
    }// valor previo obj1.intValue() = 0

    private float t(NodePathSolution x) {
        return x.getObjectives()[obj2.intValue()];
    }// valor previo obj2.intValue() = 1

    private boolean checkLabels(Long v, float c, float t) {
        /*
        TODO terminal este check
         */
        if (L.get(v) != null) {
            for (float[] costs : L.get(v)) {
                if ((costs[0] <= c) && (costs[1] <= t)) {
                    return true;
                }
            }
        }
        return false;
    }
}
