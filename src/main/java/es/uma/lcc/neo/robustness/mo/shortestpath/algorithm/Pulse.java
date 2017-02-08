package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm;

import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Created by Cintrano on 8/02/17.
 * Duque, Lozano & Medaglia, An exact method for the biobjective shortest path problem for large-scale road networks
 */
public class Pulse implements RoutingAlgorithm {

    private GraphTable graph;
    public void setGraph(GraphTable graph) {
        this.graph = graph;
    }

    public List<Node> getPath(Node start, Node end) {
        List<Node> P = new ArrayList<Node>();

        return null;
    }

    public Set pulseAlgorithm(Node start, Node end) {
        List<Node> P = new ArrayList<Node>();
        float c = 0;
        float t = 0;
        initialization(graph);
        pulse(start, c, t, P);
        return Xe;
    }

    private void pulse(Node v, float c, float t, List<Node> P) {
        if (isAcyclic(v, P)) {
            if (!checkNadirPoint(v, c, t)) {
                if (!checkEfficientSet(v, c, t)) {
                    if (!checkLabels(v, c, t)) {
                        store(c, t);
                        List<Node> Pnew = union(P, v);
///////                        for (: ) {
                                c = c +;
                                t = t + ;
                                pulse(vj, c, t, Pnew);
///////                        }
                    }
                }
            }
        }
    }

    private boolean checkNadirPoint(Node v, float c, float t) {
        boolean prune = false;
        if (c > cSup(v) || t > tSup(v)) { //~c
            prune = true;
        }
        return prune;
    }


    private boolean checkEfficientSet(Node v, float c, float t) {
        boolean prune = false;
        for (Node x : efficientSet) {
            if (c(x) <= c + cInf(v) && t(x) <= t + tInf(v)) {
                prune = true;
            }
        }
        return prune;
    }

    private boolean checkLabels(Node v, float c, float t) {
        boolean prune = false;
        
        return prune;
    }

}
