package es.uma.lcc.neo.robustness.mo.shortestpath;

import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;

import java.util.*;

/**
 * Created by Christian Cintrano on 31/01/17.
 * Two Phases Method,
 * E.L. Ulungu, J. Teghem, The two phases method: An efficient procedure to solve bi-objective combinatorial
 * optimization problems, Foundations of Computing and Decision Sciences 20 (1995) 149–165.
 */
public class TwoPhasesMethod {

    public static Set<NodePathSolution> phase1(GraphTable graph, Long[] objectives) {
        Set<NodePathSolution> sol = new HashSet<NodePathSolution>();
        int objectivesSize = objectives.length;
        float[] boundaries = new float[objectivesSize];

        Set<NodePathSolution> extremes = getLexicographicalExtremes(graph, objectives);
        sol.addAll(extremes);

        Stack<NodePathSolution[]> pairs = new Stack<NodePathSolution[]>();
        for (NodePathSolution[] elem : generatePairs(extremes, objectivesSize)) {
            pairs.push(elem);
        }

        return sol;
    }

    private static List<NodePathSolution[]> generatePairs(Set<NodePathSolution> extremes, int size) {
        List<NodePathSolution[]> out = new ArrayList<NodePathSolution[]>();
        List<NodePathSolution> aux = new ArrayList<NodePathSolution>(extremes);


        return out;
    }

    /**
     * Computation of the set Yextr in the bi-objective case (Aneja and Nair (1979))
     * @param graph
     * @param objectives
     * @return
     */
    private static Set<NodePathSolution> getLexicographicalExtremes(GraphTable graph, Long[] objectives) {
        // Algorithm 1.1
        Set<NodePathSolution> sol = new HashSet<NodePathSolution>();
        /*
        TODO hacer que sean los extremos de verdad
         */
        return sol;
    }

}
