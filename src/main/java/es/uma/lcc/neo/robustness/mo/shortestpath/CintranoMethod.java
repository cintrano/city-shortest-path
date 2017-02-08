package es.uma.lcc.neo.robustness.mo.shortestpath;

import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.Dijkstra;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;

import java.util.*;

import static java.lang.Math.min;

/**
 * Created by Christian Cintrano on 31/01/17.
 * Two Phases Method,
 * E.L. Ulungu, J. Teghem, The two phases method: An efficient procedure to solve bi-objective combinatorial
 * optimization problems, Foundations of Computing and Decision Sciences 20 (1995) 149–165.
 */
public class CintranoMethod {

    public static Set<NodePathSolution> compute(GraphTable graph, Long[] objectives) {
        Set<NodePathSolution> sol = new HashSet<NodePathSolution>();
        int objectivesSize = objectives.length;
        float[] boundaries = new float[objectivesSize];

        Set<NodePathSolution> extremes = getLexicographicalExtremes(graph, objectives, boundaries);
        sol.addAll(extremes);

        Stack<NodePathSolution[]> pairs = new Stack<NodePathSolution[]>();
        for (NodePathSolution[] elem : generatePairs(extremes, objectivesSize)) {
            pairs.push(elem);
        }

        // recorrer el stack
        while (!pairs.isEmpty()) {
            NodePathSolution[] par = pairs.pop();
            for (int i = 0; i < objectivesSize; i++) {
                boundaries[i] = min(par[0].getObjectives()[i], par[1].getObjectives()[i]);
            }
            extremes = getLexicographicalExtremes(graph, objectives, boundaries);
            if (!extremes.isEmpty()) {
                sol.addAll(extremes);
                Set<NodePathSolution> extAux = new HashSet<NodePathSolution>();
                extAux.add(par[0]);
                extAux.add(par[1]);
                extAux.addAll(extremes);
                for (NodePathSolution[] elem : generatePairs(extAux, objectivesSize)) {
                    pairs.push(elem);
                }
            }
        }

        return sol;
    }

    private static List<NodePathSolution[]> generatePairs(Set<NodePathSolution> extremes, int size) {
        List<NodePathSolution[]> out = new ArrayList<NodePathSolution[]>();
        List<NodePathSolution> aux = new ArrayList<NodePathSolution>(extremes);

        /*
        TODO
         */

        return out;
    }

    /**
     * Computation of the set Yextr in the bi-objective case (Aneja and Nair (1979))
     * @param graph
     * @param objectives
     * @param boundaries one for each objective
     * @return
     */
    private static Set<NodePathSolution> getLexicographicalExtremes(GraphTable graph, Long[] objectives, float[] boundaries) {
        // Algorithm 1.1
        Set<NodePathSolution> sol = new HashSet<NodePathSolution>();
        /*
        TODO hacer que sean los extremos de verdad
         */
        //DijkstraBounded();
        return sol;
    }

}
