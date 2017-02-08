package es.uma.lcc.neo.robustness.mo.shortestpath;

import es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.DijkstraWNDim;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;

import java.util.*;

/**
 * Created by Christian Cintrano on 2/02/17.
 * Simple method to calculate efficient solutions based on normal direction
 */
class NormalMethod {

    private static List<NodePathSolution[]> spaces = new ArrayList<NodePathSolution[]>();

    /*
     TODO: Se considera que el grafo solo tiene los pesos necesarios, ni uno mas, ni uno menos
     */
    static Set<NodePathSolution> compute(GraphTable graph, int objectives, Node start, Node end) {
        Set<NodePathSolution> sol = new HashSet<NodePathSolution>();

        System.out.println("\t__step 1__");
        List<NodePathSolution> aux = new ArrayList<NodePathSolution>();

        Long[] labs = new Long[objectives];
        for (int i = 0; i < labs.length; i++) {
            labs[i] = (long) i;
        }
        System.out.print("\t\tLABs: ");
        for (Long lab : labs) {
            System.out.print(lab + " ");
        }
        System.out.println();


        DijkstraWNDim algorithm = new DijkstraWNDim(labs);
        algorithm.setGraph(graph);

        for (int i = 0; i < objectives; i++) {
            algorithm = new DijkstraWNDim(labs);
            algorithm.setGraph(graph);

            System.out.print("\t\tPi" + i + " ");

            float[] w = new float[objectives];
            w[i] = 1f;

            for (float aW : w) {
                System.out.print(aW + " ");
            }
            System.out.println();

            algorithm.setWeights(w);

            List<Node> path = algorithm.getPath(start, end);
            //for (Node n : path) {
            //    System.out.print(n.getId() + " ");
            //}
            aux.add(new NodePathSolution(graph.getFitness(path), path));
            System.out.print("\t\t\t");
            for (Float f : aux.get(aux.size()-1).getObjectives()) {
                System.out.print(f + " ");
            }
            System.out.println();
        }

        // add solution to output set
        for (NodePathSolution s : aux) {
            sol.add(s);
        }

        System.out.println("\t\tSizes; sol: " + sol.size() + " , aux: " + aux.size());

        System.out.println("\t__step 2__");
        ///// añadir vectores faltantes para crear el plano
        float[][] vectors = new float[objectives][];
        if (sol.size() != objectives) {
            /*
            TODO preparar los vectores necesarios para crear el hiperespacio
             */
            if (sol.size() == 1) {
                return sol;
            }
            aux = repairPoints(sol, aux, objectives);
            int i = aux.size();
            int index = 0;
            while (i > 1) {
                vectors[index] = getVector(aux.get(0).getObjectives(), aux.get(index + 1).getObjectives());
                index++;
                i--;
            }
            /*
            TODO Corregir tanto el aux anterior que no es correcto, como los vectores iguales
             */
        } else {
            // aux has four elements
            vectors[0] = getVector(aux.get(0).getObjectives(), aux.get(1).getObjectives());
            vectors[1] = getVector(aux.get(0).getObjectives(), aux.get(2).getObjectives());
            vectors[2] = getVector(aux.get(0).getObjectives(), aux.get(3).getObjectives());
        }
        float[] normal = getNormal(vectors[0], vectors[1], vectors[2], objectives);


        for (float aNormal : normal) {
            System.out.print(aNormal + " ");
        }
        for (NodePathSolution s : sol) {
            System.out.println(s);
        }
        System.out.println("\t__step 3__");

        algorithm = new DijkstraWNDim(labs);
        algorithm.setGraph(graph);
        algorithm.setWeights(normal);
        System.out.println("\t__step 3.1__");
        List<Node> path = algorithm.getPath(start, end);
        System.out.println("\t__step 3.2__");
        NodePathSolution newSolution = new NodePathSolution(graph.getFitness(path), path);

        System.out.println("\t__step 3.3__");
        sol.add(newSolution);
        Stack<NodePathSolution> stack = new Stack<NodePathSolution>();
        stack.push(newSolution);
        System.out.println("\t__step 3.4__");

        int loopIt = 0;

        while (!stack.isEmpty()) {
            System.out.println("LOOP " + loopIt);
            loopIt++;
            NodePathSolution top = stack.pop();
            if (!sol.contains(top)) {
                // empty the input set
                spaces = new ArrayList<NodePathSolution[]>();
                NodePathSolution[] input = new NodePathSolution[sol.size()];
                int i = 0;
                for (NodePathSolution s : sol) {
                    input[i] = s;
                    i++;
                }
                // get combinations
                combinations(input, objectives - 1, 0, new NodePathSolution[objectives - 1]);
                // get results
                for (NodePathSolution[] space : spaces) {
                    vectors = new float[objectives][];
                    vectors[0] = getVector(top.getObjectives(), space[0].getObjectives());
                    vectors[1] = getVector(top.getObjectives(), space[1].getObjectives());
                    vectors[2] = getVector(top.getObjectives(), space[2].getObjectives());
                    normal = getNormal(vectors[0], vectors[1], vectors[2], objectives);

                    System.out.print("\t\t" + printVector(normal) + " ");

                    algorithm.setWeights(normal);
                    path = algorithm.getPath(start, end);
                    newSolution = new NodePathSolution(graph.getFitness(path), path);

                    System.out.println(printVector(newSolution.getObjectives()));

                    sol.add(newSolution);
                    stack.push(newSolution);
                }
            }
        }


        return sol;
    }

    private static String printVector(float[] vector) {
        String s = "< ";
        for (float f :
                vector) {
            s += f + " ";
        }
        s += ">";
        return s;
    }

    private static List<NodePathSolution> repairPoints(Set<NodePathSolution> sol, List<NodePathSolution> points, int dim) {
        Map<Integer, Integer> indexes = new HashMap<Integer, Integer>();
        for (int i = 0; i < points.size() - 1; i++) {
            for (int j = i + 1; j < points.size(); j++) {
                if (comparePoints(points.get(i), points.get(j))) {
                    indexes.put(i, j);
                }
            }
        }
        List<NodePathSolution> aux = new ArrayList<NodePathSolution>();
        for (NodePathSolution s : sol) {
            aux.add(s);
        }
        for (Integer k : indexes.keySet()) {
            Long[] variables = new Long[0];
            float[] objectives = new float[dim];
            /*
            TODO comprobar cuales tienen que ser los puntos
             */
            objectives[k] = 1;
            objectives[indexes.get(k)] = 1;
            NodePathSolution dummy = new NodePathSolution(objectives, variables);
            aux.add(dummy);
        }
        return aux;
    }

    private static boolean comparePoints(NodePathSolution p1, NodePathSolution p2) {
        int i = 0;
        for (int j = 0; j < p1.getObjectives().length; j++) {
            if (p1.getObjectives()[j] == p2.getObjectives()[j]) {
                i++;
            }
        }
        return i == p1.getObjectives().length;
    }

    private static void combinations(NodePathSolution[] arr, int len, int startPosition, NodePathSolution[] result){
        if (len == 0){
            NodePathSolution[] aux = new NodePathSolution[result.length];
            System.arraycopy(result, 0, aux, 0, result.length);

            spaces.add(aux);
            return;
        }
        for (int i = startPosition; i <= arr.length-len; i++){
            result[result.length - len] = arr[i];
            combinations(arr, len-1, i+1, result);
        }
    }

    private static float[] getVector(float[] p1, float[] p2) {
        float[] v = new float[p1.length];
        for (int i = 0; i < v.length; i++) {
            v[i] = p1[i] + p2[i];
        }
        return v;
    }

    private static float[] getNormal(float[] a, float[] b, float[] c, int dim) {
        /*
        TODO generalizar para cualquier dim. Actualmente está fijo para 4 dimensiones, es decir, tres vectores
         */
        float[] normal = new float[dim];
        normal[0] = (a[1] * b[2] * c[3]) + (a[2] * b[3] * c[1]) + (a[3] * b[1] * c[2]) - (
                (a[3] * b[2] * c[1]) + (a[1] * b[3] * c[2]) + (a[2] * b[1] * c[3]));
        normal[1] = -1 * ((a[0] * b[2] * c[3]) + (a[2] * b[3] * c[0]) + (a[3] * b[0] * c[2]) - (
                (a[3] * b[2] * c[0]) + (a[0] * b[3] * c[2]) + (a[2] * b[0] * c[3])));
        normal[2] = (a[0] * b[1] * c[3]) + (a[1] * b[3] * c[0]) + (a[3] * b[0] * c[1]) - (
                (a[3] * b[1] * c[0]) + (a[0] * b[3] * c[1]) + (a[1] * b[0] * c[3]));
        normal[3] = -1 * ((a[0] * b[1] * c[2]) + (a[1] * b[2] * c[0]) + (a[2] * b[0] * c[1]) - (
                (a[2] * b[1] * c[0]) + (a[0] * b[2] * c[1]) + (a[1] * b[0] * c[2])));


        return makeUnitary(normal);
    }

    private static float[] makeUnitary(float[] vector) {
        float sum = 0f;
        for (float elem : vector) {
            sum += (elem * elem);
        }
        sum = (float) Math.sqrt(sum);

        for (int i = 0; i < vector.length; i++) {
            vector[i] = vector[i] / sum;
        }
        return vector;
    }
}
