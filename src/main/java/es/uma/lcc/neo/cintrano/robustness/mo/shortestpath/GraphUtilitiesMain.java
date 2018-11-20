package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.utilities.ProcessGraph;

/**
 * Created by Christian Cintrano on 24/04/17.
 *
 */
public class GraphUtilitiesMain {

    public static void main (String[] args) {
        System.out.println("=== START EXPERIMENTS ===");
        for (String s : args) {
            System.out.print(s + " ");
        }
        System.out.println();
        //generateMalagaMapFiles();
        //generateMalagaHBEFAWeight();
        //prepareSimpleColoradoWeights();
        prepareSimpleNYWeights();

    }

    private static void prepareSimpleNYWeights() {
        GraphTable graph = ProcessGraph.parserFile("USA-road-d.NY.co");
        graph = ProcessGraph.applyArcs(graph, 11L, "USA-road-d.NY.gr");
        graph = ProcessGraph.applyArcs(graph, 0L, "USA-road-t.NY.gr");

        graph = HBEFAWeightToGraph(graph, 11L, 0L, 1L, 1);

        graph = ProcessGraph.computeRandomWeights(graph, 0L, 0.9f, 1.1f, 2L);
        graph = ProcessGraph.computeRandomWeights(graph, 1L, 0.9f, 1.1f, 3L);

        graph = ProcessGraph.applyMapping(graph, "mapping-ny.txt");

        ProcessGraph.printGraph(graph, "hbefa-ny-graph.xml");
        ProcessGraph.printWeights(graph, "ny_weights_time-hbefa_COMPLETE.xml");

        graph.getWeightsMatrix().column(11L).clear();
        //graph = ProcessGraph.applyMapping(graph, "mapping-malaga.txt");

        ProcessGraph.printWeights(graph, "ny_weights_time-hbefa.xml");

    }

    private static void prepareSimpleColoradoWeights() {
        GraphTable graph = ProcessGraph.parserFile("USA-road-d.COL.co");
        graph = ProcessGraph.applyArcs(graph, 11L, "USA-road-d.COL.gr");
        graph = ProcessGraph.applyArcs(graph, 0L, "USA-road-t.COL.gr");

        graph = HBEFAWeightToGraph(graph, 11L, 0L, 1L, 1);

        graph = ProcessGraph.computeRandomWeights(graph, 0L, 0.9f, 1.1f, 2L);
        graph = ProcessGraph.computeRandomWeights(graph, 1L, 0.9f, 1.1f, 3L);

        graph = ProcessGraph.applyMapping(graph, "mapping-col.txt");

        ProcessGraph.printGraph(graph, "hbefa-col-graph.xml");
        ProcessGraph.printWeights(graph, "col_weights_time-hbefa_COMPLETE.xml");

        graph.getWeightsMatrix().column(11L).clear();
        //graph = ProcessGraph.applyMapping(graph, "mapping-malaga.txt");

        ProcessGraph.printWeights(graph, "col_weights_time-hbefa.xml");

    }

    private static void generateMalagaHBEFAWeight() {
        // Graph
        GraphTable graph = ProcessGraph.parserFile("graph.ampliado.xml"); // Load Speed (km/h) -> 10
        graph = ProcessGraph.getMaxConnectedComponent(graph, "componentes_conexas_malaga.txt");
        //graph = ProcessGraph.applyWeights(graph, weightFilePath0);
        //graph = ProcessGraph.applyWeights(graph, "wVar0.xml");
        //graph = ProcessGraph.applyWeights(graph, "wVar1.xml");

        graph = ProcessGraph.applyWeights(graph, "weights.ampliado.xml"); // Load Distance (m) -> 11
        graph = ProcessGraph.divideWeights(graph, 11L, 10L, 0L, 3.6f); // time (s) -> 0
        graph = HBEFAWeightToGraph(graph, 10L, 1L, 1000/3600);


        graph = ProcessGraph.computeRandomWeights(graph, 0L, 0.9f, 1.1f, 2L);
        graph = ProcessGraph.computeRandomWeights(graph, 1L, 0.9f, 1.1f, 3L);
        ProcessGraph.printGraph(graph, "hbefa-malaga-graph.xml");
        ProcessGraph.printWeights(graph, "weights_time-hbefa_COMPLETE.xml");

        graph.getWeightsMatrix().column(10L).clear();
        graph.getWeightsMatrix().column(11L).clear();
        //graph = ProcessGraph.applyMapping(graph, "mapping-malaga.txt");

        ProcessGraph.printWeights(graph, "weights_time-hbefa.xml");

    }

    private static GraphTable HBEFAWeightToGraph(GraphTable graph, long speed, long res, int k) {
        float v, value;
        for (Long arc : graph.getWeightsMatrix().rowKeySet()) {
            v = graph.getWeightsMatrix().get(arc, speed) * k;
            value = (float) Math.max(0f, 2166.094 + (124.9567 * v) + ( -0.9551173 * v * v) + (0.013322 * v * v * v));
            graph.getWeightsMatrix().put(arc, res, value);
        }
        return graph;
    }

    private static GraphTable HBEFAWeightToGraph(GraphTable graph, long space, long time, long res, int k) {
        float v, value;
        for (Long arc : graph.getWeightsMatrix().rowKeySet()) {
            v = graph.getWeightsMatrix().get(arc, space) * k / graph.getWeightsMatrix().get(arc, time);
            value = (float) Math.max(0f, 2166.094 + (124.9567 * v) + ( -0.9551173 * v * v) + (0.013322 * v * v * v));
            graph.getWeightsMatrix().put(arc, res, value);
        }
        return graph;
    }

    private static void generateMalagaMapFiles() {
        // graph = prepareGraph();
        // graph = ProcessGraph.fixVertexIndex(graph);

        GraphTable graph = ProcessGraph.parserFile("graph.ampliado.xml"); // Load Speed (km/h) -> 10
        graph = ProcessGraph.getMaxConnectedComponent(graph, "componentes_conexas_malaga.txt");
        //graph = ProcessGraph.applyMapping(graph, "mapping-malaga.txt");
        //graph = ProcessGraph.fixVertexIndex(graph);

        graph = ProcessGraph.applyWeights(graph, "weights.ampliado.xml"); // Load Distance (m) -> 11
        graph = ProcessGraph.divideWeights(graph, 11L, 10L, 0L, 3.6f); // time (s) -> 0

        graph = ProcessGraph.addValuesGraph(graph, "malaga.opendata.noise.mod.ssv"); // Load Noise -> 1

        graph = ProcessGraph.computeRandomWeights(graph, 0L, 0.9f, 1.1f, 2L);
        ProcessGraph.printGraph(graph, "new-malaga-graph.xml");
        ProcessGraph.printWeights(graph, "weights_time-noise_COMPLETE.xml");

        graph.getWeightsMatrix().column(10L).clear();
        graph.getWeightsMatrix().column(11L).clear();
        ProcessGraph.printWeights(graph, "weights_time-noise.xml");

/*
        graph = prepareGraph();
        ProcessGraph.printRandomWeights(graph, "wVar0.xml", 0L, 0.9f, 1.1f, 2);
        ProcessGraph.printRandomWeights(graph, "wVar1.xml", 1L, 0.9f, 1.1f, 3);
        ProcessGraph.printMapping(graph);
        */
    }

    private static GraphTable prepareGraph() {
        // Graph
        String graphFilePath = "graph_connected.xml";
        String weightFilePath0 = "wNew.xml";
        GraphTable graph = ProcessGraph.parserFile(graphFilePath);
        System.out.println("Adding weight to the graph...");
        graph = ProcessGraph.applyWeights(graph, weightFilePath0);
        graph = ProcessGraph.applyWeights(graph, "wVar0.xml");
        graph = ProcessGraph.applyWeights(graph, "wVar1.xml");
        graph = ProcessGraph.applyMapping(graph, "mapping-malaga.txt");
        return graph;
    }
}
