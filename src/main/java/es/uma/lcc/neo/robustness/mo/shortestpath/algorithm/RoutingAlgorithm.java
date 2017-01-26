package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm;

import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;

import java.util.List;

/**
 * Created by Christian Cintrano on 26/01/17.
 * Routing Algorithm interface
 */
public interface RoutingAlgorithm {

    public void setGraph(GraphTable graph);

    public List<Node> getPath(Node start, Node end);
}
