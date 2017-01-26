package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm;

import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;

import java.util.*;


/**
 * Created by Christian Cintrano on 16/01/17.
 * A* algorithm
 */
public class Astar implements RoutingAlgorithm {

    private GraphTable graph;
    private Long target;

    public void setGraph(GraphTable graph) {
        this.graph = graph;
        //fillG();
        //fillH();
    }
    public void setTarget(Long target) {
        this.target = target;
    }

    public LinkedList<Node> getPath(Node from, Node to) {
        //System.out.println("A* getting path from:" + from + " to " + to);
        Set<Long> closedSet = new HashSet<Long>();

        Set<Long> openSet = new HashSet<Long>();
        openSet.add(from.getId());

        Map<Long, Long> cameFrom = new HashMap<Long, Long>();

        Map<Long, Float> gScore = new HashMap<Long, Float>();
        for (Long node : graph.getAdjacencyMatrix().rowKeySet()) {//graph.getIntersections().keySet()) {
            gScore.put(node, Float.MAX_VALUE);
        }

        gScore.put(from.getId(), 0F);

        Map<Long, Float> fScore = new HashMap<Long, Float>();
        for (Long node : graph.getAdjacencyMatrix().rowKeySet()) {////graph.getIntersections().keySet()) {
            fScore.put(node, Float.MAX_VALUE);
            ////System.out.println("M  " + fScore.get(node));
        }

        fScore.put(from.getId(), heuristicCostEstimate(from.getId(), to));

        Long current = null;
        Float currentFValue = Float.MAX_VALUE;
        while (!openSet.isEmpty()) {
            current = null;
            currentFValue = Float.MAX_VALUE;
            for (Long node : openSet) {
                //System.out.println("NODE  " + node + " " + currentFValue + " " + fScore.containsKey(node) + " " + fScore.get(node));
                if (currentFValue >= fScore.get(node)) {
                    current = node;
                    currentFValue = fScore.get(node);
                    //System.out.println("current");
                }
            }
            //System.out.println(current);

            if (current.equals(to.getId())) {
                return reconstructPath(cameFrom, current);
            }

            openSet.remove(current);

            closedSet.add(current);
            Float tentativeGScore;
            for (Long neighbor: graph.getAdjacencyMatrix().row(current).keySet()) {
                //System.out.println("Neighbor " + neighbor);
                if (!closedSet.contains(neighbor)) { // changed from Wikipedia because is a road shortest path problem
                    //System.out.println("contains");
                    // break this loop it

                    tentativeGScore = gScore.get(current) + distBetween(current, neighbor);
                    boolean control = true;
                    if (!openSet.contains(neighbor)) {
                        //System.out.println("!openSet.contains(neighbor) " + neighbor);
                        openSet.add(neighbor);
                    } else if (tentativeGScore >= gScore.get(neighbor)) {
                        // break this loop it
                        control = false;
//System.out.println("neighbor " + neighbor);
                    }
                    if(control) {
                        cameFrom.put(neighbor, current);
                        gScore.put(neighbor, tentativeGScore);
                        fScore.put(neighbor, gScore.get(neighbor) + heuristicCostEstimate(neighbor, to));
                    }
                }
            }
        }
        //System.out.println("No path found");
        return new LinkedList<Node>();


        /*
        openList = new LinkedList<Node>();
        closedList = new LinkedList<Node>();
        openList.add(from); // add starting node to open list

        boolean done = false;
        Node current;
        while (!done) {
            current = lowestFInOpen(); // get node with lowest fCosts from openList
            closedList.add(current); // add current node to closed list
            openList.remove(current); // delete current node from open list

            if (current.equals(to)) { // found goal
                return calcPath(from, current);
            }

            // for all adjacent nodes:
            //List<Node> adjacentNodes = getAdjacent(current);
            List<Long> adjacentNodes = new ArrayList<Long>(graph.getAdjacencyMatrix().row(current.getId()).keySet());
            for (int i = 0; i < adjacentNodes.size(); i++) {
                Long currentAdj = adjacentNodes.get(i);
                if (!openList.contains(currentAdj)) { // node is not in openList
                    currentAdj.setPrevious(current); // set current node as previous for this node
                    currentAdj.sethCosts(to); // set h costs of this node (estimated costs to goal)
                    currentAdj.setgCosts(current); // set g costs of this node (costs from start to this node)
                    openList.add(currentAdj); // add node to openList
                } else { // node is in openList
                    if (g.get(currentAdj) > currentAdj.calculategCosts(current)) { // costs from current node are cheaper than previous costs
                        currentAdj.setPrevious(current); // set current node as previous for this node
                        currentAdj.setgCosts(current); // set g costs of this node (costs from start to this node)
                    }
                }
            }

            if (openList.isEmpty()) { // no path exists
                return new LinkedList<Node>(); // return empty list
            }
        }
        */
    }

    private Float distBetween(Long current, Long neighbor) {
        return (graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(current, neighbor), 0L) * 0.5F) +
                (graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(current, neighbor), 1L) * 0.5F);
    }

    private LinkedList<Node> reconstructPath(Map<Long, Long> cameFrom, Long current) {
        LinkedList<Node> totalPath = new LinkedList<Node>();
        totalPath.add(graph.getIntersections().get(current));
//System.out.print("current ");
        while (cameFrom.keySet().contains(current)) {
//            System.out.print(current + " ");
            current = cameFrom.get(current);

            totalPath.add(graph.getIntersections().get(current));
        }
        LinkedList<Node> output = new LinkedList<Node>();
        for (int i = (totalPath.size() - 1); i >= 0; i--) {
            output.add(totalPath.get(i));
        }
//        System.out.println();
        return output;
    }

    private Float heuristicCostEstimate(Long from, Node to) {
        Node fromNode = graph.getIntersections().get(from);
        Double distance = Math.sqrt((to.getLatitude().doubleValue() - fromNode.getLatitude().doubleValue()) * (to.getLatitude().doubleValue() - fromNode.getLatitude().doubleValue()) +
                (to.getLongitude().doubleValue() - fromNode.getLongitude().doubleValue()) * (to.getLongitude().doubleValue() - fromNode.getLongitude().doubleValue()));
        Double speed = findSpeed(graph, from);
        distance = distance / speed;//graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(from, to.getId()), 5L);
        return distance.floatValue();
    }


    private Double findSpeed(GraphTable graph, Long from) {
        //JMetalRandom random = JMetalRandom.getInstance();
        List<Long> aux = new ArrayList(graph.getAdjacencyMatrix().row(from).keySet());
        int index = 0;//random.nextInt(0, aux.size()-1);
        Long to = aux.get(index);
        return  new Double(graph.getWeightsMatrix().get(graph.getAdjacencyMatrix().get(from, to), 5L));
    }


}
