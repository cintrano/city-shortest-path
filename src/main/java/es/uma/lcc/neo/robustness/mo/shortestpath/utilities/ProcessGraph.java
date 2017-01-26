package es.uma.lcc.neo.robustness.mo.shortestpath.utilities;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.robustness.mo.shortestpath.model.graph.guava.Node;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.*;
import java.util.*;


/**
 * Created by Christian Cintrano on 19/01/17.
 * Support class to handle with graphs
 */
public class ProcessGraph {

    /**
     * Get graph object from a XML file
     * @param file File
     * @return a GraphTable object
     */
    public static GraphTable parserFile(String file) {
        String extension = file.split(".")[file.split(".").length - 1];
        if (extension.equals("xml")) {
            return parserXMLFile(file);
        }
        if (extension.equals("xml")) {
            return parserCOFile(file);
        }
        return null;
    }

    private static GraphTable parserXMLFile(String file) {
        GraphTable graph = new GraphTable();
        try {
            File fXmlFile = new File(file);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(fXmlFile);
            doc.getDocumentElement().normalize();

            NodeList nList = doc.getElementsByTagName("node");
            System.out.println("Load nodes: " + nList.getLength());
            for (int temp = 0; temp < nList.getLength(); temp++) {
                org.w3c.dom.Node nNode = nList.item(temp);
                Element eElement = (Element) nNode;
                graph.getIntersections().put(
                        Long.parseLong(eElement.getAttribute("id")),
                        new Node(
                                Long.parseLong(eElement.getAttribute("id")),
                                eElement.getAttribute("lat"),
                                eElement.getAttribute("lon"))
                );
            }

            nList = doc.getElementsByTagName("arc");
            System.out.println("Load arcs: " + nList.getLength());
            for (int temp = 0; temp < nList.getLength(); temp++) {
                org.w3c.dom.Node nNode = nList.item(temp);
                Element eElement = (Element) nNode;
                graph.getAdjacencyMatrix().put(
                        Long.parseLong(eElement.getAttribute("from")),
                        Long.parseLong(eElement.getAttribute("to")),
                        Long.parseLong(eElement.getAttribute("arcid"))
                );
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return graph;
    }

    private static GraphTable parserCOFile(String file) {
        GraphTable graph = new GraphTable();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
            String line = br.readLine();

            while (line != null) {
                String[] array = line.split(" ");
                if (array[0].equals("v")) {
                    Long id = new Long(array[1]);
                    graph.getIntersections().put(id,
                            new Node(id,
                                    Double.parseDouble(array[3]) / 1000000d,
                                    Double.parseDouble(array[2]) / 1000000d));
                }
                line = br.readLine();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return graph;
    }


    public static GraphTable applyArcs(GraphTable graph, Long type, String file) {
        Long maxArcId = getMaxArcId(graph);
        BufferedReader br = null;
        Long fromId;
        Long toId;
        try {
            br = new BufferedReader(new FileReader(file));
            String line = br.readLine();

            while (line != null) {
                String[] array = line.split(" ");
                if (array[0].equals("a")) {
                    fromId = Long.parseLong(array[1]);
                    toId = Long.parseLong(array[2]);
                    if (graph.getAdjacencyMatrix().get(fromId, toId) == null) {
                        maxArcId++;
                        graph.getAdjacencyMatrix().put(fromId, toId, maxArcId);
                    }
                    graph.getWeightsMatrix().put(graph.getAdjacencyMatrix().get(fromId, toId), type, Float.parseFloat(array[3]));
                }
                line = br.readLine();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return graph;
    }

    private static Long getMaxArcId(GraphTable graph) {
        Long max = 0L;
        if (!graph.getWeightsMatrix().rowKeySet().isEmpty()) {
            List<Long> list = new ArrayList(graph.getWeightsMatrix().rowKeySet());
            Collections.sort(list);
            max = list.get(list.size() - 1);
        }
        return max;
    }

    /**
     * Add the arcs labels to the graph pass to the reference from a XML file
     * @param graph GraphTable object
     * @param file XML file with the arcs labels
     * @return the graph object
     */
    public static GraphTable applyWeights(GraphTable graph, String file) {
        try {
            File fXmlFile = new File(file);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(fXmlFile);
            doc.getDocumentElement().normalize();

            NodeList nList = doc.getElementsByTagName("weight");

            for (int temp = 0; temp < nList.getLength(); temp++) {
                org.w3c.dom.Node nNode = nList.item(temp);
                if (nNode.getNodeType() == org.w3c.dom.Node.ELEMENT_NODE) {
                    Element eElement = (Element) nNode;
                    //System.out.println("weight: " + Long.parseLong(eElement.getAttribute("arcid")) + " " +
                    //         Long.parseLong(eElement.getAttribute("type")) + " " + Float.parseFloat(eElement.getAttribute("value")));
                    graph.getWeightsMatrix().put(
                            Long.parseLong(eElement.getAttribute("arcid")),
                            Long.parseLong(eElement.getAttribute("type")),
                            Float.parseFloat(eElement.getAttribute("value"))
                    );
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return graph;
    }

    private static void getConnectedComponents(GraphTable graph) {
        Map<Long, Integer> ids = new HashMap<Long, Integer>(); // <key, nodeID>
        //Long[] ids = new Long[graph.getIntersections().size()];
        Long[] labels = new Long[graph.getIntersections().size()];
        int i = 0;
        for (Long l : graph.getIntersections().keySet()) {
            ids.put(l, i);
            labels[i] = l;
            i++;
        }
        for (long j = 0; j < graph.getAdjacencyMatrix().size(); j++) {
            Set<Long> toSet = graph.getAdjacencyMatrix().row(j).keySet();
            for (Long l : toSet) {
                changeKey(labels[ids.get(l)], j, labels);
            }
        }

        for (int j = 0; j < labels.length; j++) {
            System.out.println("L: " + labels[j]);
        }

        for (Long node : ids.keySet()) {
            System.out.println("<" + node + " " + ids.get(node) + ">");
        }

        Map<Long, Integer> results = new HashMap<Long, Integer>(); // <nodeId, numElementInComponent>
        for (int j = 0; j < labels.length; j++) {
            if (results.get(labels[j]) == null) {
                results.put(labels[j], 1);
            } else {
                results.put(labels[j], results.get(labels[j]) + 1);
            }
        }
        System.out.println("Connected Components");
        for (Long node : results.keySet()) {
            System.out.println("<" + node + " " + results.get(node) + ">");
        }
    }

    /**
     *
     * @param oldLabel
     * @param newLabel
     * @param labels
     */
    private static void changeKey(Long oldLabel, Long newLabel, Long[] labels) {
        System.out.print("Changed");
        for (int i = 0; i < labels.length; i++) {
            if (labels[i].equals(oldLabel)) {
                labels[i] = newLabel;
            }
        }
    }


    // The main function that finds and prints all strongly
    // connected components
    public static void printSCCs(GraphTable graph) {
        Map<Integer, Long> ids = new HashMap<Integer, Long>(); // <key, nodeID>
        Map<Long, Integer> nodes = new HashMap<Long, Integer>(); // <key, nodeID>
        //Long[] ids = new Long[graph.getIntersections().size()];
        Long[] labels = new Long[graph.getIntersections().size()];
        int i = 0;
        for (Long l : graph.getIntersections().keySet()) {
            ids.put(i, l);
            nodes.put(l, i);
            labels[i] = l;
            i++;
        }

        Stack stack = new Stack();

        // Mark all the vertices as not visited (For first DFS)
        boolean visited[] = new boolean[graph.getIntersections().size()];
        for(i = 0; i < graph.getIntersections().size(); i++)
            visited[i] = false;

        // Fill vertices in stack according to their finishing
        // times
        for (i = 0; i < visited.length; i++)
            if (visited[i] == false)
                fillOrder(i, visited, stack, graph, ids, nodes, true);

        // Create a reversed graph
        //Graph gr = getTranspose();

        // Mark all the vertices as not visited (For second DFS)
        for (i = 0; i < graph.getIntersections().size(); i++)
            visited[i] = false;

        // Now process all vertices in order defined by Stack
        while (stack.empty() == false) {
            // Pop a vertex from stack
            int v = (Integer) stack.pop();

            // Print Strongly connected component of the popped vertex
            if (visited[v] == false) {
                DFSUtil(v, visited, graph, ids, nodes, false);
                System.out.println();
            }
        }
    }

    // A recursive function to print DFS starting from v
    static void DFSUtil(int v, boolean visited[], final GraphTable graph, final Map<Integer, Long> ids,
                        final Map<Long, Integer> nodes, boolean order) {
        // Mark the current node as visited and print it
        visited[v] = true;
        System.out.print(ids.get(v) + " ");

        //int n;

        // Recur for all the vertices adjacent to this vertex
        if (order) {
            graph.getAdjacencyMatrix().row(ids.get(v));
            for (Long n : graph.getAdjacencyMatrix().row(ids.get(v)).keySet()) {
                if (nodes.get(n)!= null && !visited[nodes.get(n)])
                    DFSUtil(nodes.get(n), visited, graph, ids, nodes, order);
            }
        } else {
            graph.getAdjacencyMatrix().column(ids.get(v));
            for (Long n : graph.getAdjacencyMatrix().column(ids.get(v)).keySet()) {
                if (nodes.get(n)!= null && !visited[nodes.get(n)])
                    DFSUtil(nodes.get(n), visited, graph, ids, nodes, order);
            }
        }
        /*
        Iterator<Integer> i =adj[v].iterator();
        while (i.hasNext())
        {
            n = i.next();
            if (!visited[n])
                DFSUtil(n,visited);
        }
        */
    }


    static void fillOrder(int v, boolean visited[], Stack stack, final GraphTable graph, final Map<Integer,
            Long> ids, final Map<Long, Integer> nodes, boolean order) {
        // Mark the current node as visited and print it
        visited[v] = true;

        // Recur for all the vertices adjacent to this vertex
        if (order) {
            graph.getAdjacencyMatrix().row(ids.get(v));
            for (Long n : graph.getAdjacencyMatrix().row(ids.get(v)).keySet()) {
                //System.out.println(v + "->" + nodes.get(n) + "..." + ids.get(v) + "->" + n); // idArray ... idNodes
                //for (Long aux :
                //        graph.getAdjacencyMatrix().row(ids.get(v)).keySet()) {
                //    System.out.print(aux + " ");
                //}
                //System.out.println();
                if (nodes.get(n)!= null && !visited[nodes.get(n).intValue()])
                    fillOrder(nodes.get(n), visited, stack, graph, ids, nodes, order);
            }
        } else {
            graph.getAdjacencyMatrix().column(ids.get(v));
            for (Long n : graph.getAdjacencyMatrix().column(ids.get(v)).keySet()) {
                if (nodes.get(n)!= null && !visited[nodes.get(n).intValue()])
                    fillOrder(nodes.get(n), visited, stack, graph, ids, nodes, order);
            }
        }
        /*
        Iterator<Integer> i = adj[v].iterator();
        while (i.hasNext())
        {
            int n = i.next();
            if(!visited[n])
                fillOrder(n, visited, stack);
        }
*/
        // All vertices reachable from v are processed by now,
        // push v to Stack
        stack.push(v);
    }



    public static GraphTable getMaxConnectedComponent(GraphTable graph, String filepath) {
        System.out.println("=== making graph strong connected ===");
        String[] bestCC = new String[0];
        String[] auxCC;

        String line;
        File file = new File(filepath);
        FileInputStream fis = null;
        BufferedReader br = null;
        InputStreamReader isr = null;
        try {
            fis = new FileInputStream(file);
            isr = new InputStreamReader(fis);
            br = new BufferedReader(isr);
            while ((line = br.readLine()) != null) {
                auxCC = line.split(" ");
                if (auxCC.length > bestCC.length) {
                    bestCC = auxCC;
                }
            }

            List<Long> convexComponent = new ArrayList<Long>();
            for (int i = 0; i < bestCC.length; i++) {
                convexComponent.add(Long.parseLong(bestCC[i]));
            }

            // remove node
            Iterator iter = graph.getIntersections().keySet().iterator();
            while (iter.hasNext()) {
                if (!convexComponent.contains(iter.next())) {
                    iter.remove();
                }
            }
            Table<Long, Long, Long> adjacencyMatrix = HashBasedTable.create();
            // remove arcs
            for (Long r : graph.getAdjacencyMatrix().rowKeySet()) {
                if (convexComponent.contains(r)) {
                    for (Long c : graph.getAdjacencyMatrix().row(r).keySet()) {
                        if (convexComponent.contains(c)) {
                            adjacencyMatrix.put(r, c, graph.getAdjacencyMatrix().get(r, c));
                        }
                    }
                }
            }
            graph.setAdjacencyMatrix(adjacencyMatrix);

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null) br.close();
                if (isr != null) isr.close();
                if (fis != null) fis.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        System.out.println("=== return graph strong connected ===");
        return graph;
    }


    /**
     * Make a file with the data in a GraphTable object
     * @param graph GraphTable object
     * @param filepath File output name
     */
    public static void printGraph(GraphTable graph, String filepath) {
        System.out.println("writing file...");
        Writer writer = null;

        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filepath), "utf-8"));

            writer.write("<graph>\n");
            String line;
            // arcs
            writer.write("<arcs>\n");
            for (Long r : graph.getAdjacencyMatrix().rowKeySet()) {
                for (Long c : graph.getAdjacencyMatrix().row(r).keySet()) {
                    line = "<arc " +
                            "arcid='" + graph.getAdjacencyMatrix().get(r, c) +
                            "' from='" + r +
                            "' to='" + c +
                            "' />\n";

                    writer.write(line);
                }
            }
            writer.write("</arcs>\n");
            // nodes
            writer.write("<nodes>\n");
            for (Long n : graph.getIntersections().keySet()) {
                System.out.println(graph.getIntersections().get(n));
                line = "<node " +
                        "id='" + n +
                        "' lat='" + graph.getIntersections().get(n).getLatitude().toString() +
                        "' lon='" + graph.getIntersections().get(n).getLongitude().toString() +
                        "' />\n";

                writer.write(line);
            }
            writer.write("</nodes>\n");


            writer.write("</graph>\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {writer.close();} catch (Exception ex) {ex.printStackTrace();}
        }

        System.out.println("...wrote file");
    }

    public static void printWeights(GraphTable graph, String filepath) {
        System.out.println("writing file " + filepath + "...");
        Writer writer = null;

        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filepath), "utf-8"));

            String line;
            // arcs
            writer.write("<weights>\n");
            for (Long r : graph.getWeightsMatrix().rowKeySet()) {
                for (Long c : graph.getWeightsMatrix().row(r).keySet()) {
                    line = "<weight " +
                            "arcid='" + r +
                            "' value='" + graph.getWeightsMatrix().get(r, c) +
                            "' type='" + c +
                            "' />\n";

                    writer.write(line);
                }
            }
            writer.write("</weights>\n");

            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {
                assert writer != null;
                writer.close();} catch (Exception ex) {ex.printStackTrace();}
        }

        System.out.println("...wrote file");
    }


    /**
     * Make a random labels weights file to a graph
     * @param graph GraphTable object, it is not necessary to have previous weights
     * @param filepath File output name
     * @param lowerBound lowerBound
     * @param upperBound upperBound
     * @param type type attribute for the label
     */
    public static void printRandomWeights(GraphTable graph, String filepath, float lowerBound, float upperBound, int type) {
        Random random = new Random(type);
        System.out.println("writing file " + filepath + "...");
        Writer writer = null;

        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filepath), "utf-8"));

            String line;
            // arcs
            writer.write("<weights>\n");
            float value;
            for (Long arcid : graph.getWeightsMatrix().rowKeySet()) {
                value = ((( upperBound - lowerBound) * random.nextFloat()) + lowerBound) * graph.getWeightsMatrix().get(arcid, 4L);
                line = "<weight " +
                        "arcid='" + arcid +
                        "' value='" + value +
                        "' type='" + type +
                        "' />\n";

                writer.write(line);
            }
            writer.write("</weights>\n");

            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {
                assert writer != null;
                writer.close();} catch (Exception ex) {ex.printStackTrace();}
        }

        System.out.println("...wrote file");
    }

    public static void printRandomWeights(GraphTable graph, String filepath, Long baseType, float lowerBound, float upperBound, int type) {
        Random random = new Random(type);
        System.out.println("writing file " + filepath + "...");
        Writer writer = null;

        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filepath), "utf-8"));

            String line;
            // arcs
            writer.write("<weights>\n");
            float value;
            for (Long arcid : graph.getWeightsMatrix().rowKeySet()) {
                value = ((( upperBound - lowerBound) * random.nextFloat()) + lowerBound) * graph.getWeightsMatrix().get(arcid, baseType);
                line = "<weight " +
                        "arcid='" + arcid +
                        "' value='" + value +
                        "' type='" + type +
                        "' />\n";

                writer.write(line);
            }
            writer.write("</weights>\n");

            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {
                assert writer != null;
                writer.close();} catch (Exception ex) {ex.printStackTrace();}
        }

        System.out.println("...wrote file");
    }

    public static void computeBounds(GraphTable graph) {
        Float value;
        for (Long type : graph.getWeightsMatrix().columnKeySet()) {
            Float min = Float.MAX_VALUE;
            Float max = Float.MIN_VALUE;
            for (Long arc : graph.getWeightsMatrix().column(type).keySet()) {
                value = graph.getWeightsMatrix().get(arc, type);
                if (value < min) {
                    min = value;
                }
                if (value > max) {
                    max = value;
                }
            }
            graph.getLowerBound().put(type, min);
            graph.getUpperBound().put(type, max);
        }
    }

    public static GraphTable normalizate(GraphTable graph) {
        Float value;
        Float min = Float.MAX_VALUE;
        for (Long type : graph.getWeightsMatrix().columnKeySet()) {
            System.out.println("analizing type: " + type);
            Float avg = 0F;
            Float sd = 0F;
            Float sumSquare = 0F;
            Float n = 0F;
            for (Long arc : graph.getWeightsMatrix().column(type).keySet()) {
                value = graph.getWeightsMatrix().get(arc, type);
                avg += value;
                sumSquare += value * value;
                n++;
            }
            avg = avg / n;
            sumSquare = sumSquare / n;
            sd = new Float(Math.sqrt(sumSquare.doubleValue() - (avg.doubleValue() * avg.doubleValue())));
            System.out.println("normalizing: " + type);
            for (Long arc : graph.getWeightsMatrix().column(type).keySet()) {
                graph.getWeightsMatrix().put(arc, type, (graph.getWeightsMatrix().get(arc, type) - avg) / sd);
            }

        }

        for (Long r : graph.getWeightsMatrix().rowKeySet()) {
            for (Long c : graph.getWeightsMatrix().row(r).keySet()) {
                if (min > graph.getWeightsMatrix().get(r, c)) {
                    min = graph.getWeightsMatrix().get(r, c);
                }
            }
        }

        System.out.println("moving distributions");
        min = -1 * min;//min;
        for (Long r : graph.getWeightsMatrix().rowKeySet()) {
            for (Long c : graph.getWeightsMatrix().row(r).keySet()) {
                graph.getWeightsMatrix().put(r, c, graph.getWeightsMatrix().get(r, c) + min);
            }
        }
        System.out.println("normalization process ENS");
        return graph;
    }

    public static GraphTable computeNewWeight(GraphTable graph, Long type, Long inverseType, Long w1, Long w2) {
        // w1 / w2
        Float value;
        for (Long arc : graph.getWeightsMatrix().rowKeySet()) {
            value = graph.getWeightsMatrix().get(arc, w1) / graph.getWeightsMatrix().get(arc, w2);
            graph.getWeightsMatrix().put(arc, type, -1 * value);
            graph.getWeightsMatrix().put(arc, inverseType, value);
        }
        return graph;
    }
}
