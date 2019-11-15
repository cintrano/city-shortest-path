package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.utilities;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.TlLogic;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.*;
import java.math.BigDecimal;
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
        String extension = file.split("\\.")[file.split("\\.").length - 1];
        if (extension.equals("xml")) {
            return parserXMLFile(file);
        }
        if (extension.equals("co")) {
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
            int size = nList.getLength();
            for (int temp = 0; temp < size; temp++) {
                org.w3c.dom.Node nNode = nList.item(temp);
                Element eElement = (Element) nNode;
                graph.getIntersections().put(
                        Long.parseLong(eElement.getAttribute("id")),
                        new Node(
                                Long.parseLong(eElement.getAttribute("id")),
                                eElement.getAttribute("lat"),
                                eElement.getAttribute("lon"))
                );
                graph.getInverseIntersections().put(graph.getIntersections().get(Long.parseLong(eElement.getAttribute("id"))),
                        Long.parseLong(eElement.getAttribute("id"))
                );
            }

            nList = doc.getElementsByTagName("arc");
            System.out.println("Load arcs: " + nList.getLength());
            //org.w3c.dom.Node nNode;
            Element eElement;
            size = nList.getLength();
            for (int temp = 0; temp < size; temp++) {
                //nNode = nList.item(temp);

                eElement = (Element) nList.item(temp);//nNode;

                graph.getAdjacencyMatrix().put(
                        Long.parseLong(eElement.getAttribute("from")),
                        Long.parseLong(eElement.getAttribute("to")),
                        Long.parseLong(eElement.getAttribute("arcid"))
                );
                /*
                if (eElement.getAttribute("type") != null && !eElement.getAttribute("type").equals("")) {
                    System.out.println(eElement.getAttribute("type"));
                    graph.getWeightsMatrix().put(
                            Long.parseLong(eElement.getAttribute("arcid")),
                            10L,
                            getSpeed(eElement.getAttribute("type"))
                    );
                }
                */
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return graph;
    }

    // http://wiki.openstreetmap.org/wiki/OSM_tags_for_routing/Maxspeed
    private static float getSpeed(String type) {
        if (Objects.equals(type, "motorway")) {
            return 120f;
        }
        if (Objects.equals(type, "motorway_link") || Objects.equals(type, "trunk")) {
            return 100f;
        }
        if (Objects.equals(type, "trunk_link")) {
            return 80f;
        }
        if (Objects.equals(type, "living_street") || Objects.equals(type, "residential")) {
            return 20f;
        }
        return 50f;
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
        System.out.print("Adding weight to the graph...");
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
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        System.out.println("added");
        return graph;
    }

    private static Long getMaxArcId(GraphTable graph) {
        Long max = 0L;
        if (!graph.getWeightsMatrix().rowKeySet().isEmpty()) {
            List<Long> list = new ArrayList<>(graph.getWeightsMatrix().rowKeySet());
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
        System.out.print("Adding weight to the graph...");
        try {
            File fXmlFile = new File(file);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(fXmlFile);
            doc.getDocumentElement().normalize();

            NodeList nList = doc.getElementsByTagName("weight");
            int size = nList.getLength();
            Element eElement;
            for (int temp = 0; temp < size; temp++) {
                    eElement = (Element) nList.item(temp);
                    graph.getWeightsMatrix().put(
                            Long.parseLong(eElement.getAttribute("arcid")),
                            Long.parseLong(eElement.getAttribute("type")),
                            Float.parseFloat(eElement.getAttribute("value"))
                    );

            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("added");
        return graph;
    }

    public static GraphTable readWeights(GraphTable graph, String file) {
        System.out.print("Adding weight to the graph...");
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
            String line = br.readLine();

            while (line != null && line.length() > 0) {
                String[] array = line.split(" ");
                graph.getWeightsMatrix().put(
                        Long.parseLong(array[0]),
                        Long.parseLong(array[2]),
                        Float.parseFloat(array[1])
                );
                line = br.readLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        System.out.println("added");
        return graph;
    }

    public static void readTlLogics(GraphTable graph, String filename) {
        JSONParser parser = new JSONParser();
        try {

            Object obj = parser.parse(new FileReader(filename));

            JSONArray jsonRootObject = (JSONArray) obj;
            for (JSONObject item : (Iterable<JSONObject>) jsonRootObject) {
                // TODO change to use the mapping of the nodes
                // TODO String -> Long
                Long nodeFrom = Long.parseLong((String) item.get("node_from"));
                Long nodeTo = Long.parseLong((String) item.get("node_to"));
                TlLogic tl = new TlLogic();
                JSONArray phases = (JSONArray) item.get("phases");
                int i = 0;
                char[] aux  = new char[3];
                char current_light;
                for (JSONObject phase : (Iterable<JSONObject>) phases) {
                    current_light = ((String) phase.get("state")).charAt(0);
                    //System.out.println(current_light);
                    if (tl.getTime(current_light) == 0) {
                        //System.out.println(i + " " + phase);
                        aux[i] = current_light;
                        tl.setType(i, aux[i]);
                        i++;
                    }
                    tl.addTime(current_light, Integer.parseInt((String) phase.get("duration")));
                }
                tl.setType(aux);
                if (graph.getMapping().get(nodeFrom) != null && graph.getMapping().get(nodeTo) != null) {
                    graph.getTlMatrix().put(graph.getMapping().get(nodeFrom), graph.getMapping().get(nodeTo), tl);
                } else {
                    System.out.print(".");
                }
            }
            System.out.println(graph.getTlMatrix().rowKeySet().size());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void getConnectedComponents(GraphTable graph) {
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

        for (Long label : labels) {
            System.out.println("L: " + label);
        }

        for (Long node : ids.keySet()) {
            System.out.println("<" + node + " " + ids.get(node) + ">");
        }

        Map<Long, Integer> results = new HashMap<Long, Integer>(); // <nodeId, numElementInComponent>
        for (Long label : labels) {
            if (results.get(label) == null) {
                results.put(label, 1);
            } else {
                results.put(label, results.get(label) + 1);
            }
        }
        System.out.println("Connected Components");
        for (Long node : results.keySet()) {
            System.out.println("<" + node + " " + results.get(node) + ">");
        }
    }

    /**
     * Rename a label
     * @param oldLabel old label name
     * @param newLabel new label name
     * @param labels list of labels to be changed
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
        //Long[] labels = new Long[graph.getIntersections().size()];
        int i = 0;
        for (Long l : graph.getIntersections().keySet()) {
            ids.put(i, l);
            nodes.put(l, i);
            //labels[i] = l;
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
            if (!visited[i])
                fillOrder(i, visited, stack, graph, ids, nodes, true);

        // Create a reversed graph
        //Graph gr = getTranspose();

        // Mark all the vertices as not visited (For second DFS)
        for (i = 0; i < graph.getIntersections().size(); i++)
            visited[i] = false;

        // Now process all vertices in order defined by Stack
        while (!stack.empty()) {
            // Pop a vertex from stack
            int v = (Integer) stack.pop();

            // Print Strongly connected component of the popped vertex
            if (!visited[v]) {
                DFSUtil(v, visited, graph, ids, nodes, false);
                System.out.println();
            }
        }
    }

    // A recursive function to print DFS starting from v
    private static void DFSUtil(int v, boolean visited[], final GraphTable graph, final Map<Integer, Long> ids,
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


    private static void fillOrder(int v, boolean visited[], Stack stack, final GraphTable graph, final Map<Integer,
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
                if (nodes.get(n)!= null && !visited[nodes.get(n)])
                    fillOrder(nodes.get(n), visited, stack, graph, ids, nodes, order);
            }
        } else {
            graph.getAdjacencyMatrix().column(ids.get(v));
            for (Long n : graph.getAdjacencyMatrix().column(ids.get(v)).keySet()) {
                if (nodes.get(n)!= null && !visited[nodes.get(n)])
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
            for (String aBestCC : bestCC) {
                convexComponent.add(Long.parseLong(aBestCC));
            }

            // remove node
            Iterator iterator = graph.getIntersections().keySet().iterator();
            while (iterator.hasNext()) {
                if (!convexComponent.contains(iterator.next())) {
                    iterator.remove();
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
            try {
                if (writer != null) {
                    writer.close();
                }
            } catch (Exception ex) {ex.printStackTrace();}
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
            for (Long arcId : graph.getWeightsMatrix().rowKeySet()) {
                value = ((( upperBound - lowerBound) * random.nextFloat()) + lowerBound)
                        * graph.getWeightsMatrix().get(arcId, 4L);
                line = "<weight " +
                        "arcid='" + arcId +
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
            for (Long arcId : graph.getWeightsMatrix().rowKeySet()) {
                value = ((( upperBound - lowerBound) * random.nextFloat()) + lowerBound)
                        * graph.getWeightsMatrix().get(arcId, baseType);
                line = "<weight " +
                        "arcid='" + arcId +
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

    public static GraphTable computeRandomWeights(GraphTable graph, Long baseType, float lowerBound, float upperBound, Long type) {
        Random random = new Random(type);
        float value;
        for (Long arcId : graph.getWeightsMatrix().rowKeySet()) {
            value = ((( upperBound - lowerBound) * random.nextFloat()) + lowerBound)
                    * graph.getWeightsMatrix().get(arcId, baseType);
            graph.getWeightsMatrix().put(arcId, type, value);
        }
        return graph;
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

    public static GraphTable normalize(GraphTable graph) {
        System.out.print("Normalizing...");
        Float value;
        Float min = Float.MAX_VALUE;
        for (Long type : graph.getWeightsMatrix().columnKeySet()) {
            //System.out.println("analizing type: " + type);
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
            sd = (float) Math.sqrt(sumSquare.doubleValue() - (avg.doubleValue() * avg.doubleValue()));
            //System.out.println("normalizing: " + type);
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

        //System.out.println("moving distributions");
        min = -1 * min;//min;
        for (Long r : graph.getWeightsMatrix().rowKeySet()) {
            for (Long c : graph.getWeightsMatrix().row(r).keySet()) {
                graph.getWeightsMatrix().put(r, c, graph.getWeightsMatrix().get(r, c) + min);
            }
        }
        System.out.println("end");
        return graph;
    }

    public static GraphTable computeNewWeight(GraphTable graph, Long type, Long inverseType, Long w1, Long w2) {
        // w1 / w2
        Float value;
        for (Long arc : graph.getWeightsMatrix().rowKeySet()) {
            value = graph.getWeightsMatrix().get(arc, w1) / graph.getWeightsMatrix().get(arc, w2);
            graph.getWeightsMatrix().put(arc, type, value);
            graph.getWeightsMatrix().put(arc, inverseType, -1 * value);
        }
        return graph;
    }

    public static GraphTable applyMapping(GraphTable graph, String file) {
        System.out.print("Adding mapping to the graph...");
        BufferedReader br = null;
        Map<Long, Long> mapping = new HashMap<Long, Long>();
        try {
            br = new BufferedReader(new FileReader(file));
            String line = br.readLine();

            while (line != null && line.length() > 0) {
                String[] array = line.split(" ");
                mapping.put(Long.parseLong(array[0]), Long.parseLong(array[1]));
                line = br.readLine();
            }
            graph.setMapping(mapping);
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        System.out.println("added");
        return graph;
    }

    public static void printMapping(GraphTable graph) {
        Map<Long, Long> mapping = new HashMap<>();
        Long index = 1L;
        for (Long l : graph.getIntersections().keySet()) {
            mapping.put(l, index);
            index++;
        }
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("mapping-malaga.txt"));
            String line;
            for (Long k : mapping.keySet()) {
                line = k + " " + mapping.get(k) + "\n";
                out.write(line);
            }

        } catch (IOException e) {
            // error processing code
        } finally {
            if (out != null) {
                try {
                    out.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * Apply the mapping to the intersections and adjacencyMatrix
     * @param graph the graph
     * @return the graph
     */
    public static GraphTable fixVertexIndex(GraphTable graph) {
        Table<Long, Long, Long> newAdjacencyMatrix = HashBasedTable.create();
        Map<Long, Node> newIntersections = new HashMap<Long, Node>();

        // ya tengo mapping
        for (Long l : graph.getIntersections().keySet()) {
            newIntersections.put(graph.getMapping().get(l), graph.getIntersections().get(l));
        }
        for (Long r : graph.getAdjacencyMatrix().rowKeySet()) {
            for (Long c : graph.getAdjacencyMatrix().row(r).keySet()) {
//                System.err.println("(" + graph.getMapping().get(r) + "," + graph.getMapping().get(c) + ")->" + graph.getAdjacencyMatrix().get(r, c) + ",  ");
                newAdjacencyMatrix.put(
                        graph.getMapping().get(r),
                        graph.getMapping().get(c),
                        graph.getAdjacencyMatrix().get(r, c)
                );
            }
        }
        graph.setIntersections(newIntersections);
        graph.setAdjacencyMatrix(newAdjacencyMatrix);
        return graph;
    }

    public static List<Long> nodeToLong(GraphTable graph, List<Node> pathN) {
        List<Long> path = new ArrayList<>();
        for (Node n : pathN) {
            path.add(graph.getMapping().get(n.getId()));
        }
        return path;
    }

    public static GraphTable divideWeights(GraphTable graph, long w1, long w2, long res, float k) {
        for (Long arc : graph.getWeightsMatrix().rowKeySet()) {
            graph.getWeightsMatrix().put(arc, res,
                    k * graph.getWeightsMatrix().get(arc, w1) / graph.getWeightsMatrix().get(arc, w2));
        }
        return graph;
    }

    public static GraphTable addValuesGraph(GraphTable graph, String filename) {
        // load data
        BufferedReader br = null;
        List<Float[]> input = new ArrayList<>();
        try {
            br = new BufferedReader(new FileReader(filename));
            String line = br.readLine();

            while (line != null) {
                String[] a = line.split(" ");
                input.add(new Float[]{
                        Float.parseFloat(a[0]),
                        Float.parseFloat(a[1]),
                        Float.parseFloat(a[2]),
                        Float.parseFloat(a[3])});
                line = br.readLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        // compute missing information
        //asignateFixInformation(graph, input);

        for (Long node : graph.getIntersections().keySet()) {
            calculateInverseMeanNoise(graph.getIntersections().get(node), input);
        }

        // add to graph
        for (Long r : graph.getAdjacencyMatrix().rowKeySet()) {
            for (Long c : graph.getAdjacencyMatrix().row(r).keySet()) {
                graph.getWeightsMatrix().put(graph.getAdjacencyMatrix().get(r, c), 1L,
                        (graph.getIntersections().get(r).getNoiseMean() + graph.getIntersections().get(c).getNoiseMean())/2);
                graph.getWeightsMatrix().put(graph.getAdjacencyMatrix().get(r, c), 3L,
                        (graph.getIntersections().get(r).getNoiseSD() + graph.getIntersections().get(c).getNoiseSD())/2);
            }
        }

        return graph;
    }

    private static void calculateInverseMeanNoise(Node node, List<Float[]> input) {
        int dim = input.size();
        float[] distances = new float[dim];
        float sum = 0.0f;
        for (int i = 0; i < dim; i++) {
            distances[i] = invDistance(node.getLatitude(), node.getLongitude(), input.get(i)[0], input.get(i)[1]);
            if (distances[i] == -1) {
                node.setNoiseMean(input.get(i)[2]);
                node.setNoiseSD(input.get(i)[3]);
                return;
            }
            sum += distances[i];
        }

        float mean = 0.0f, sd = 0.0f;
        for (int i = 0; i < dim; i++) {
            mean += input.get(i)[2]*distances[i]/sum;
            sd += input.get(i)[3]*distances[i]/sum;
        }
        node.setNoiseMean(mean/dim);
        node.setNoiseSD(sd/dim);
    }

    private static float invDistance(BigDecimal latitude, BigDecimal longitude, Float lat, Float lon) {
        double dist = distance(latitude.doubleValue(), longitude.doubleValue(), lat, lon, "K");
        if (dist <= 0.001) { // 1 meter
            return -1;
        }
        return (float) dist;// -1 if is the same point
    }

    private static double distance(double lat1, double lon1, double lat2, double lon2, String unit) {
        double theta = lon1 - lon2;
        double dist = Math.sin(deg2rad(lat1)) * Math.sin(deg2rad(lat2)) + Math.cos(deg2rad(lat1)) * Math.cos(deg2rad(lat2)) * Math.cos(deg2rad(theta));
        dist = Math.acos(dist);
        dist = rad2deg(dist);
        dist = dist * 60 * 1.1515;
        if (unit.equals("K")) {
            dist = dist * 1.609344;
        } else if (unit.equals("N")) {
            dist = dist * 0.8684;
        }

        return (dist);
    }
    /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	/*::	This function converts decimal degrees to radians						 :*/
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    private static double deg2rad(double deg) {
        return (deg * Math.PI / 180.0);
    }

    /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	/*::	This function converts radians to decimal degrees						 :*/
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    private static double rad2deg(double rad) {
        return (rad * 180 / Math.PI);
    }

    public static GraphTable prepareGraph(String graphFilePath, String weightFilePath0, String mapping, String tag) {
        // Graph
        GraphTable graph = ProcessGraph.parserFile(graphFilePath);
        if (tag.equals("COL")) {
            readWeights(graph, weightFilePath0);
        } else {
            applyWeights(graph, weightFilePath0);
        }

        assert graph != null;
        graph.getWeightsMatrix().column(10L).clear();
        ProcessGraph.applyMapping(graph, mapping);
        return graph;
    }
}
