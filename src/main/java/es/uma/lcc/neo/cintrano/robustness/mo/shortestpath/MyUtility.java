package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;

import java.io.*;
import java.util.List;
import java.util.Set;

public class MyUtility {

    static Long[] readPoints(int id, String city) {
        Long[] points = null;
        BufferedReader br = null;
        try {
            String input_file = null;
            switch (city) {
                case "Malaga": input_file = "input-MAL.data"; break;
                case "Colorado": input_file = "input-COL.data"; break;
                case "NY": input_file = "input-NY.data"; break;
            }
            assert input_file != null;
            br = new BufferedReader(new FileReader(input_file));
            String line = br.readLine();

            int count = 0;
            while (line != null && count <= id) {
                if (count == id) {
                    String[] array = line.split(" ");
                    points = new Long[]{Long.parseLong(array[0]), Long.parseLong(array[1])};
                }
                count++;
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
        return points;
    }

    static void printSummary(long timeInit, long timeStart, Set<NodePathSolution> solutions, long timeEnd) {
        System.out.println("NUMBER OF SOLUTIONS: " + solutions.size());
        System.out.println("... method end");

        System.out.println("Initialization time: \t" + (timeStart - timeInit));
        System.out.println("Execution time: \t" + (timeEnd - timeStart));
        System.out.println("... method end");

        System.out.println("................");
        System.out.println("\nPrint solutions FUN:");
    }


    static void printSolutions(GraphTable graph, NodePathSolution s, long t1, long t2, String algorithm, String i) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("resultsFUN.ssv", true));
            String line = graph.getMapping().get(s.getVariables()[0]) + " ";
            line += graph.getMapping().get(s.getVariables()[s.getVariables().length-1]) + " ";
            line += t1 + " ";
            line += t2 + " ";
            line += algorithm + " ";
            line += i + " ";
            line += s.getObjectives()[0] + " ";
            line += s.getObjectives()[1] + " ";
            line += s.getObjectives()[2] + " ";
            line += s.getObjectives()[3] + " ";
            line += "\n";
            out.write(line);

        } catch (IOException e) {
            // error processing code
        } finally {
            if (out != null) try {
                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    static void printSolutions(GraphTable graph, NodePathSolution s, long t1, long t2, String algorithm, String i, int seed) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("resultsFUN.ssv", true));
            String line = graph.getMapping().get(s.getVariables()[0]) + " ";
            line += graph.getMapping().get(s.getVariables()[s.getVariables().length-1]) + " ";
            line += seed + " ";
            line += t1 + " ";
            line += t2 + " ";
            line += algorithm + " ";
            line += i + " ";
            line += s.getObjectives()[0] + " ";
            line += s.getObjectives()[1] + " ";
            line += s.getObjectives()[2] + " ";
            line += s.getObjectives()[3] + " ";
            line += "\n";
            out.write(line);

        } catch (IOException e) {
            // error processing code
        } finally {
            if (out != null) try {
                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    static void printSolutions(GraphTable graph, Set<NodePathSolution> paths, long t1, long t2, String algorithm, String i, int seed) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("resultsFUN.ssv", true));
            int index = 0;
            for (NodePathSolution s : paths) {
                String line = graph.getMapping().get(s.getVariables()[0]) + " ";
                line += graph.getMapping().get(s.getVariables()[s.getVariables().length - 1]) + " ";
                line += seed + " ";
                line += t1 + " ";
                line += t2 + " ";
                line += algorithm + " ";
                line += i + " ";
                line += s.getObjectives()[0] + " ";
                line += s.getObjectives()[1] + " ";
//                line += s.getObjectives()[2] + " ";
//                line += s.getObjectives()[3] + " ";
                line += index + " ";
                line += "\n";
                out.write(line);
                index ++;
            }

        } catch (IOException e) {
            // error processing code
        } finally {
            if (out != null) try {
                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    static void printSolutions(Set<NodePathSolution> solutions, long t1, long t2, String algorithm, String tag) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("resultsFUN.ssv", true));
            int i = 0;
            String line;
            for (NodePathSolution s : solutions) {
                line = s.getVariables()[0] + " ";
                line += s.getVariables()[s.getVariables().length-1] + " ";
                line += t1 + " ";
                line += t2 + " ";
                line += algorithm + " ";
                line += i + "-" + tag + " ";
                line += s.getObjectives()[0] + " ";
                line += s.getObjectives()[1] + " ";
                line += s.getObjectives()[2] + " ";
                line += s.getObjectives()[3] + " ";
                line += "\n";
                out.write(line);
                i++;
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




    static void printSolutions(Set<NodePathSolution> solutions, long t1, long t2, String algorithm) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("resultsFUN.ssv", true));
            int i = 0;
            String line;
            for (NodePathSolution s : solutions) {
                line = s.getVariables()[0] + " ";
                line += s.getVariables()[s.getVariables().length-1] + " ";
                line += t1 + " ";
                line += t2 + " ";
                line += algorithm + " ";
                line += i + " ";
                line += s.getObjectives()[0] + " ";
                line += s.getObjectives()[1] + " ";
                line += s.getObjectives()[2] + " ";
                line += s.getObjectives()[3] + " ";
                line += "\n";
                out.write(line);
                i++;
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

    public static float variance(List<Double> samples, double mean) {
        if (samples.size() == 1) return 0;
        float var = 0;
        for (Double sample : samples) {
            var += Math.pow(sample - mean, 2f);
        }
        var = var / (samples.size() -1);
        return var;
    }
}
