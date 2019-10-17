package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.Node;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.TlLogic;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
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
            line += s.getTl() + " ";
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
                line += (t1 + t2) + " ";
                line += s.getTl() + " ";
                line += algorithm + " ";
                line += i + " ";
                line += s.getObjectives()[0] + " ";
                line += s.getObjectives()[1] + " ";
                line += s.getObjectives()[2] + " ";
                line += s.getObjectives()[3] + " ";
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

        try {
            out = new BufferedWriter(new FileWriter("resultsVAR.ssv", true));
            int index = 0;
            for (NodePathSolution s : paths) {
                String line = graph.getMapping().get(s.getVariables()[0]) + " ";
                line += graph.getMapping().get(s.getVariables()[s.getVariables().length - 1]) + " ";
                line += seed + " ";
                line += t1 + " ";
                line += t2 + " ";
                line += algorithm + " ";
                line += i + " ";
                for (Long l : s.getVariables()) {
                    line += l + " ";
                }
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


    public static float[] fitnessTL(GraphTable graph, long[] s, int numObjectives, int MAX_SAMPLES, Random rand) {
        float[] fit = new float[numObjectives];
        for (int objective = 0; objective < numObjectives; objective += 2) {
            int tlCount = 0;
            long n1, n2;
            float amount = 0;
            List<Double> samples = new ArrayList<>();
            for (int sample = 0; sample < MAX_SAMPLES; sample++) {
                double cost = 0;
                for (int i = 0; i < s.length - 1; i++) {
                    n1 = graph.getMapping().get(s[i]);
                    n2 = graph.getMapping().get(s[i + 1]);
                    long arc = graph.getAdjacencyMatrix().get(n1, n2);
                    float mu = graph.getWeightsMatrix().get(arc, (long) objective);
                    float sigma = graph.getWeightsMatrix().get(arc, (long) objective + 1);
                    float tentativeTime = Math.abs((float) rand.nextGaussian() * sigma + mu);
                    cost += tentativeTime; // TODO Add drive profile
                    // TrafficLight phase
                    TlLogic tl = graph.getTlMatrix().get(n1, n2);
                    if (tl != null) {
                        tlCount++;
                        int time = tl.calculateTimeStop(Math.round(cost + 0.5d));
                        cost += time;
                    }
                }
                samples.add(cost);
                amount += cost;
            }
            float mean = amount / (float) MAX_SAMPLES;
            float variance = MyUtility.variance(samples, mean);
            System.out.println("F " + mean + " " + variance + " " + tlCount);
            fit[objective] = mean;
            fit[objective+1] = variance;
        }
        return fit;
    }
}
