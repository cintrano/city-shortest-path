package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath;

import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.MyDoubleSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.GraphTable;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.NodePathSolution;
import es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.model.graph.guava.TlLogic;
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.solution.DoubleSolution;

import java.io.*;
import java.nio.charset.StandardCharsets;
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

    static void printStaticFinalLog(Algorithm<List<es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.NodePathSolution>> algorithm, long time, String label) {
        Writer writer = null;
        String filename = "result_F.csv";
        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filename), StandardCharsets.UTF_8));
            int i, tls;
            long iteration = 111L;
            i = 0;
            for (es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.NodePathSolution sol : algorithm.getResult()) {
                System.out.print("SOL:: ");
                System.out.println(sol);
                System.out.println(sol.getAttribute("tl"));
                tls = (Integer) sol.getAttribute("tl");
                writer.write(iteration + " " + i + " "+ time + " " + tls + " ");
                for (double e : sol.getObjectives()) {
                    writer.write(e + " ");
                }
                i++;
                writer.write("\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {
                assert writer != null;
                writer.close();} catch (Exception ex) {ex.printStackTrace();}
        }
        filename = "result_V.csv";
        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filename), StandardCharsets.UTF_8));
            int i = 0;
            long iteration = 111L;
            for (es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics.NodePathSolution sol : algorithm.getResult()) {
                writer.write(iteration + " " + i + " ");
                for (Long e : sol.getVariables()) {
                    writer.write(e + " ");
                }
                i++;
                writer.write("\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {writer.close();} catch (Exception ex) {ex.printStackTrace();}
        }
    }

    static void printStaticFinalLog(Algorithm<List<MyDoubleSolution>> algorithm, long time) {
        Writer writer = null;
        String filename = "result_F.csv";
        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filename), StandardCharsets.UTF_8));
            int i, tls;
            long iteration = 111L;
            i = 0;
            for (MyDoubleSolution sol : algorithm.getResult()) {
                tls = (Integer) sol.getAttribute("tl");
                writer.write(iteration + " " + i + " "+ time + " " + tls + " ");
                for (double e : sol.getObjectives()) {
                    writer.write(e + " ");
                }
                i++;
                writer.write("\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {
                assert writer != null;
                writer.close();} catch (Exception ex) {ex.printStackTrace();}
        }
        filename = "result_V.csv";
        try {
            writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(filename), StandardCharsets.UTF_8));
            int i = 0;
            long iteration = 111L;
            for (DoubleSolution sol : algorithm.getResult()) {
                writer.write(iteration + " " + i + " ");
                for (int j = 0; j < sol.getNumberOfVariables(); j++) {
                    writer.write(sol.getVariableValue(j) + " ");
                }
                i++;
                writer.write("\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            // report
        } finally {
            try {writer.close();} catch (Exception ex) {ex.printStackTrace();}
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


    public static float[] fitnessTL(GraphTable graph, long[] s, int numObjectives, int MAX_SAMPLES, boolean robustFlag, Random rand) {
        float[] fit = new float[numObjectives + 1]; // To return the num of TL as the last one
        int tlCount = 0;
        long n1, n2;
        float[] mu, sigma;
        float[] sum = new float[2];
        float[] sumSq = new float[2];

        for (int sample = 0; sample < MAX_SAMPLES; sample++) {
            float[] cost = new float[numObjectives/2];
            for (int i = 0; i < s.length - 1; i++) {
                n1 = graph.getMapping().get(s[i]);
                n2 = graph.getMapping().get(s[i + 1]);
                long arc = graph.getAdjacencyMatrix().get(n1, n2);
                mu = new float[numObjectives/2];
                sigma = new float[numObjectives/2];

                for (int objective = 0; objective < numObjectives/2; objective ++) {
                    mu[objective] = graph.getWeightsMatrix().get(arc, (long) objective);
                    sigma[objective] = graph.getWeightsMatrix().get(arc, (long) objective + 1);
                }
                float tentativeTime = Math.abs((float) rand.nextGaussian() * sigma[0] + mu[0]); // time
                float tentativePollution = Math.abs((float) rand.nextGaussian() * sigma[1] + mu[1]);//mu[0] * (tentativeTime/mu[0]); // CO2
                cost[0] = tentativeTime; // TODO Add drive profile
                cost[1] = tentativePollution;
                // TrafficLight phase
                TlLogic tl = graph.getTlMatrix().get(n1, n2);
                if (tl != null) {
                    tlCount++;
                    int time = tl.calculateTimeStop(Math.round(cost[0] + 0.5d));
                    cost[0] += time;
                    float pollution = 2166.094f * (float) time; // c0 * time
                    cost[1] += pollution;
                }
            }
            // MEANS
            fit[0] = (fit[0] * sample + cost[0]) / (sample + 1);
            fit[2] = (fit[2] * sample + cost[1]) / (sample + 1);
            if (robustFlag) {
                // Variances
                sum[0] += cost[0];
                sumSq[0] += cost[0] * cost[0];
                sum[1] += cost[1];
                sumSq[1] += cost[1] * cost[1];
            }
        }
        if (robustFlag) {
            fit[1] = (sumSq[0] - (sum[0] * sum[0]) / (float) MAX_SAMPLES) / (float) MAX_SAMPLES;
            fit[3] = (sumSq[1] - (sum[1] * sum[1]) / (float) MAX_SAMPLES) / (float) MAX_SAMPLES;
        } else {
            fit[1] = 0;
            fit[3] = 0;
        }

        fit[fit.length - 1] = tlCount;
        return fit;
    }

    public static float[] fitnessTLMeta(GraphTable graph, long[] s, int numObjectives, int MAX_SAMPLES, boolean robustFlag, Random rand) {
        float[] fit = new float[numObjectives + 1]; // To return the num of TL as the last one
        int tlCount = 0;
        long n1, n2;
        float[] mu, sigma;
        float[] sum = new float[2];
        float[] sumSq = new float[2];

        for (int sample = 0; sample < MAX_SAMPLES; sample++) {
            float[] cost = new float[numObjectives/2];
            for (int i = 0; i < s.length - 1; i++) {
                n1 = s[i];
                n2 = s[i + 1];
                long arc = graph.getAdjacencyMatrix().get(n1, n2);
                mu = new float[numObjectives/2];
                sigma = new float[numObjectives/2];

                for (int objective = 0; objective < numObjectives/2; objective ++) {
                    mu[objective] = graph.getWeightsMatrix().get(arc, (long) objective);
                    sigma[objective] = graph.getWeightsMatrix().get(arc, (long) objective + 1);
                }
                float tentativeTime = Math.abs((float) rand.nextGaussian() * sigma[0] + mu[0]); // time
                float tentativePollution = Math.abs((float) rand.nextGaussian() * sigma[1] + mu[1]);//mu[0] * (tentativeTime/mu[0]); // CO2
                cost[0] = tentativeTime; // TODO Add drive profile
                cost[1] = tentativePollution;
                // TrafficLight phase
                TlLogic tl = graph.getTlMatrix().get(n1, n2);
                if (tl != null) {
                    tlCount++;
                    int time = tl.calculateTimeStop(Math.round(cost[0] + 0.5d));
                    cost[0] += time;
                    float pollution = 2166.094f * (float) time; // c0 * time
                    cost[1] += pollution;
                }
            }
            // MEANS
            fit[0] = (fit[0] * sample + cost[0]) / (sample + 1);
            fit[2] = (fit[2] * sample + cost[1]) / (sample + 1);
            if (robustFlag) {
                // Variances
                sum[0] += cost[0];
                sumSq[0] += cost[0] * cost[0];
                sum[1] += cost[1];
                sumSq[1] += cost[1] * cost[1];
            }
        }
        if (robustFlag) {
            fit[1] = (sumSq[0] - (sum[0] * sum[0]) / (float) MAX_SAMPLES) / (float) MAX_SAMPLES;
            fit[3] = (sumSq[1] - (sum[1] * sum[1]) / (float) MAX_SAMPLES) / (float) MAX_SAMPLES;
        } else {
            fit[1] = 0;
            fit[3] = 0;
        }

        fit[fit.length - 1] = tlCount;
        return fit;
    }
}
