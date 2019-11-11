package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.metaheuristics;

import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.solution.Solution;

import java.util.HashMap;
import java.util.Map;

public class MyDoubleSolution implements DoubleSolution {
    private double[] objectives;
    private Double[] variables;
    private Map<String, Integer> attributes;

    MyDoubleSolution(double[] objectives, Double[] variables) {
        this.objectives = objectives;
        this.variables = variables;
        attributes = new HashMap<>() ;
    }

    @Override
    public Double getLowerBound(int k) {
        int i = k  % variables.length;
        double d = variables[0];
        for (int j = 0; j < variables.length; j++) {
            if (d > variables[i]) d = variables[i];
        }
        return d;
    }

    @Override
    public Double getUpperBound(int k) {
        int i = k  % variables.length;
        double d = variables[0];
        for (int j = 0; j < variables.length; j++) {
            if (d < variables[i]) d = variables[i];
        }
        return d;
    }

    @Override
    public void setObjective(int i, double v) {
        this.objectives[i] = v;
    }

    @Override
    public double getObjective(int i) {
        return this.objectives[i];
    }

    @Override
    public double[] getObjectives() {
        return this.objectives;
    }

    public Double[] getVariables() {
        return this.variables;
    }

    @Override
    public Double getVariableValue(int i) {
        return this.variables[i % variables.length];
    }

    @Override
    public void setVariableValue(int i, Double aDouble) {
        this.variables[i % variables.length] = aDouble;
    }

    public void setVariables(Double[] aDouble) {
        this.variables = aDouble;
    }


    @Override
    public String getVariableValueString(int i) {
        return this.variables[i % variables.length].toString();
    }

    @Override
    public int getNumberOfVariables() {
        return variables.length;
    }

    @Override
    public int getNumberOfObjectives() {
        return objectives.length;
    }

    @Override
    public Solution<Double> copy() {
        double[] objectivesN = objectives.clone();
        Double[] variablesN = variables.clone();
        MyDoubleSolution aux = new MyDoubleSolution(objectivesN, variablesN);
        for (String k : attributes.keySet()) {
            aux.setAttribute(k, attributes.get(k));
        }
        return aux;
    }


    @Override
    public void setAttribute(Object id, Object value) {
        attributes.put((String) id, (Integer) value);
    }

    @Override
    public Object getAttribute(Object id) {
        return attributes.get(id) ;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder("[");
        for (Double var : variables) {
            s.append(var).append(",");
        }
        s.append("] ");

        for (Double var : objectives) {
            s.append(var).append(",");
        }

        return s.toString();
    }
}

