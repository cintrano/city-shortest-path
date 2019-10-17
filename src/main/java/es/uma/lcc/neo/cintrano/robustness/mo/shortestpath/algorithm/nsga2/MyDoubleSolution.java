package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2;

import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.solution.Solution;

import java.util.HashMap;
import java.util.Map;

public class MyDoubleSolution implements DoubleSolution {
    double[] objectives;
    Double[] variables;
    protected Map<Object, Object> attributes ;

    public MyDoubleSolution(double[] objectives, Double[] variables) {
        this.objectives = objectives;
        this.variables = variables;
        attributes = new HashMap<Object, Object>() ;
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

    public void setVariables(Double[] variables) {
        this.variables = variables;
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
        return new MyDoubleSolution(objectivesN, variablesN);
    }


    @Override
    public void setAttribute(Object id, Object value) {
        attributes.put(id, value) ;
    }

    @Override
    public Object getAttribute(Object id) {
        return attributes.get(id) ;
    }

    public Map<Object, Object> getAttributes() {
        return attributes;
    }
}

