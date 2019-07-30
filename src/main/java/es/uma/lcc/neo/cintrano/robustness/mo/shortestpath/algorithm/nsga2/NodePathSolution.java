package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath.algorithm.nsga2;

import org.uma.jmetal.solution.Solution;

import java.util.HashMap;
import java.util.Map;


/**
 * Created by christian on 5/12/16.
 * Path solution in the graph
 */
public class NodePathSolution implements Solution<Long> {
//public class NodePathSolution implements RankingSolution, CrowdingDistanceSolution {
//public class NodePathSolution implements SolutionAttribute<NodePathSolution, V> {
    private double[] objectives;
    private Long[] variables;

    public NodePathSolution(double[] objectives, Long[] variables) {
        this.objectives = objectives;
        this.variables = variables;
    }

    public void setObjective(int i, double v) {
        objectives[i] = v;
    }

    public double getObjective(int i) {
        return objectives[i];
    }

    public void setVariables(Long[] variables) {
        this.variables = variables;
    }

    public Long getVariableValue(int i) {
        return variables[i];
    }

    public void setVariableValue(int i, Long t) {
        variables[i] = t;
    }

    public String getVariableValueString(int i) {
        return variables[i].toString();
    }

    public int getNumberOfVariables() {
        return variables.length;
    }

    public int getNumberOfObjectives() {
        return objectives.length;
    }

    public NodePathSolution copy() {
        return this;
    }

    @Override
    public void setAttribute(Object o, Object o1) {

    }

    @Override
    public Object getAttribute(Object o) {
        return null;
    }

    private Map<Object, Object> attributes = new HashMap<>();

    public Integer getAttribute(NodePathSolution solution) {
        return (Integer) solution.getAttribute(getAttributeID());
    }

    public void setAttribute(NodePathSolution solution, Integer value) {
        solution.setAttribute(getAttributeID(), value);
    }

    public Object getAttributeID() {
        return this.getClass() ;
    }

    public Long[] getVariables() {
        return variables;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder("[");
        for (Long var : variables) {
            s.append(var).append(",");
        }
        s.append("] ");

        for (Double var : objectives) {
            s.append(var).append(",");
        }

        return s.toString();
    }

    public double[] getObjectives() {
        return objectives;
    }
}
