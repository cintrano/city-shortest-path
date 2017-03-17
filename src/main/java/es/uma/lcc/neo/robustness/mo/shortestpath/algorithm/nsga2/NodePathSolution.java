package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm.nsga2;

import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.solutionattribute.SolutionAttribute;

import java.util.HashMap;
import java.util.Map;


/**
 * Created by christian on 5/12/16.
 * Path solution in the graph
 */
public class NodePathSolution implements Solution<Long>, SolutionAttribute<NodePathSolution, Integer> {
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
        /*
        System.out.print("Node " + variables.length +" " + i + " --> ");
        for (Long l : variables) {
            System.out.print(" " + l);
        }
        System.out.println();
        */
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

    private Map<Object, Object> attributes = new HashMap<Object, Object>();

    public void setAttribute(Object id, Object value) {
        attributes.put(id, value);
    }

    public Object getAttribute(Object id) {
        return attributes.get(id);
    }
/*
    public void setAttribute(Object o, Object o1) {

    }

    public Object getAttribute(Object o) {
        return null;
    }

    public Integer getAttribute(NodePathSolution solution) {
        return solution.getAttribute(getAttributeID());
    }

    public void setAttribute(NodePathSolution solution, Object value) {
        solution.setAttribute(getAttributeID(), value);
    }

    public Object getAttributeID() {
        return this.getClass() ;
    }
*/

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
        String s = "[";
        for (Long var : variables) {
            s += var + ",";
        }
        s += "] ";

        for (Double var : objectives) {
            s += var + ",";
        }

        return s;
    }

    public double[] getObjectives() {
        return objectives;
    }
}
