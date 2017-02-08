package es.uma.lcc.neo.robustness.mo.shortestpath.algorithm;

/**
 * Created by cintrano on 7/02/17.
 */
public class MyPair implements Comparable{
    private Float k;
    private Long v;


    public MyPair(Float k, Long v) {

        this.k = k;
        this.v = v;
    }

    public Float getK() {
        return k;
    }

    public void setK(Float k) {
        this.k = k;
    }

    public Long getV() {
        return v;
    }

    public void setV(Long v) {
        this.v = v;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MyPair myPair = (MyPair) o;

        if (Float.compare(myPair.k, k) != 0) return false;
        return Float.compare(myPair.v, v) == 0;
    }

    @Override
    public int hashCode() {
        int result = (k != +0.0f ? Float.floatToIntBits(k) : 0);
        result = 31 * result + (v != +0.0f ? Float.floatToIntBits(v) : 0);
        return result;
    }

    public int compareTo(Object o) {

        return k.compareTo(((MyPair) o).getK());
    }

}
