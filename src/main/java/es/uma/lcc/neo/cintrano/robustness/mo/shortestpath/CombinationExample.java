package es.uma.lcc.neo.cintrano.robustness.mo.shortestpath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Code from http://stackoverflow.com/questions/127704/algorithm-to-return-all-combinations-of-k-elements-from-n
 * 2/02/17.
 */
public class CombinationExample {
    private static List<String[]> out = new ArrayList<String[]>();
    public static void main(String[] args){
        String[] arr = {"A","B","C","D","E","F"};
        combinations2(arr, 3, 0, new String[3]);

        System.out.println("---");
        for (String[] e:
             out) {
            System.out.println(Arrays.toString(e));
        }
    }

    static void combinations2(String[] arr, int len, int startPosition, String[] result){
        if (len == 0){
            System.out.println(Arrays.toString(result));
            String[] aux = new String[result.length];
            for (int i = 0; i < result.length; i++) {
                aux[i] = result[i];
            }
            out.add(aux);
            return;
        }
        for (int i = startPosition; i <= arr.length-len; i++){
            result[result.length - len] = arr[i];
            combinations2(arr, len-1, i+1, result);
        }
    }
}
