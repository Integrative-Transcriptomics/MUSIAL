/**
 * 
 */
package utilities;

import java.util.Arrays;
import java.util.Map;

/**
 * @author Alexander Seitz
 *
 */
public class Statistics 
{
    double[] data;
    int size;
    
    public Statistics(Map<Integer, Double> coverage){
    	this.data = new double[coverage.size()];
    	int i=0;
    	for(int key: coverage.keySet()){
    		this.data[i] = coverage.get(key);
    		i++;
    	}
    	this.size = coverage.size();
    }

    public Statistics(double[] data) 
    {
        this.data = data;
        size = data.length;
    }   

    public double getMean()
    {
        double sum = 0.0;
        for(double a : data)
            sum += a;
        return sum/size;
    }

    public double getVariance()
    {
        double mean = getMean();
        double temp = 0;
        for(double a :data)
            temp += (a-mean)*(a-mean);
        return temp/size;
    }

    public double getStdDev()
    {
        return Math.sqrt(getVariance());
    }

    public double median() 
    {
       Arrays.sort(data);

       if (data.length % 2 == 0) 
       {
          return (data[(data.length / 2) - 1] + data[data.length / 2]) / 2.0;
       } 
       return data[data.length / 2];
    }
}
