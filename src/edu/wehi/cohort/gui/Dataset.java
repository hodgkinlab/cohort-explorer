package edu.wehi.cohort.gui; // Don't forget to change


import java.util.List;
import java.io.IOException;
import org.ejml.simple.SimpleMatrix;

/**
 * one input excel file corresponds to one instance of this class; this is an
 * internal data structure independent of data format in Excel
 *
 * @author akan
 */
public class Dataset {

    final String name;

    // in these matrices:
    // rows: time points, then generations
    // columns: conditions, then replicates
    final SimpleMatrix counts;
    final SimpleMatrix cohorts;
    final SimpleMatrix missing;

    // 1st dim: time point
    // 2nd dim: conditions
    // 3rd dim: replicates
    final public double[][][] totalNumbers;

    final public double[][][] cohortGenreps;
    final public double[][][] cohortSumreps;
    final public double[][][] sumExcGen0reps;

    final public double[] cohortSum0r;
    final public double[] cohortSum1r;

    final public double[][][] meanDiveach;
    final public double[][][] totCohorteach;
    final public double[][][] mdnExcGen0reps;
    final public double[][][] totexcTP0each;
    
    final public boolean[][][] missingStatreps;
    final public boolean[][][] missingEx0reps;
    
    // cohort stats
    // rows: time points
    // columns: conditions
    final public double[][] meanDivNum;
    final public double[][] mdnExcGen0;
    final public double[][] totCohort;
    final public double[][] totCohortnotnorm;
    final public double[][] totExcTP0;
    final public boolean[][] missingStat;
    final public boolean[][] missingEx0;
    

    final SimpleMatrix timeStamps; // list of time stamps
    final List<String> conditionNames; // list of conditions names

    final public int numTimePoints;
    final public int lastGeneration;
    final public int numConditions;
    final public int numReplicates;

    public Dataset(
            String name,
            double[][][] totalNumbers,
            SimpleMatrix counts,
            SimpleMatrix cohorts,
            SimpleMatrix missing,
            SimpleMatrix timeStamps,
            List<String> conditionNames)
            throws IOException {

        this.name = name;
        this.totalNumbers = totalNumbers;
        this.counts = counts;
        this.cohorts = cohorts;
        this.missing = missing;
        this.timeStamps = timeStamps;
        this.conditionNames = conditionNames;

        this.numTimePoints = timeStamps.getNumElements();
        double dNumGen = (double) cohorts.numRows() / (double) numTimePoints;
        if ((dNumGen != Math.floor(dNumGen)) || Double.isInfinite(dNumGen)) {
            throw new IOException("data matrix dimensions mismatch");
        }
        this.lastGeneration = (int) dNumGen - 1;

        this.numConditions = conditionNames.size();
        double dNumRep = (double) cohorts.numCols() / (double) numConditions;
        if ((dNumRep != Math.floor(dNumRep)) || Double.isInfinite(dNumRep)) {
            throw new IOException("data matrix dimensions mismatch");
        }
        this.numReplicates = (int) dNumRep;

        if (totalNumbers.length != this.numTimePoints) {
            throw new IOException("tot. counts matrix num. t.p. mismatch");
        }
        if (totalNumbers[0].length != this.numConditions) {
            throw new IOException("tot. counts matrix num. cond. mismatch");
        }
        if (totalNumbers[0][0].length != this.numReplicates) {
            throw new IOException("tot. counts matrix num. rep. mismatch");
        }

        this.meanDivNum = new double[numTimePoints][numConditions];
        this.mdnExcGen0 = new double[numTimePoints][numConditions];
        this.totCohort = new double[numTimePoints][numConditions];
        this.totCohortnotnorm = new double[numTimePoints][numConditions];
        this.totExcTP0 = new double[numTimePoints][numConditions];
        this.missingStat = new boolean[numTimePoints][numConditions];
        this.missingEx0 = new boolean[numTimePoints][numConditions];
        this.missingStatreps = new boolean[numTimePoints][numConditions][numReplicates];
        this.missingEx0reps = new boolean[numTimePoints][numConditions][numReplicates];

        this.meanDiveach = new double[numTimePoints][numConditions][numReplicates];
        this.totCohorteach = new double[numTimePoints][numConditions][numReplicates];
        this.mdnExcGen0reps = new double[numTimePoints][numConditions][numReplicates];
        this.totexcTP0each = new double[numTimePoints][numConditions][numReplicates];

        this.cohortSum0r = new double[numReplicates];
        this.cohortSum1r = new double[numReplicates];

        this.cohortGenreps = new double[numTimePoints][numConditions][numReplicates];
        this.cohortSumreps = new double[numTimePoints][numConditions][numReplicates];
        this.sumExcGen0reps = new double[numTimePoints][numConditions][numReplicates];
        
        //This section is with repetition
        for (int icnd = 0; icnd < numConditions; icnd++) {

            for (int itpt = 0; itpt < numTimePoints; itpt++) {
                for (int irep = 0; irep < numReplicates; irep++) {
                    for (int igen = 0; igen <= lastGeneration; igen++) {
                        
                        missingStatreps[itpt][icnd][irep] = false;

                        int row = igen + itpt * (lastGeneration + 1);
                        int col = irep + icnd * numReplicates;

                        if (missing.get(row, col) == 0) {
                            cohortGenreps[itpt][icnd][irep] += cohorts.get(row, col) * igen; 
                            //sum through generations, times cohort weighted by gen
                            cohortSumreps[itpt][icnd][irep] += cohorts.get(row, col); 
                            //Sum of cohort number
                        }
                        else{
                            missingStatreps[itpt][icnd][irep] = true;
                            //flags missing data
                        }

                        if (igen > 0) {
                            sumExcGen0reps[itpt][icnd][irep] += cohorts.get(row, col); 
                            //sum of cohort through gen, excluding 0 gen
                        }
                    }

                    //Finished generation loops

                    if (itpt == 0 && !missingStatreps[itpt][icnd][irep]) {
                        cohortSum0r[irep] = cohortSumreps[itpt][icnd][irep];
                        //Sum of only first time point
                    }

                    if (itpt == 1 && !missingStatreps[itpt][icnd][irep]) {
                        cohortSum1r[irep] = cohortSumreps[itpt][icnd][irep];
                        //sum of only second time point onwards
                    } 

                    if (cohortSumreps[itpt][icnd][irep] > 0 && !missingStatreps[itpt][icnd][irep]) {
                        meanDiveach[itpt][icnd][irep] = cohortGenreps[itpt][icnd][irep] / cohortSumreps[itpt][icnd][irep];
                        //Cohort average, the average #times it has divided
                    }

                    if (sumExcGen0reps[itpt][icnd][irep] > 0 && !missingStatreps[itpt][icnd][irep]) {
                        mdnExcGen0reps[itpt][icnd][irep] = cohortGenreps[itpt][icnd][irep] / sumExcGen0reps[itpt][icnd][irep];
                        //Cohort average, the average #times it has divided excluding generation 0
                    }
                    else if (sumExcGen0reps[itpt][icnd][irep] == 0) {
                        missingEx0reps[itpt][icnd][irep] = true;
                    }


                    if (cohortSum0r[irep] > 0 && !missingStatreps[itpt][icnd][irep] && !missingStatreps[0][icnd][irep]) {
                        totCohorteach[itpt][icnd][irep] = cohortSumreps[itpt][icnd][irep] / cohortSum0r[irep];
                        //Normalising with first time point sum
                    }

                    if (cohortSum1r[irep] > 0 && !missingStatreps[itpt][icnd][irep] && !missingStatreps[1][icnd][irep]) {
                        totexcTP0each[itpt][icnd][irep] = cohortSumreps[itpt][icnd][irep] / cohortSum1r[irep];
                        //Normalising with second time point sum
                    }
                }

            }

        }
        
        //This Section is the original, without repitition

        for (int icnd = 0; icnd < numConditions; icnd++) {
            double cohortSum0 = 0;
            double cohortSum1 = 0;

            for (int itpt = 0; itpt < numTimePoints; itpt++) {
                double cohortSum = 0;
                double sumExcGen0 = 0;
                double cohortGen = 0;
                missingStat[itpt][icnd] = false;
                missingEx0[itpt][icnd] = false;

                for (int igen = 0; igen <= lastGeneration; igen++) {
                    double sumRep = 0;
                    int numRep = 0;

                    for (int irep = 0; irep < numReplicates; irep++) {

                        int row = igen + itpt * (lastGeneration + 1);
                        int col = irep + icnd * numReplicates;
                        if (missing.get(row, col) == 0) {
                            sumRep += cohorts.get(row, col);
                            numRep++;
                        }
                    }

                    if (numRep == 0) {//If there are no reps, time and cond that is no data
                        missingStat[itpt][icnd] = true;
                        missingEx0[itpt][icnd] = false;
                    } else {
                        double meanRep = sumRep / numRep; //This is just an average of replicate cohorts 
                        cohortGen += meanRep * igen; //This is looping through generations not replicates
                        cohortSum += meanRep;
                        if (igen > 0) {
                            sumExcGen0 += meanRep;
                        }
                    }
                }

                //Looping through time starts here
                //The missing stat thing here means !False=True (Only for zeroth gen)
                if ((itpt == 0) && !missingStat[itpt][icnd]) {
                    cohortSum0 = cohortSum;
                }

                //If no missing stat, cohort sum of time =1
                if ((itpt == 1) && !missingStat[itpt][icnd]) {
                    cohortSum1 = cohortSum;
                }

                //If no missing stat, cohort sum calculate MeanDivNum
                if (!missingStat[itpt][icnd]) {
                    if (cohortSum > 0) {
                        meanDivNum[itpt][icnd] = cohortGen / cohortSum;
                    } else {
                        missingStat[itpt][icnd] = true;
                    }
                }

                if (!missingEx0[itpt][icnd]) {
                    if (sumExcGen0 > 0) {
                        mdnExcGen0[itpt][icnd] = cohortGen / sumExcGen0;
                    } else {
                        missingEx0[itpt][icnd] = true;
                    }
                }
                if ((cohortSum0 > 0) && !missingStat[itpt][icnd]) {
                    totCohort[itpt][icnd] = cohortSum / cohortSum0;
                    totCohortnotnorm[itpt][icnd] = cohortSum;
                }
                // WHAT I THINK IT SUPPOSED TO BE DOING (wrong flag??)
                if ((cohortSum1 > 0) && !missingStat[itpt][icnd]) {
                    totExcTP0[itpt][icnd] = cohortSum / cohortSum1;
                }
                // OLD ONE
//                if ((cohortSum1 > 0) && !missingEx0[itpt][icnd]) {
//                    totExcTP0[itpt][icnd] = cohortSum / cohortSum1;
//                }
            }
        }
    }

    public String getName() {
        return name;
    }

    public SimpleMatrix getTimeStamps() {
        return timeStamps;
    }

    public SimpleMatrix getCohorts() {
        return cohorts;
    }

    public SimpleMatrix getMissing() {
        return missing;
    }

    public List<String> getConditionNames() {
        return conditionNames;
    }

    public SimpleMatrix getCounts() { return counts; }

    public double getCohort(int itpt, int icnd, int irep, int igen) {
        int row = igen + itpt * (lastGeneration + 1);
        int col = irep + icnd * numReplicates;
        return cohorts.get(row, col);
    }

    public double getMissing(int itpt, int icnd, int irep, int igen) {
        int row = igen + itpt * (lastGeneration + 1);
        int col = irep + icnd * numReplicates;
        return missing.get(row, col);
    }

    public boolean[][] getmissingStat() {
        return missingStat;
    }

    public boolean[][][] getMissingStatreps() {
        return missingStatreps;
    }
}
