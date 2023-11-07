package edu.wehi.cohort.format;


import java.io.*;
import java.util.Scanner;
import java.nio.charset.StandardCharsets;
import java.util.Random;

import javax.swing.*;
import javax.swing.border.Border;

import de.erichseifert.gral.data.DataTable;

import edu.wehi.cohort.gui.Dataset;


/**
 * @author HoChan 6/7/17.
 */

public class Bootstrap {

    final private int numIter;
    final private double[][][] AStar;
    final private double[][][] muStar;
    final private double[][][] sigmaStar;

    final private Random randGen = new Random();

    public Bootstrap(Dataset DATA,
                     double[][][] origX, double[][][][] origY,
                     boolean[][] missingStat,
                     int[][][] numRepDataCounter,
                     JTextArea text,
                     JProgressBar progressBar,
                     File inputFile){

        /* ***************************************************************************
         * Reads 'settings.txt' file for bootstrap iterations
         *   If not exists create one -> set to default iteration
         *   Exists but unexpected format -> set to default iteration with error msgs
         *****************************************************************************/
        boolean fileExists = false;
        boolean foundNumber = false;
        double foundBootIter = 0;
        int defaultBootIter = 5;

        try {
            File textfile = new File(inputFile.getParent() + "/settings.txt");

            // create new 'settings.txt' file in the working directory
            if (!textfile.exists()) {
                text.append("     ***'settings.txt' file not found! creating one... setting default bootstrap iteration number\n");

                Writer out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(textfile), StandardCharsets.UTF_8));
                out.write("=================================================================================\r\n" +
                        "Set number of bootstrap iterations in this file!\r\n" +
                        "[Figure 2] - Discrete Gaussian model: \r\n" +
                        "   Total iterations = (# time points) x (# conditions) x (# bootstrap iterations) \r\n" +
                        "[Figure 4A] - Cumulative distribution function model: \r\n" +
                        "   Total iterations = (# conditions) x (# bootstrap iterations) \r\n" +
                        "=================================================================================\r\n" +
                        "Bootstrap iterations per condition: 5");
                out.flush();
                out.close();
            }
            /* ************************************************************************************************************
             * if exists, attempt to read and scan contents
             * if bootstrap iteration numbers are not in an expected format (i.e. integer) following exceptions are handled
             *   1. double value -> convert to integer
             *   2. string -> error msg and continue with default setting
             *   3. mixed -> error msg and continue with default setting
             *   4. empty lines and other junks -> the program only seeks for a line with ":" otherwise all ignored
             **************************************************************************************************************/
            else {
                FileInputStream fis = new FileInputStream(textfile);
                Scanner in = new Scanner(fis);

                while (in.hasNextLine()) {
                    fileExists = true;
                    String currentLine = in.nextLine();
                    if (currentLine.contains(":")) {
                        String words[] = currentLine.split(":");

                        for (String str : words) {
                            try {
                                foundBootIter = Double.parseDouble(str);
                                // Negative number exception -> simply report to the user and use default setting
                                if (foundBootIter < 0) {
                                    text.append("     ***Found negative bootstrap iteration number! using default setting...\n");
                                    foundBootIter = defaultBootIter;
                                }
                                // Notice that array declarations in DataManger for 95% CIs requires non-zero size of array. So guide user to have at least 1 iteration per condition.
                                else if (foundBootIter == 0) {
                                    text.append("     ***Found 0 bootstrap iteration number! The program needs at least 1 iteration!\n");
                                    foundBootIter = 1;
                                }
                                // foundNumber would not even reach this if the loop passes string varaible
                                foundNumber = true;
                            } catch (NumberFormatException nfe) {
                                // string is not a number
                            }
                        }
                    }
                }
                in.close();
                fis.close();
            }
        } catch(IOException e){
            e.printStackTrace();
        }

        // provide feedback to the user
        if (fileExists && foundNumber) {
            this.numIter = (int) foundBootIter;
            text.append("...Discrete Gaussian : Bootstrapping...");
        }
        else if (fileExists && !foundNumber) {
            this.numIter = defaultBootIter;
            text.append("     ***Found 'settings.txt' file, but no iteration numbers specified! using default setting...\n");
            text.append("...Discrete Gaussian : Bootstrapping...");
        }
        else {
            this.numIter = defaultBootIter;
            text.append("...Discrete Gaussian : Bootstrapping...");
        }


        // BOOTSTRAP BEGINS
        randGen.setSeed(38298723);

        this.AStar = new double[DATA.numTimePoints][DATA.numConditions][numIter];
        this.muStar = new double[DATA.numTimePoints][DATA.numConditions][numIter];
        this.sigmaStar = new double[DATA.numTimePoints][DATA.numConditions][numIter];

        this.AmpCDFStar = new double[DATA.numConditions][numIter];
        this.muCDFStar = new double[DATA.numConditions][numIter];
        this.sigmaCDFStar = new double[DATA.numConditions][numIter];

        double[][] DataPerGen = new double[DATA.lastGeneration][DATA.numReplicates];

        int progressCounter = 0;
        double totalProcess = DATA.numTimePoints*DATA.numConditions*numIter;

        for (int itpt = 0; itpt < DATA.numTimePoints; itpt++) {
            for (int icnd = 0; icnd < DATA.numConditions; icnd++) {
                // Check if there's an empty rows (i.e. condition)
                if (!missingStat[itpt][icnd]) {
                    for (int biter = 0; biter < numIter; biter++) {

                        int counter = 0;
                        int countZeros = 0;
                        int xCordMax = 0;
                        double max = 0.0;
                        double min = Double.MAX_VALUE;

                        for (int igen = 1; igen <= DATA.lastGeneration; igen++) {
                            for (int irep = 0; irep < numRepDataCounter[itpt][icnd][igen - 1]; irep++) {
                                DataPerGen[igen - 1][irep] = origY[itpt][icnd][igen - 1][irep];
                                counter++;
                            }
                        }

                        double[] RESAMPLED = Resample(DATA, itpt, icnd, DataPerGen, numRepDataCounter, counter);
                        for (int count = 0; count < counter; count++) {
                            if (RESAMPLED[count] > max) {
                                max = RESAMPLED[count];
                                xCordMax = count;
                            }

                            if (RESAMPLED[count] < min) {
                                min = RESAMPLED[count];
                            }

                            if (RESAMPLED[count] == 0.0) {
                                countZeros++;
                            }
                        }

                        Object[] CARUANA = FitData.CaruanaFit(origX[itpt][icnd], RESAMPLED);
                        DataTable reEstimate = FitData.lmGaussFit(DATA, itpt, origX[itpt][icnd], RESAMPLED, counter, countZeros, xCordMax, max, min, (double) CARUANA[0], (double) CARUANA[1], (double) CARUANA[2]);

                        double A = Double.valueOf(reEstimate.get(2, 0).toString());
                        double mu = Double.valueOf(reEstimate.get(3, 0).toString());
                        double sigma = Double.valueOf(reEstimate.get(4, 0).toString());

                        AStar[itpt][icnd][biter] = A;
                        muStar[itpt][icnd][biter] = mu;
                        sigmaStar[itpt][icnd][biter] = sigma;

                        progressCounter++;
                        Border titledBorder = BorderFactory.createTitledBorder("Discrete Gaussian : Bootstrapping...["+progressCounter+"/"+(int)totalProcess+"]");
                        double currProcess = 100 * (progressCounter) / totalProcess;
                        progressBar.setBorder(titledBorder);
                        progressBar.setValue((int) currProcess);
                    }
                }
                else {
                    progressCounter += numIter;
                    Border titledBorder = BorderFactory.createTitledBorder("Discrete Gaussian :  Bootstrapping...["+progressCounter+"/"+(int)totalProcess+"]");
                    double currProcess = 100 * (progressCounter) / totalProcess;
                    progressBar.setBorder(titledBorder);
                    progressBar.setValue((int) currProcess);
                }

            }
        }

        text.append("done\n");
    }

    private double[] Resample(Dataset DATA, int itpt, int icnd, double[][] targetData, int[][][] numRepDataCounter, int numData) {

        double[] resData = new double[numData];

        int counter = 0;
        for (int igen = 1; igen <= DATA.lastGeneration; igen++) {
            for (int irep = 0; irep < numRepDataCounter[itpt][icnd][igen-1]; irep++) {
                resData[counter] = targetData[igen-1][randGen.nextInt(numRepDataCounter[itpt][icnd][igen-1])];
                counter++;
            }
        }

        return resData;
    }

    final private double[][] AmpCDFStar;
    final private double[][] muCDFStar;
    final private double[][] sigmaCDFStar;

    public void BootstrapCDF(Dataset DATA,
                             double[][] origX, double[][][] origY,
                             boolean[][] missingStat,
                             int[][] numRepDataCounter,
                             JTextArea text,
                             JProgressBar progressBar) {

        int iters = getBootIter();

        text.append("...Cumulative Distribution : Bootstrapping...");
        int progressCounter = 0;
        double totalProcess = DATA.numConditions*iters;

        double[][] DataPerTP = new double[DATA.numTimePoints][DATA.numReplicates];
        for (int icnd = 0; icnd < DATA.numConditions; icnd++) {
            for (int biter = 0; biter < iters; biter++) {

                int dataCounter = 0;
                int countZeros = 0;
                int xCordMax = 0;
                double max = 0.0;
                double min = Double.MAX_VALUE;

                for (int itpt = 0; itpt < DATA.numTimePoints; itpt++) {

                    if (!missingStat[itpt][icnd]) {
                        for (int irep = 0; irep < numRepDataCounter[itpt][icnd]; irep++) {
                            DataPerTP[itpt][irep] = origY[itpt][icnd][irep];
                            dataCounter++;
                        }
                    }
                }

                double[] RESAMPLED = Resample2(DATA, icnd, DataPerTP, numRepDataCounter, dataCounter);

                for (int count = 0; count < dataCounter; count++) {
                    if (RESAMPLED[count] > max) {
                        max = RESAMPLED[count];
                        xCordMax = count;
                    }

                    if (RESAMPLED[count] < min) {
                        min = RESAMPLED[count];
                    }

                    if (RESAMPLED[count] == 0.0) {
                        countZeros++;
                    }
                }

                Object[] reEstimate = FitData.lmCDF(origX[icnd], RESAMPLED, dataCounter, countZeros, xCordMax, max, min);

                DataTable GAUSS = (DataTable) reEstimate[0];
//                DataTable CDF = (DataTable) reEstimate[1];

                double Amp = Double.valueOf(GAUSS.get(2, 0).toString());
                double mu = Double.valueOf(GAUSS.get(3, 0).toString());
                double sigma = Double.valueOf(GAUSS.get(4, 0).toString());

                AmpCDFStar[icnd][biter] = Amp;
                muCDFStar[icnd][biter] = mu;
                sigmaCDFStar[icnd][biter] = sigma;

                progressCounter++;
                Border titledBorder = BorderFactory.createTitledBorder("Cumulative Distribution : Bootstrapping...["+progressCounter+"/"+(int)totalProcess+"]");
                progressBar.setBorder(titledBorder);
                double currProgress = 100*(progressCounter/totalProcess);
                progressBar.setValue((int)currProgress);
            }
        }

        text.append("done\n");
    }

    private double[] Resample2(Dataset DATA, int icnd, double[][] targetData, int[][] numRepDataCounter, int numData) {

        double[] resData = new double[numData];

        int dataCounter = 0;
        for (int itpt = 0; itpt < DATA.numTimePoints; itpt++) {
            for (int irep = 0; irep < numRepDataCounter[itpt][icnd]; irep++) {
                resData[dataCounter] = targetData[itpt][randGen.nextInt(numRepDataCounter[itpt][icnd])];
                dataCounter++;
            }
        }

        return resData;
    }

    public int getBootIter() { return numIter; }

    public double[][][] getAstar() { return AStar; }
    public double[][][] getMuStar() { return muStar; }
    public double[][][] getSigmaStar() { return sigmaStar; }

    public double[][] getAmpCDFStar() { return AmpCDFStar; }
    public double[][] getMuCDFStar() { return muCDFStar; }
    public double[][] getSigmaCDFStar() { return sigmaCDFStar; }
}
