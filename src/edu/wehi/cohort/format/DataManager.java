/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.wehi.cohort.format;


import java.io.File;
import java.lang.Object;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.text.DecimalFormat;
import javax.swing.JTextArea;
import javax.swing.JProgressBar;
import javax.swing.BorderFactory;
import javax.swing.border.Border;

import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.data.DataSeries;
import de.erichseifert.gral.data.DataSource;
import de.erichseifert.gral.data.comparators.Ascending;
import org.ejml.simple.SimpleMatrix;

import edu.wehi.cohort.gui.Dataset;

/**
 *
 * @author akan
 */
public class DataManager {

    final public int numConditions;
    final public int numReplicates;
    final public int numTimePoints;
    final public int lastGen;
    
    // Source data for plots (simply 'add' to the XYPlot object)
    final public List<List <DataSource>> CohortvTime;
    final public List<List <DataSource>> fittedGaussCurve;

    final public List<DataSource> MDNvTime;
    final public List<DataSource> fittedPieceWise;
    final public List<DataSource> AmplitudevTime;
    final public List<DataSource> FittedCrash;

    final public List<DataSource> SigmaTime;
    final public List<DataSource> MuSigma;
    final public List<DataSource> gen0;

    final public List<DataSource> CellvTime;
    final public List<DataSource> ArithmeticMDNvTime;
    final public List<DataSource> ArithmeticCohortvTime;
    final public List<DataSource> ArithmeticCrash;

    final public List<DataSource> ArithmeticCohortProportion;

    final public List<DataSource> fittedCohortProportionCDF;
    final public List<DataSource> complementaryGaussCurve;

    /* Arrys of fitted parameters via Caruana and Levenberg-Marquardt algorithms
     * Caruana fit is only used for Gaussian fit as it does NOT converge for Logistic func
     * 1st dim: time point
     * 2nd dim: conditions
     * 3rd dim:
     *    if GaussianFit: [0]=Amplitude   [1]=Mu   [2]=Sigma
     *    if LogisticFit: [0]=A   [1]=B   [2]=C   [3]=D
     * 
     * NB: carParams[itpt][icnd][3] - 3 parameters for Gaussian
     *     lmParamsGauss[itpt][icnd][3] - 3 parameters for Gaussian
     *     lmParamsLogistic[icnd][4] - 4 parameters for Logistic
     */
    final public double[][][] carParams;
    final public double[][][] lmParamsGauss;
    final public double[][] lmParamsPieceWise;
    final public double[][] lmCDF;

    final public int BOOT_ITERATIONS;

    final public double[][] lowA;
    final public double[][] highA;
    final public double[][] lowMu;
    final public double[][] highMu;
    final public double[][] lowSigma;
    final public double[][] highSigma;

    final public double[] lowtDD;
    final public double[] hightDD;
    final public double[] lowDD;
    final public double[] highDD;
    final public double[] lowSlope;
    final public double[] highSlope;
    final public double[] lowtFirstDiv;
    final public double[] hightFirstDiv;

    final public double[] lowAmpCDF;
    final public double[] highAmpCDF;
    final public double[] lowMuCDF;
    final public double[] highMuCDF;
    final public double[] lowSigmaCDF;
    final public double[] highSigmaCDF;
    
    public DataManager(Dataset DATA,
                       JTextArea text,
                       JProgressBar progressBar,
                       File inputFile) {

        this.numConditions = DATA.numConditions;
        this.numReplicates = DATA.numReplicates;
        this.numTimePoints = DATA.numTimePoints;
        this.lastGen = DATA.lastGeneration;

        this.CellvTime = new ArrayList<>();
        this.ArithmeticCohortvTime = new ArrayList<>();
        this.ArithmeticMDNvTime = new ArrayList<>();
        this.ArithmeticCrash = new ArrayList<>();
        this.ArithmeticCohortProportion = new ArrayList<>();

        this.fittedCohortProportionCDF = new ArrayList<>();
        this.complementaryGaussCurve = new ArrayList<>();

        this.CohortvTime = new ArrayList<>();
        this.fittedGaussCurve = new ArrayList<>();
        this.MDNvTime = new ArrayList<>();
        this.fittedPieceWise = new ArrayList<>();
        this.AmplitudevTime = new ArrayList<>();
        this.SigmaTime = new ArrayList<>();
        this.MuSigma = new ArrayList<>();
        this.FittedCrash = new ArrayList<>();

        // Fitted Parameters
        this.carParams = new double[numTimePoints][numConditions][3];
        this.lmParamsGauss = new double[numTimePoints][numConditions][3];
        this.lmParamsPieceWise = new double[numConditions][4];
        this.lmCDF = new double[numConditions][3];


        SimpleMatrix timeStamps = DATA.getTimeStamps();
        List<String> condNames = DATA.getConditionNames();
        SimpleMatrix cohorts = DATA.getCohorts();
        boolean[][] missingStat = DATA.getmissingStat();
        boolean[][][] missingStatReps = DATA.getMissingStatreps();

        /* ********************************************************************
        * Load data form Dataset object 'DATA' & sort them for fittings later on
        * *********************************************************************/
        text.append("\n...Sorting data...");
        Object[][] sortedDATA = loadPreComputedData(DATA, condNames, timeStamps, cohorts, missingStat, missingStatReps);

        // Cohort per generation : sorted dataset for discrete Gaussian fit (Fig2)
        int[][] XCordMax = (int[][])sortedDATA[0][6];
        int[][] CountZeros = (int[][])sortedDATA[0][5];
        int[][] NumDataCounter = (int[][])sortedDATA[0][3];
        int[][][] NumRepDataCounter = (int[][][])sortedDATA[0][4];
        double[][] MaxDATA = (double[][])sortedDATA[0][7];
        double[][] MinDATA = (double[][])sortedDATA[0][8];
        double[][][] XData = (double[][][])sortedDATA[0][0];
        double[][][] YData = (double[][][])sortedDATA[0][1];
        double[][][][] YDataRep = (double[][][][])sortedDATA[0][2];

        // Proportion of total cohort number (excl. div 0) for CDF fit (Fig4)
        int[] XCordMaxCDF = (int[])sortedDATA[1][6];
        int[] CountZerosCDF = (int[])sortedDATA[1][5];
        int[] NumDataCounterCDF = (int[])sortedDATA[1][3];
        int[][] NumRepDataCounterCDF = (int[][])sortedDATA[1][4];
        double[] MaxDataCDF = (double[])sortedDATA[1][7];
        double[] MinDataCDF = (double[])sortedDATA[1][8];
        double[][] XDataCDF = (double[][])sortedDATA[1][0];
        double[][] YDataCDF = (double[][])sortedDATA[1][1];
        double[][][] YDataRepCDF = (double[][][])sortedDATA[1][2];
        text.append("done\n");

        /* *****************************
        * Run first bootstrap (Figure 2)
        * ******************************/
        this.lowA = new double[numTimePoints][numConditions];
        this.highA = new double[numTimePoints][numConditions];
        this.lowMu = new double[numTimePoints][numConditions];
        this.highMu = new double[numTimePoints][numConditions];
        this.lowSigma = new double[numTimePoints][numConditions];
        this.highSigma = new double[numTimePoints][numConditions];

        Bootstrap BOOT = new Bootstrap(DATA, XData, YDataRep, missingStat, NumRepDataCounter, text, progressBar, inputFile);

        this.BOOT_ITERATIONS = BOOT.getBootIter();

        double CI = 0.05;
        double[][][] ASTAR = BOOT.getAstar();
        double[][][] MUSTAR = BOOT.getMuStar();
        double[][][] SIGMASTAR = BOOT.getSigmaStar();


        /* *********************************************
        * Fitting section 1 : Discrete Gaussian Function
        * **********************************************/
        int progressCounter = 0;
        double totalProcess = numTimePoints*numConditions;
        Border titledBorder1 = BorderFactory.createTitledBorder("Fitting discrete Gaussian model...");
        progressBar.setBorder(titledBorder1);
        text.append("...Fitting discrete Gaussian model...");

        for (int itpt = 0; itpt < numTimePoints; itpt++) {

            // Initiate discrete Gaussian fitting
            List<DataSource> LM_DISCRETE_GAUSS_OUTPUT = new ArrayList<>();
            for (int icnd = 0; icnd < numConditions; icnd++) {

                if (!missingStat[itpt][icnd]) {
                    // Run Caruana's algorithm
                    Object[] CARUANA = FitData.CaruanaFit(XData[itpt][icnd], YData[itpt][icnd]);

                    // Store Caruana parameters
                    carParams[itpt][icnd][0] = (double) CARUANA[0];
                    carParams[itpt][icnd][1] = (double) CARUANA[1];
                    carParams[itpt][icnd][2] = (double) CARUANA[2];

                    // Run Levenberg-Marquardt algorithm with Caruana's parameters as initial guesses
                    DataTable LM = FitData.lmGaussFit(DATA, itpt,
                            XData[itpt][icnd], YData[itpt][icnd],
                            NumDataCounter[itpt][icnd], CountZeros[itpt][icnd],
                            XCordMax[itpt][icnd], MaxDATA[itpt][icnd], MinDATA[itpt][icnd],
                            carParams[itpt][icnd][0],
                            carParams[itpt][icnd][1],
                            carParams[itpt][icnd][2]);

                    // Store LM params in double[][][] lmParamsGauss array for MDN plots
                    double A = Double.valueOf(LM.get(2, 0).toString());
                    double mu = Double.valueOf(LM.get(3, 0).toString());
                    double sigma = Double.valueOf(LM.get(4, 0).toString());
                    lmParamsGauss[itpt][icnd][0] = A;
                    lmParamsGauss[itpt][icnd][1] = mu;
                    lmParamsGauss[itpt][icnd][2] = sigma;

                    // Get bootstrap re-estimated parameters for plotting & 95% confidence interval
                    DataTable Astar = new DataTable(Double.class);
                    DataTable Mustar = new DataTable(Double.class);
                    DataTable Sigmastar = new DataTable(Double.class);
                    for (int bootIDX = 0; bootIDX < BOOT_ITERATIONS; bootIDX++) {
                        Astar.add(ASTAR[itpt][icnd][bootIDX]);
                        Mustar.add(MUSTAR[itpt][icnd][bootIDX]);
                        Sigmastar.add(SIGMASTAR[itpt][icnd][bootIDX]);
                    }

                    // Sort in ascending order to cut 95% above and below the resampled distributions
                    Astar.sort(new Ascending(0));
                    Mustar.sort(new Ascending(0));
                    Sigmastar.sort(new Ascending(0));

                    double AlowQ = (CI / 2.0) * Astar.getRowCount();
                    double AhighQ = (1.0 - CI / 2.0) * Astar.getRowCount();
                    lowA[itpt][icnd] = Double.valueOf(Astar.get(0, (int) AlowQ).toString());
                    highA[itpt][icnd] = Double.valueOf(Astar.get(0, (int) AhighQ).toString());

                    double MulowQ = (CI / 2.0) * Mustar.getRowCount();
                    double MuhighQ = (1.0 - CI / 2.0) * Mustar.getRowCount();
                    lowMu[itpt][icnd] = Double.valueOf(Mustar.get(0, (int) MulowQ).toString());
                    highMu[itpt][icnd] = Double.valueOf(Mustar.get(0, (int) MuhighQ).toString());

                    double SigmalowQ = (CI / 2.0) * Sigmastar.getRowCount();
                    double SigmahighQ = (1.0 - CI / 2.0) * Sigmastar.getRowCount();
                    lowSigma[itpt][icnd] = Double.valueOf(Sigmastar.get(0, (int) SigmalowQ).toString());
                    highSigma[itpt][icnd] = Double.valueOf(Sigmastar.get(0, (int) SigmahighQ).toString());

                    // beautify the plot legends!
                    DecimalFormat df = new DecimalFormat("#0.0");
                    if (Double.isNaN(mu) && Double.isNaN(sigma)) {
                        DataSeries fittedData = new DataSeries("\u00B5=NaN\n" + "\u03C3=NaN", LM);
                        LM_DISCRETE_GAUSS_OUTPUT.add(fittedData);
                    } else if (Double.isNaN(mu)) {
                        DataSeries fittedData = new DataSeries("\u00B5=NaN\n" + "\u03C3=" + df.format(lmParamsGauss[itpt][icnd][2]), LM);
                        LM_DISCRETE_GAUSS_OUTPUT.add(fittedData);
                    } else if (Double.isNaN(sigma)) {
                        DataSeries fittedData = new DataSeries("\u00B5=" + df.format(lmParamsGauss[itpt][icnd][1]) + "\n" + "\u03C3=NaN", LM);
                        LM_DISCRETE_GAUSS_OUTPUT.add(fittedData);
                    } else {
                        // TODO: figure out a better way to put CI in the legend
                        if (lmParamsGauss[itpt][icnd][1] < 0 || lowMu[itpt][icnd] < 0) {
                            DataSeries fittedData = new DataSeries(
                                    "\u00B5=" + df.format(lmParamsGauss[itpt][icnd][1]) + "\n" +
                                            "(" + df.format(lowMu[itpt][icnd]) + ", " + df.format(highMu[itpt][icnd]) + ")\n" +
                                            "\u03C3=" + df.format(lmParamsGauss[itpt][icnd][2]) + "\n" +
                                            "(" + df.format(lowSigma[itpt][icnd]) + ", " + df.format(highSigma[itpt][icnd]) + ")", LM);
                            LM_DISCRETE_GAUSS_OUTPUT.add(fittedData);
                        } else {
                            DataSeries fittedData = new DataSeries(
                                    "\u00B5=" + df.format(lmParamsGauss[itpt][icnd][1]) +
                                            "(" + df.format(lowMu[itpt][icnd]) + ", " + df.format(highMu[itpt][icnd]) + ")\n" +
                                            "\u03C3=" + df.format(lmParamsGauss[itpt][icnd][2]) +
                                            "(" + df.format(lowSigma[itpt][icnd]) + ", " + df.format(highSigma[itpt][icnd]) + ")", LM);
                            LM_DISCRETE_GAUSS_OUTPUT.add(fittedData);
                        }
                    }
                }
                progressCounter++;
                double currProcess = 100*(progressCounter)/totalProcess;
                progressBar.setValue((int)currProcess);
            }
            fittedGaussCurve.add(LM_DISCRETE_GAUSS_OUTPUT);
        }
        text.append("done\n");

        /* ************************************************************
         * Fitting section 2 : Piece-wise Linear Functions
         * Call bootstrap iteration number
         *   -> construct [MDN vs. Time] arrays per bootstrap iteration
         *   -> run Piece-wise fit for each of it
         *   -> compute 95% CI
         * *************************************************************/
        progressCounter = 0;
        double totalProcess2 = (numConditions*numTimePoints*BOOT_ITERATIONS + numConditions*numTimePoints);
        Border titledBorder2 = BorderFactory.createTitledBorder("Fitting piece-wise linear model...");
        progressBar.setBorder(titledBorder2);
        text.append("...Fitting piece-wise linear model...");

        double[] bootTimePer = new double[numTimePoints];
        double[] bootMDNPer = new double[numTimePoints];
        double[][] WEIGHT = new double[numConditions][numTimePoints];

        double[][] tDDPer = new double[numConditions][BOOT_ITERATIONS];
        double[][] DDPer = new double[numConditions][BOOT_ITERATIONS];
        double[][] slopePer = new double[numConditions][BOOT_ITERATIONS];
        double[][] tFirstDiv = new double[numConditions][BOOT_ITERATIONS];

        this.lowtDD = new double[numConditions];
        this.hightDD = new double[numConditions];
        this.lowDD = new double[numConditions];
        this.highDD = new double[numConditions];
        this.lowSlope = new double[numConditions];
        this.highSlope = new double[numConditions];
        this.lowtFirstDiv = new double[numConditions];
        this.hightFirstDiv = new double[numConditions];

        // Compute WEIGHT for mu per time point
        double[] tempMDN = new double[numTimePoints*BOOT_ITERATIONS];
        for (int icnd = 0; icnd < numConditions; icnd++) {
            int counter = 0;
            int counter2 = 0;
            for (int itpt = 0; itpt < numTimePoints; itpt++) {
                if (!missingStat[itpt][icnd]) {
                    for (int biter = 0; biter < BOOT_ITERATIONS; biter++) {
                        tempMDN[counter] = MUSTAR[itpt][icnd][biter];
                        counter++;
                    }

                    double sum = 0.0;
                    for (int i = 0; i < BOOT_ITERATIONS; i++) {
                        sum += tempMDN[counter - i - 1];
                    }
                    double mean = sum / BOOT_ITERATIONS;
                    double tss = 0.0;
                    for (int i = 0; i < BOOT_ITERATIONS; i++) {
                        tss += Math.pow(tempMDN[counter - i - 1] - mean, 2);
                    }

                    if (tss == 0.0) {
                        WEIGHT[icnd][counter2] = 1.0;
                    } else {
                        WEIGHT[icnd][counter2] = 1.0 / Math.sqrt((tss / (BOOT_ITERATIONS - 1)));
                    }
                    counter2++;
                }
            }
        }

        // Compute 95% confidence interval
        for (int icnd = 0; icnd < numConditions; icnd++) {
            for (int biter = 0; biter < BOOT_ITERATIONS; biter++) {
                int dataCOUNTER = 0;
                for (int itpt = 0; itpt < numTimePoints; itpt++) {
                    if (!missingStat[itpt][icnd]) {
                        bootTimePer[dataCOUNTER] = timeStamps.get(itpt);
                        bootMDNPer[dataCOUNTER] = MUSTAR[itpt][icnd][biter];

                        dataCOUNTER++;
                    }
                    progressCounter++;
                    double currProgress = 100*(progressCounter/totalProcess2);
                    progressBar.setValue((int)currProgress);
                }
                DataTable PieceWise = FitData.lmPiecewise(bootTimePer, bootMDNPer, WEIGHT[icnd], dataCOUNTER);
                double timeDD = Double.valueOf(PieceWise.get(2, 0).toString());
                double yInters = Double.valueOf(PieceWise.get(3, 0).toString());
                double DR = Double.valueOf(PieceWise.get(4, 0).toString());

                tDDPer[icnd][biter] = timeDD;
                DDPer[icnd][biter] = yInters + DR*timeDD;
                slopePer[icnd][biter] = 1.0/DR;
                tFirstDiv[icnd][biter] = (1.0 - yInters)/DR;
            }
            Arrays.sort(tDDPer[icnd]);
            Arrays.sort(DDPer[icnd]);
            Arrays.sort(slopePer[icnd]);
            Arrays.sort(tFirstDiv[icnd]);

            double tDDLowQ = (CI/2.0)*tDDPer[icnd].length;
            double tDDHighQ = (1-CI/2.0)*tDDPer[icnd].length;

            double DDLowQ = (CI/2.0)*DDPer[icnd].length;
            double DDHighQ = (1-CI/2.0)*DDPer[icnd].length;

            double slopeLowQ = (CI/2.0)*slopePer[icnd].length;
            double slopeHighQ = (1-CI/2.0)*slopePer[icnd].length;

            double tFirstDivLowQ = (CI/2.0)*tFirstDiv[icnd].length;
            double tFirstDivHighQ = (1-CI/2.0)*tFirstDiv[icnd].length;

            lowtDD[icnd] = tDDPer[icnd][(int)tDDLowQ];
            hightDD[icnd] = tDDPer[icnd][(int)tDDHighQ];

            lowDD[icnd] = DDPer[icnd][(int)DDLowQ];
            highDD[icnd] = DDPer[icnd][(int)DDHighQ];

            lowSlope[icnd] = slopePer[icnd][(int)slopeLowQ];
            highSlope[icnd] = slopePer[icnd][(int)slopeHighQ];

            lowtFirstDiv[icnd] = tFirstDiv[icnd][(int)tFirstDivLowQ];
            hightFirstDiv[icnd] = tFirstDiv[icnd][(int)tFirstDivHighQ];
        }

        int[] NaNcounter = new int[numConditions];
        double[] origTime = new double[numTimePoints];
        double[] origMDN = new double[numTimePoints];
        double[] origSIGMA = new double[numTimePoints];
//        double[] bootTime = new double[numTimePoints*BOOT_ITERATIONS];
//        double[] bootMDN = new double[numTimePoints*BOOT_ITERATIONS];
        for (int icnd = 0; icnd < numConditions; icnd++) {

            DataTable MDN = new DataTable(Double.class, Double.class, Double.class, Double.class);
            DataTable SIGMATIME = new DataTable(Double.class, Double.class);
            DataTable MUSIGMA = new DataTable(Double.class, Double.class);
            DataTable AmplitudeData = new DataTable(Double.class, Double.class);
            DataTable fitCRASH = new DataTable(Double.class, Double.class);

            int dataCOUNTER = 0;
            for (int itpt = 0; itpt < numTimePoints; itpt++) {
                if (!missingStat[itpt][icnd]) {
                    double time = timeStamps.get(itpt);
                    double A = lmParamsGauss[itpt][icnd][0];
                    double mu = lmParamsGauss[itpt][icnd][1];
                    double sigma = lmParamsGauss[itpt][icnd][2];

                    // Bring all bootrapped parameters
//                for (int biter = 0; biter < BOOT_ITERATIONS; biter++) {
//                    bootTime[counter] = time;
//                    bootMDN[counter] = MUSTAR[itpt][icnd][biter];
//                    counter++;
//
//                    progressCounter++;
//                    double currProgress = 100*(progressCounter/totalProcess2);
//                    progressBar.setValue((int)currProgress);
//                }
                    origTime[dataCOUNTER] = time;
                    origMDN[dataCOUNTER] = mu;
                    origSIGMA[dataCOUNTER] = sigma;

                    double ERRTOP = Math.abs(highMu[itpt][icnd] - origMDN[dataCOUNTER]);
                    double ERRBOT = Math.abs(origMDN[dataCOUNTER] - lowMu[itpt][icnd]);

                    MDN.add(origTime[dataCOUNTER], origMDN[dataCOUNTER], ERRBOT, ERRTOP);
                    SIGMATIME.add(origTime[dataCOUNTER], Math.pow(origSIGMA[dataCOUNTER], 2));
                    MUSIGMA.add(origMDN[dataCOUNTER], Math.pow(origSIGMA[dataCOUNTER], 2));

                    // NaN counter for Log-norm fit which determines dimension of useful data
                    double pred = FitData.DGaussFunc(true, false, 0, A, mu, sigma);
                    if (!Double.isNaN(pred)) {
                        NaNcounter[icnd] += 1;
                    }

                    AmplitudeData.add(origTime[dataCOUNTER], A * Math.sqrt(2 * Math.PI * Math.pow(sigma, 2)));
                    fitCRASH.add(origMDN[dataCOUNTER], A * Math.sqrt(2 * Math.PI * Math.pow(sigma, 2)));

                    dataCOUNTER++;
                }

                progressCounter++;
                double currProgress = 100*(progressCounter/totalProcess2);
                progressBar.setValue((int)currProgress);
            }

            DataSeries finalMDN = new DataSeries(condNames.get(icnd), MDN);
            MDNvTime.add(finalMDN);

            DataTable PieceWise = FitData.lmPiecewise(origTime, origMDN, WEIGHT[icnd], dataCOUNTER);
            double timeDD = Double.valueOf(PieceWise.get(2, 0).toString());  // breakpoint
            double yInters = Double.valueOf(PieceWise.get(3, 0).toString()); // y-intercept
            double DR = Double.valueOf(PieceWise.get(4, 0).toString());      // b
            lmParamsPieceWise[icnd][0] = timeDD;                // time to division destiny
            lmParamsPieceWise[icnd][1] = yInters + DR*timeDD;   // division desitny
            lmParamsPieceWise[icnd][2] = 1.0/DR;                // time to enter k to k+1 division
            lmParamsPieceWise[icnd][3] = (1.0 - yInters)/DR;

            DecimalFormat df = new DecimalFormat("#0.0");
            DataSeries finalPieceWise = new DataSeries(
                    "m=" + df.format(1.0/DR) +
                            " (" + df.format(lowSlope[icnd]) + ", " +
                            df.format(highSlope[icnd]) + ")" + "\n" +
                    "DD=" + df.format(yInters + DR*timeDD) +
                            " ("+df.format(lowDD[icnd]) + ", " +
                            df.format(highDD[icnd]) + ")"+ "\n" +
                    "tFD=" + df.format((1.0 - yInters)/DR) +
                            " ("+df.format(lowtFirstDiv[icnd]) + ", " +
                            df.format(hightFirstDiv[icnd]) + ")" + "\n" +
                        "tDD=" + df.format(timeDD) +
                            " ("+df.format(lowtDD[icnd]) + ", " +
                             df.format(hightDD[icnd]) + ")", PieceWise);
            fittedPieceWise.add(finalPieceWise);

            DataSeries finalSIGMATIME = new DataSeries(condNames.get(icnd), SIGMATIME);
            SigmaTime.add(finalSIGMATIME);
            DataSeries finalMUSIGMA = new DataSeries(condNames.get(icnd), MUSIGMA);
            MuSigma.add(finalMUSIGMA);

            DataSeries AmplitudeDataz = new DataSeries(condNames.get(icnd), AmplitudeData);
            AmplitudevTime.add(AmplitudeDataz);

            DataSeries fitCRASHz = new DataSeries(condNames.get(icnd), fitCRASH);
            FittedCrash.add(fitCRASHz);
        }
        text.append("done\n");
        this.gen0 = loadCohortGen0(DATA, missingStat, missingStatReps);


        /* *******************************
        * Run second bootstrap (Figure 4A)
        * ********************************/
        BOOT.BootstrapCDF(DATA, XDataCDF, YDataRepCDF, missingStat, NumRepDataCounterCDF, text, progressBar);

        progressCounter = 0;
        double totalProcess3 = numConditions;
        Border titledBorder3 = BorderFactory.createTitledBorder("Fitting CDF model...");
        progressBar.setBorder(titledBorder3);
        text.append("...Fitting cumulative distribution model...");

        this.lowAmpCDF = new double[numConditions];
        this.highAmpCDF = new double[numConditions];
        this.lowMuCDF = new double[numConditions];
        this.highMuCDF = new double[numConditions];
        this.lowSigmaCDF = new double[numConditions];
        this.highSigmaCDF = new double[numConditions];

        double[][] AmpCDF = BOOT.getAmpCDFStar();
        double[][] muCDF = BOOT.getMuCDFStar();
        double[][] sigmaCDF = BOOT.getSigmaCDFStar();

        /* **********************************************
        * Fitting section 3 : CDF fit (Figure 4A)
        * ************************************************/
        for (int icnd = 0; icnd < numConditions; icnd++){
            Arrays.sort(AmpCDF[icnd]);
            Arrays.sort(muCDF[icnd]);
            Arrays.sort(sigmaCDF[icnd]);

            double lowAmpCDFQ = (CI/2.0)*AmpCDF[icnd].length;
            double highAmpCDFQ = (1-CI/2.0)*AmpCDF[icnd].length;

            double lowMuCDFQ = (CI/2.0)*muCDF[icnd].length;
            double highMuCDFQ = (1-CI/2.0)*muCDF[icnd].length;

            double lowSigmaCDFQ = (CI/2.0)*sigmaCDF[icnd].length;
            double highSigmaCDFQ = (1-CI/2.0)*sigmaCDF[icnd].length;

            lowAmpCDF[icnd] = AmpCDF[icnd][(int)lowAmpCDFQ];
            highAmpCDF[icnd] = AmpCDF[icnd][(int)highAmpCDFQ];

            lowMuCDF[icnd] = muCDF[icnd][(int)lowMuCDFQ];
            highMuCDF[icnd] = muCDF[icnd][(int)highMuCDFQ];

            lowSigmaCDF[icnd] = sigmaCDF[icnd][(int)lowSigmaCDFQ];
            highSigmaCDF[icnd] = sigmaCDF[icnd][(int)highSigmaCDFQ];

            Object[] RESULT = FitData.lmCDF(XDataCDF[icnd], YDataCDF[icnd], NumDataCounterCDF[icnd], CountZerosCDF[icnd], XCordMaxCDF[icnd], MaxDataCDF[icnd], MinDataCDF[icnd]);

            DataTable LM_GAUSS_CURVE = (DataTable) RESULT[0];
            DataTable LM_CDF_CURVE = (DataTable) RESULT[1];

            double Amp = Double.valueOf(LM_GAUSS_CURVE.get(2, 0).toString());
            double mu = Double.valueOf(LM_GAUSS_CURVE.get(3, 0).toString());
            double sigma = Double.valueOf(LM_GAUSS_CURVE.get(4, 0).toString());

            lmCDF[icnd][0] = Amp;
            lmCDF[icnd][1] = mu;
            lmCDF[icnd][2] = sigma;

            DecimalFormat df = new DecimalFormat("#0.0");
//            DecimalFormat df2 = new DecimalFormat("#0.00");
//            DataSeries fitted = new DataSeries("A=" + df2.format(Amp) +"\n" +
//                    "\u00B5=" + df.format(mu) + "\n" +
//                    "\u03C3=" + df.format(sigma), LM_GAUSS_CURVE);
            DataSeries fittedGauss = new DataSeries("\u00B5=" + df.format(mu) + " (" + df.format(lowMuCDF[icnd]) + ", " + df.format(highMuCDF[icnd]) +")\n" +
                    "\u03C3=" + df.format(sigma) + "(" + df.format(lowSigmaCDF[icnd]) + ", " + df.format(highSigmaCDF[icnd]) + ")", LM_GAUSS_CURVE);
            DataSeries fittedCDF = new DataSeries(condNames.get(icnd), LM_CDF_CURVE);

            complementaryGaussCurve.add(fittedGauss);
            fittedCohortProportionCDF.add(fittedCDF);

            progressCounter++;
            double currProgress = 100*(progressCounter/totalProcess3);
            progressBar.setValue((int)currProgress);
        }
        text.append("done\n");
    }

    private Object[][] loadPreComputedData(Dataset DATA,
                                         List<String> condNames,
                                         SimpleMatrix timeStamps,
                                         SimpleMatrix cohorts,
                                         boolean[][] missingStat,
                                         boolean[][][] missingStatReps) {

        int NdataCDF = numTimePoints*numReplicates;

        double[][] xDataCDF = new double[numConditions][NdataCDF];
        double[][] yDataCDF = new double[numConditions][NdataCDF];

        // this is for bootstrap resampling
        double[][][] yDataRepCDF = new double[numTimePoints][numConditions][numReplicates];
        int[][] numRepDATAcounterCDF = new int[numTimePoints][numConditions];

        // flagging missing data & rest is for initial guesses
        int[] numDATAcounterCDF = new int[numConditions];
        int[] countZerosCDF = new int[numConditions];
        int[] xCordMaxCDF = new int[numConditions];
        double[] maxDATACDF = new double[numConditions];
        double[] minDATACDF = new double[numConditions];

        for (int icnd = 0; icnd < numConditions; icnd++) {

            // Total Cell Number plot related
            DataTable TOTALCELLS = new DataTable(Double.class, Double.class, Double.class);

            // Arithmetic plots related
            DataTable COHORTSUMREPS = new DataTable(Double.class, Double.class, Double.class);
            DataTable MDNeach = new DataTable(Double.class, Double.class, Double.class);
            DataTable crashMDNCOHORT = new DataTable(Double.class, Double.class, Double.class);
            DataTable COHORTPROPS = new DataTable(Double.class, Double.class, Double.class);

            int dataIDX = 0;
            double max = 0;
            double min = Double.MAX_VALUE;

            for (int itpt = 0; itpt < numTimePoints; itpt++) {

                if (!missingStat[itpt][icnd]) {
                    // Simply number of cells
                    double numCELLS = 0;

                    // Arithmeticallu computed variables
                    double sumCOHORTSUMREPS = 0;
                    double sumMDNeach = 0;

                    double sumCOHORTProps = 0;

                    double arithMDN = DATA.meanDivNum[itpt][icnd];

                    int numRepCells = 0;
                    int numRepCohortSum = 0;
                    int numRepMDN = 0;
                    int numCohortProps = 0;
                    int numRepCDF = 0;

                    for (int irep = 0; irep < numReplicates; irep++) {

                        if (!missingStatReps[itpt][icnd][irep]) {
                            numCELLS += DATA.totalNumbers[itpt][icnd][irep];
                            numRepCells++;

                            sumCOHORTSUMREPS += DATA.cohortSumreps[itpt][icnd][irep];
                            numRepCohortSum++;

                            sumMDNeach += DATA.meanDiveach[itpt][icnd][irep];
                            numRepMDN++;

                            sumCOHORTProps += DATA.sumExcGen0reps[itpt][icnd][irep]/DATA.cohortSumreps[itpt][icnd][irep];
                            numCohortProps++;

                            double props = DATA.sumExcGen0reps[itpt][icnd][irep] / DATA.cohortSumreps[itpt][icnd][irep];
                            xDataCDF[icnd][dataIDX] = timeStamps.get(itpt);
                            yDataCDF[icnd][dataIDX] = props;
                            yDataRepCDF[itpt][icnd][numRepCDF] = props;

                            if (yDataCDF[icnd][dataIDX] > max) {
                                max = yDataCDF[icnd][dataIDX];
                                maxDATACDF[icnd] = max;
                                xCordMaxCDF[icnd] = dataIDX;
                            }

                            if (yDataCDF[icnd][dataIDX] < min) {
                                min = yDataCDF[icnd][dataIDX];
                                minDATACDF[icnd] = min;
                            }

                            if (yDataCDF[icnd][dataIDX] == 0.0) {
                                countZerosCDF[icnd] += 1;
                            }

                            numDATAcounterCDF[icnd] += 1;
                            numRepDATAcounterCDF[itpt][icnd] += 1;

                            dataIDX++;
                            numRepCDF++;
                        }
                    }

                    // TODO: error propagation for Crash plot -> horizontal error bars
                    // Compute statistics: average, total sum square(tss)
                    double meanCELLS = numCELLS / numRepCells;
                    double tssNUMCELLS = 0;
                    double avgCOHORTSUMREPS = sumCOHORTSUMREPS / numRepCohortSum;
                    double tssCOHORTSUMREPS = 0;
                    double avgMDNeach = sumMDNeach / numRepMDN;
                    double tssMDNeach = 0;
                    double avgCOHORTProps = sumCOHORTProps / numCohortProps;
                    double tssCOHORTProps = 0;
                    for (int irep = 0; irep < numReplicates; irep++) {
                        if (!missingStatReps[itpt][icnd][irep]) {
                            tssNUMCELLS += Math.pow(DATA.totalNumbers[itpt][icnd][irep] - meanCELLS, 2);
                            tssCOHORTSUMREPS += Math.pow(DATA.cohortSumreps[itpt][icnd][irep] - avgCOHORTSUMREPS, 2);
                            tssMDNeach += Math.pow(DATA.meanDiveach[itpt][icnd][irep] - avgMDNeach, 2);
                            tssCOHORTProps += Math.pow(DATA.sumExcGen0reps[itpt][icnd][irep]/DATA.cohortSumreps[itpt][icnd][irep] - avgCOHORTProps, 2);
                        }
                    }

                    // Compute standard error
                    double sigmaCELLS = 0;
                    double stderrCELLS = 0;
                    if (numRepCells > 1) {
                        sigmaCELLS = Math.sqrt(tssNUMCELLS / (numRepCells - 1));
                        stderrCELLS = sigmaCELLS / Math.sqrt(numRepCells);
                    }
                    else {
                        stderrCELLS = 0.0;
                    }
                    TOTALCELLS.add(timeStamps.get(itpt), meanCELLS, stderrCELLS);

                    double sigmaCOHORTSUMREPS = 0;
                    double stderrCOHORTSUMREPS = 0;
                    if (numRepCohortSum > 1) {
                        sigmaCOHORTSUMREPS = Math.sqrt(tssCOHORTSUMREPS / (numRepCohortSum - 1));
                        stderrCOHORTSUMREPS = sigmaCOHORTSUMREPS / Math.sqrt(numRepCohortSum);
                    }
                    else {
                        stderrCOHORTSUMREPS = 0.0;
                    }
                    COHORTSUMREPS.add(timeStamps.get(itpt), avgCOHORTSUMREPS, stderrCOHORTSUMREPS);
                    crashMDNCOHORT.add(arithMDN, avgCOHORTSUMREPS, stderrCOHORTSUMREPS);

                    double sigmaMDNeach = 0;
                    double stderrMDNeach = 0;
                    if (numRepMDN > 1) {
                        sigmaMDNeach = Math.sqrt(tssMDNeach / (numRepMDN - 1));
                        stderrMDNeach = sigmaMDNeach / Math.sqrt(numRepMDN);
                    }
                    else {
                        stderrMDNeach = 0.0;
                    }
                    MDNeach.add(timeStamps.get(itpt), avgMDNeach, stderrMDNeach);

                    double sigmaEx0Cohort = 0;
                    double stderrEx0Cohort = 0;
                    if (numCohortProps > 1) {
                        sigmaEx0Cohort = Math.sqrt(tssCOHORTProps / (numCohortProps - 1));
                        stderrEx0Cohort = sigmaEx0Cohort / Math.sqrt(numCohortProps);
                    }
                    else {
                        stderrEx0Cohort = 0.0;
                    }
                    COHORTPROPS.add(timeStamps.get(itpt), avgCOHORTProps, stderrEx0Cohort);
                }
            }

            // give each of dataset a label for legends
            DataSeries finalTOTALCELLS = new DataSeries(condNames.get(icnd), TOTALCELLS);
            DataSeries finalCOHORTSUMREPS = new DataSeries(condNames.get(icnd), COHORTSUMREPS);
            DataSeries finalMDNeach = new DataSeries(condNames.get(icnd), MDNeach);
            DataSeries finalCrashMDNCOHORT = new DataSeries(condNames.get(icnd), crashMDNCOHORT);
            DataSeries finalCOHORTPROPS = new DataSeries(condNames.get(icnd), COHORTPROPS);

            CellvTime.add(finalTOTALCELLS);
            ArithmeticCohortvTime.add(finalCOHORTSUMREPS);
            ArithmeticMDNvTime.add(finalMDNeach);
            ArithmeticCrash.add(finalCrashMDNCOHORT);
            ArithmeticCohortProportion.add(finalCOHORTPROPS);
        }

        /* Following loop returns :
         *   1. Average value of cohorts (totalCohort / # of replicates) per generation
         *   2. Standard error of the sample (symmetric top & bottom)
         */
        int Ndata = lastGen * numReplicates;

        double[][][] xData = new double[numTimePoints][numConditions][Ndata];
        double[][][] yData = new double[numTimePoints][numConditions][Ndata];

        // this is for bootstrap resampling
        double[][][][] yDataRep = new double[numTimePoints][numConditions][lastGen][numReplicates];

        // flagging missing data & rest is for initial guesses
        int[][] numDATAcounter = new int[numTimePoints][numConditions];
        int[][][] numRepDATAcounter = new int[numTimePoints][numConditions][lastGen];
        int[][] countZeros = new int[numTimePoints][numConditions];
        int[][] xCordMax = new int[numTimePoints][numConditions];
        double[][] maxDATA = new double[numTimePoints][numConditions];
        double[][] minDATA = new double[numTimePoints][numConditions];

        for (int itpt = 0; itpt < numTimePoints; itpt++) {

            List<DataSource> COHORT_OUTPUT = new ArrayList<>();

            for (int icnd = 0; icnd < numConditions; icnd++) {

                DataTable COHORT = new DataTable(Integer.class, Double.class, Double.class);

                if(!missingStat[itpt][icnd]) {
                    int IDX = 0;
                    double max = 0;
                    double min = Double.MAX_VALUE;

                    for (int igen = 0; igen <= lastGen; igen++) {
                        double totalCohort = 0;
                        int numRepCohort = 0;
                        for (int irep = 0; irep < numReplicates; irep++) {
                            if (!missingStatReps[itpt][icnd][irep]) {
                                int dataRow = igen + itpt * (lastGen + 1);
                                int dataCol = irep + icnd * numReplicates;

                                double cohort = cohorts.get(dataRow, dataCol);

                                totalCohort += cohort;

                                // Excluding generation 0 for fit analysis
                                if (igen > 0) {

                                    xData[itpt][icnd][IDX] = igen;
                                    yData[itpt][icnd][IDX] = cohort;
                                    yDataRep[itpt][icnd][igen-1][numRepCohort] = cohort;

                                    if (yData[itpt][icnd][IDX] > max) {
                                        max = yData[itpt][icnd][IDX];
                                        maxDATA[itpt][icnd] = max;
                                        xCordMax[itpt][icnd] = IDX;
                                    }

                                    if (yData[itpt][icnd][IDX] < min) {
                                        min = yData[itpt][icnd][IDX];
                                        minDATA[itpt][icnd] = min;
                                    }

                                    if (yData[itpt][icnd][IDX] == 0.0) {
                                        countZeros[itpt][icnd] += 1;
                                    }

                                    numDATAcounter[itpt][icnd] += 1;
                                    numRepDATAcounter[itpt][icnd][igen-1] += 1;

                                    IDX++;
                                }
                                numRepCohort++;
                            }
                        }

                        double meanCohort = totalCohort / numRepCohort;

                        // Compute standard error
                        double tssCOHORT = 0;
                        for (int irep = 0; irep < numReplicates; irep++) {
                            int dataRow = igen + itpt * (lastGen + 1);
                            int dataCol = irep + icnd * numReplicates;

                            if (!missingStatReps[itpt][icnd][irep]) {
                                tssCOHORT += Math.pow(cohorts.get(dataRow, dataCol) - meanCohort, 2);
                            }
                        }

                        double sigmaCOHORT = 0;
                        double stderrCOHORT = 0;
                        if (numRepCohort > 1) {
                            sigmaCOHORT = Math.sqrt(tssCOHORT / (numRepCohort - 1));
                            stderrCOHORT = sigmaCOHORT / Math.sqrt(numRepCohort);
                        }
                        else {
                            stderrCOHORT = 0.0;
                        }


                        COHORT.add(igen, meanCohort, stderrCOHORT);
                    }
                    DataSeries dataz = new DataSeries(condNames.get(icnd), COHORT);
                    COHORT_OUTPUT.add(dataz);
                }
            }
            CohortvTime.add(COHORT_OUTPUT);
        }

        return new Object[][]
            {
                // Cohort per generation : Discrete Gaussian fit data
                {xData, yData, yDataRep, numDATAcounter, numRepDATAcounter, countZeros, xCordMax, maxDATA, minDATA},
                // Prop. of Div0 cohort : CDF fit data
                {xDataCDF, yDataCDF, yDataRepCDF, numDATAcounterCDF, numRepDATAcounterCDF, countZerosCDF, xCordMaxCDF, maxDATACDF, minDATACDF}
            };
    }

    private List<DataSource> loadCohortGen0(Dataset DATA,
                                            boolean[][] missingStat,
                                            boolean[][][] missingStatRep) {

        List<DataSource> COHORT_GEN0_OUTPUT = new ArrayList<>();

        SimpleMatrix cohorts = DATA.getCohorts();
        SimpleMatrix timeStamps = DATA.getTimeStamps();
        List<String> condNames = DATA.getConditionNames();

        double igen = 0;
        double pred = 0;
        double[][][] Params = getlmParamsGauss();

        for (int icnd = 0; icnd < numConditions; icnd++) {

            DataTable resGen0 = new DataTable(Double.class, Double.class, Double.class);

            for (int itpt = 0; itpt < numTimePoints; itpt++) {

                if (!missingStat[itpt][icnd]) {
                    double residue = 0;

                    int dataCounter = 0;
                    for (int irep = 0; irep < numReplicates; irep++) {
                        if (!missingStatRep[itpt][icnd][irep]) {
                            int dataRow = (int) igen + itpt * (lastGen + 1);
                            int dataCol = irep + icnd * numReplicates;

                            double A = Params[itpt][icnd][0];
                            double mu = Params[itpt][icnd][1];
                            double sigma = Params[itpt][icnd][2];

                            pred = FitData.DGaussFunc(true, false, igen, A, mu, sigma);
                            residue += cohorts.get(dataRow, dataCol) - pred;
                            dataCounter++;
                        }
                    }

                    double meanResidue = residue / dataCounter;

                    // Compute standard error for samples of replicates
                    double diffSqr = 0;
                    for (int irep = 0; irep < numReplicates; irep++) {
                        if (!missingStatRep[itpt][icnd][irep]) {
                            int dataRow = (int) igen + itpt * (lastGen + 1);
                            int dataCol = irep + icnd * numReplicates;

                            diffSqr += Math.pow(cohorts.get(dataRow, dataCol) - pred - meanResidue, 2);
                        }
                    }

                    double StdDevResidue = 0.0;
                    double StdError = 0.0;
                    if (dataCounter > 1) {
                        StdDevResidue = Math.sqrt(diffSqr / (dataCounter - 1));
                        StdError = StdDevResidue / Math.sqrt(dataCounter);
                    }
                    else {
                        StdError = 0.0;
                    }

                    resGen0.add(timeStamps.get(itpt), meanResidue, StdError);
                }
            }
            DataSeries finalResGen0 = new DataSeries(condNames.get(icnd), resGen0);
            COHORT_GEN0_OUTPUT.add(finalResGen0);
        }

        return COHORT_GEN0_OUTPUT;
    }

    public double[][][] getCarParams() { return carParams; }

    public int getBOOT_ITERATIONS() { return BOOT_ITERATIONS; }
    
    public double[][][] getlmParamsGauss() { return lmParamsGauss; }
    public double[][] getlowA() { return lowA; }
    public double[][] gethighA() { return highA; }
    public double[][] getlowMu() { return lowMu; }
    public double[][] gethighMu() { return highMu; }
    public double[][] getlowSigma() { return lowSigma; }
    public double[][] gethighSigma() { return highSigma; }

    public double[][] getlmParamsPieceWise() { return lmParamsPieceWise; }
    public double[] getLowtDD() { return lowtDD; }
    public double[] getHightDD() { return hightDD; }
    public double[] getLowSlope() { return lowSlope; }
    public double[] getHighSlope() { return highSlope; }
    public double[] getLowDD() { return lowDD; }
    public double[] getHighDD() { return highDD; }
    public double[] getlowtFirstDiv() { return lowtFirstDiv; }
    public double[] getHightFirstDiv() { return hightFirstDiv; }

    public double[][] getlmCDF() { return lmCDF; }
    public double[] getLowAmpCDF() { return lowAmpCDF; }
    public double[] getHighAmpCDF() { return highAmpCDF; }
    public double[] getLowMuCDF() { return lowMuCDF; }
    public double[] getHighMuCDF() { return highMuCDF; }
    public double[] getLowSigmaCDF() { return lowSigmaCDF; }
    public double[] getHighSigmaCDF() { return highSigmaCDF; }
}
