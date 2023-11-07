/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.wehi.cohort.format;


import de.erichseifert.gral.data.DataTable;
import org.apache.commons.math3.special.Erf;
import org.ejml.simple.SimpleMatrix;
import org.ejml.data.SingularMatrixException;

import edu.wehi.cohort.gui.Dataset;

/**
 *
 * @author HoChan
 */

public class FitData {

    public static Object[] CaruanaFit (double[] xData, double[] yData) throws SingularMatrixException {
        /* **************************************************************************
         * Initiate fitting: 
         *    Solve m*X=n for X vector
         *    Since our model function becomes quadratic equation upon linearlisation,
         *    algorithm only requires 3x3 matrix
         ****************************************************************************/
        double[][] m = new double[3][3];
        double[] n = new double[3];

        for (int i = 0; i < yData.length; i++) {
            /* ***************************************************************
             * Check for zero values and replace it with small positive number
             * This is necessary as the algorithm takes log of the data point
             *****************************************************************/
            if (yData[i] == 0) {
                yData[i] = 1e-50;
            }

            // Guo's algorithm: weighted (with yData^2) Caruana
            m[0][0] += Math.pow(yData[i], 2);
            m[0][1] += xData[i] * Math.pow(yData[i], 2);
            m[1][0] += xData[i] * Math.pow(yData[i], 2);
            m[0][2] += Math.pow(xData[i], 2) * Math.pow(yData[i], 2);
            m[1][1] += Math.pow(xData[i], 2) * Math.pow(yData[i], 2);
            m[2][0] += Math.pow(xData[i], 2) * Math.pow(yData[i], 2);
            m[1][2] += Math.pow(xData[i], 3) * Math.pow(yData[i], 2);
            m[2][1] += Math.pow(xData[i], 3) * Math.pow(yData[i], 2);
            m[2][2] += Math.pow(xData[i], 4) * Math.pow(yData[i], 2);

            n[0] += Math.pow(yData[i], 2) * Math.log(yData[i]);
            n[1] += xData[i]*Math.pow(yData[i], 2) * Math.log(yData[i]);
            n[2] += Math.pow(xData[i], 2) * Math.pow(yData[i], 2) * Math.log(yData[i]);
        }

        // Copy m,n matrix to M,N SimepleMatrix objects
        SimpleMatrix M = new SimpleMatrix(m);
        SimpleMatrix N = new SimpleMatrix(3, 1, true, n);

        // Solve M*x=N for vector x; x = M^(-1) * N
        try {
            SimpleMatrix RESULTS = M.solve(N);
            
            double a = RESULTS.get(0);
            double b = RESULTS.get(1);
            double c = RESULTS.get(2);

            // Transform the values back to original form
            double A = Math.exp(a - Math.pow(b,2)/(4*c));
            double mu = -b/(2*c);
            double sigma = Math.sqrt(-1/(2*c));

            return new Object[]{A, mu, sigma};
        }
        catch (SingularMatrixException e) {

            double A = Double.NaN;
            double mu = Double.NaN;
            double sigma = Double.NaN;

            return new Object[]{A, mu, sigma};
        }
    }



    public static double DGaussFunc(boolean predict, boolean lastGen,
                                    double x, double A, double mu, double sigma) {

        double NormFactor = A*Math.sqrt(2.0*Math.PI*Math.pow(sigma,2));

        double z1 = (x+1-mu)/(sigma*Math.sqrt(2));
        double z2 = (x-mu)/(sigma*Math.sqrt(2));

        if (lastGen) {
            int genBeyond = 30;
            double total = 0;
            for (int i=1; i < genBeyond; i++ ) {
                z1 = (x+i-mu)/(sigma*Math.sqrt(2));
                z2 = (x+i-1-mu)/(sigma*Math.sqrt(2));
                
                total += Math.pow(2, x+i-1) * 0.5 * ( Erf.erf(z1) - Erf.erf(z2) );
            }

            total = NormFactor * total / Math.pow(2, x);
            return total;
        }
        else if (predict && x==0.0) {
            return NormFactor * (0.5 + 0.5 * Erf.erf(z1));
        }
        else {
            return NormFactor * 0.5 * ( Erf.erf(z1) - Erf.erf(z2) );
        }
    }

    // Numerical differentiation
    private static double Diff(boolean lastGen, int param,
                               double x, double A, double mu, double sigma) {
        double slope = 0;
        double eps = 1E-5;
        
        // Quadrature differentiation : 4 points interpolation
        if (param == 0) {
            slope = ( DGaussFunc(false, lastGen, x, A-2*eps, mu, sigma)
                    - 8*DGaussFunc(false, lastGen, x, A-eps, mu, sigma)
                    + 8*DGaussFunc(false, lastGen, x, A+eps, mu, sigma)
                    - DGaussFunc(false, lastGen, x, A+2*eps, mu, sigma) ) / (12*eps);
        }
        else if (param == 1) {
            slope = ( DGaussFunc(false, lastGen, x, A, mu-2*eps, sigma)
                    - 8*DGaussFunc(false, lastGen, x, A, mu-eps, sigma)
                    + 8*DGaussFunc(false, lastGen, x, A, mu+eps, sigma)
                    - DGaussFunc(false, lastGen, x, A, mu+2*eps, sigma) ) / (12*eps);
        }
        else if (param == 2) {
            slope = ( DGaussFunc(false, lastGen, x, A, mu, sigma-2*eps)
                    - 8*DGaussFunc(false, lastGen, x, A, mu, sigma-eps)
                    + 8*DGaussFunc(false, lastGen, x, A, mu, sigma+eps)
                    - DGaussFunc(false, lastGen, x, A, mu, sigma+2*eps) ) / (12*eps);
        }
        return slope;
    }
    
    public static DataTable lmGaussFit(Dataset DATA, int itpt,
                                       double[] xData, double[] yData,
                                       int Ndata, int countZeros, int xCordMax,
                                       double max, double min,
                                       double A, double mu, double sigma) {
        // Small/Big number flags
        final double SMALLVAL = 1E-4;
        final double BIGVAL = 1E4;
        
        // Initial guesses from Caruana
        double A_est = A;
        double mu_est = mu;
        double sigma_est = sigma;

        // Check Caruana fitted parameters for NaN, Big, or Small values
        if (Double.isNaN(A) || Double.isInfinite(A) || A > BIGVAL || A < min || itpt < 3) {
                // Maximum Y value in the DATA
            A_est = max;
        }
        if (Double.isNaN(mu) || Double.isInfinite(mu) || mu < 0 || mu > DATA.lastGeneration || itpt < 3) {
                // X coordinate at maximum Y
            mu_est = xData[xCordMax];
        }
        if (Double.isNaN(sigma) || Double.isInfinite(sigma) || sigma < SMALLVAL || sigma > DATA.lastGeneration) {
            sigma_est = 1.0;
        }

        // NO DATA!
        if (countZeros == Ndata || (double)countZeros/Ndata > 0.9) {
            mu_est = 0.0;
            sigma_est = Double.NaN;
        }
        else if ((double)countZeros/Ndata > 0.5) {
            A_est = max;
            mu_est = xData[xCordMax];
            sigma_est = 1.0;
        }

        // LM configuration
        int Niter = 3000;
        int Nparams = 3;
        double eps = 1E-4;
        boolean updateJacobian = true;
        boolean lastGen = false;
        
        // LM variables
        double lambda = 1E-4;
        double error = 0.0;
        
        double[] yData_est = new double[Ndata];
        double[] yData_est_lm = new double[Ndata];
        
        SimpleMatrix Jacobian = new SimpleMatrix(Ndata, Nparams);
        SimpleMatrix dError = new SimpleMatrix(Ndata, 1);
        SimpleMatrix dError_lm = new SimpleMatrix(Ndata, 1);
        SimpleMatrix Hessian = new SimpleMatrix(Nparams, Nparams);

        // Initiate Levenberg-Marquardt algorithm
        for (int iter = 0; iter < Niter; iter++) {
            
            // updateJacobian
            if (updateJacobian) {
                
                Jacobian = new SimpleMatrix(Ndata, Nparams);
                
                // Populate Jacobian matrix
                for (int row = 0; row < Ndata; row++) {
                    for (int col = 0; col < Nparams; col++) {
                        if (row == Ndata-1) lastGen = true;
                        Jacobian.set(row, col, -Diff(lastGen, col, xData[row], A_est, mu_est, sigma_est));
                    }

                    // estimate y with estimated parameters
                    yData_est[row] = DGaussFunc(false, lastGen, xData[row], A_est, mu_est, sigma_est);
                    lastGen = false;

                    // calculate difference y_estimate and y_data
                    dError.set(row, 0, yData[row] - yData_est[row]);
                }

                // compute approximate Hessian
                Hessian = Jacobian.transpose().mult(Jacobian);

                // compute initial error
                if (iter == 0) {
                    error = dError.dot(dError);
                }
            } // End of updateJacobian

            SimpleMatrix Hessian_lm = SimpleMatrix.identity(Nparams).scale(lambda).plus(Hessian);

            /* *****************************************
             *  Stopping criteria 1:
             *    Breaks the loop if Hessian is singular
             *******************************************/
            if (Hessian_lm.conditionP2() > 1E15 || Hessian_lm.determinant() < 1E-4) {
                break;
            }

            SimpleMatrix params = Hessian_lm.invert().negative().mult(Jacobian.transpose().mult(dError));

            // update parameters
            double A_lm = A_est + params.get(0);
            double mu_lm = mu_est + params.get(1);
            double sigma_lm = sigma_est + params.get(2);

            /* ******************************************************************************
             *  Stopping criteria 2
             *    Tolerence check: Norm vector of difference of previous & current parameters
             ********************************************************************************/
            double diffParams = Math.pow(A_lm-A_est, 2) + Math.pow(mu_lm-mu_est, 2) + Math.pow(sigma_lm-sigma_est, 2);
            double diffNorm = Math.sqrt(diffParams);
            double Norm = Math.sqrt(Math.pow(A_est, 2) + Math.pow(mu_est, 2) + Math.pow(sigma_est, 2));
            if (diffNorm <= eps*(Norm+eps)) {
                break;
            }

            for (int row = 0; row < Ndata; row++) {
                if (row == Ndata-1) lastGen = true;
                yData_est_lm[row] = DGaussFunc(false, lastGen, xData[row], A_lm, mu_lm, sigma_lm);
                lastGen = false;
                dError_lm.set(row, 0, yData[row] - yData_est_lm[row]);
            }

            double error_lm = dError_lm.dot(dError_lm);

            // check if lambda needs to be decreased or increased
            if (error_lm < error) {
                updateJacobian = true;
                lambda /= 2.0;
                
                A_est = A_lm;
                mu_est = mu_lm;
                sigma_est = sigma_lm;
                
                error = error_lm;
            }
            else {
                updateJacobian = false;
                lambda *= 2.0;
            }
        } // End of Levenberg-Marquardt algorithm

        DataTable fittedGaussianCurve
                = new DataTable(Integer.class, Double.class, Double.class, Double.class, Double.class);

        for (int i = 0; i <= DATA.lastGeneration; i++) {
            double yVals = 0.0;
            if (i == 0) {
                yVals = DGaussFunc(true, false, (double)i, A_est, mu_est, sigma_est);
            }
            else{
                yVals = DGaussFunc(false, false, (double)i, A_est, mu_est, sigma_est);
            }

            if (yVals < 1E-50 && yVals >= 0) {
                yVals = 0.0;
            }

            fittedGaussianCurve.add(i, yVals, A_est, mu_est, sigma_est);
        }
        return fittedGaussianCurve;
    }



    private static double piecewiseLinearFunc(double x,
                                              double bp, double a, double b,
                                              double D1, double D2) {
        /* *********************************************************
         *  This model has 3 parameters to be estimated:
         *    1) a : y-intercept
         *    2) b : slope of linear function at x < bp (breakpoint)
         *    3) bp : x value at which two functions connect
         *
         *    Note that D1, D2 are NOT parameters
         *    {D1, D2} = {1, 0} for x <= bp
         *    {D1, D2} = {0, 1} for x > bp
         *
         *  So following quantities can be calculated:
         *     Division Destiny = a + b*bp
         *     slope (overall Division Rate) = b
         *     Division Rate per unit MDN = 1/b
         **********************************************************/
        return a + b*x*D1 + b*bp*D2;
    }

    private static double dfda () { return 1.0; }

    private static double dfdb (double x, double bp, double D1, double D2) { return x*D1 + bp*D2; }

    private static double dfdbp (double b, double D2) { return b*D2; }

    public static DataTable lmPiecewise(double[] xData, double[] yData, double[] WEIGHT, int Ndata) {

        Object[] LinReg =  SimpleLinearRegression(xData, yData, Ndata);

        // Initial guesses
        double totTime = 0.0;
        for (int i = 0; i < Ndata; i++) {
            totTime += xData[i];
        }

        double bp_est = totTime/Ndata;
        double b_est = (double)LinReg[0];
        double a_est = (double)LinReg[1];

        // LM configuration
        int Niter = 3000;
        int Nparams = 3;
        boolean updateJacobian = true;

        // LM variables
        double D1 = 0.0;
        double D2 = 0.0;
        double lambda = 1E-6;
        double error = 0;
        double eps = 1E-6;
        double[] yData_est = new double[Ndata];
        double[] yData_est_lm = new double[Ndata];
        SimpleMatrix Jacobian = new SimpleMatrix(Ndata, Nparams);
        SimpleMatrix dError = new SimpleMatrix(Ndata, 1);
        SimpleMatrix dError_lm = new SimpleMatrix(Ndata, 1);
        SimpleMatrix Hessian = new SimpleMatrix(Nparams, Nparams);

        SimpleMatrix diagWEIGHT = new SimpleMatrix(Ndata, Ndata);
        for (int row = 0; row < Ndata; row++) {
            for (int col = 0; col < Ndata; col++) {
                if (row == col) {
                    diagWEIGHT.set(row, col, WEIGHT[row]);
                }
            }
        }

        // Initiate Levenberg-Marquardt algorithm
        for (int iter = 0; iter < Niter; iter++) {

            // updateJacobian
            if (updateJacobian) {
                Jacobian = new SimpleMatrix(Ndata, Nparams);

                // Populate Jacobian matrix
                for (int row = 0; row < Ndata; row++) {
                    if (xData[row] <= bp_est) {
                        D1 = 1.0;
                        D2 = 0.0;
                    }
                    else {
                        D1 = 0.0;
                        D2 = 1.0;
                    }
                    for (int col = 0; col < Nparams; col++) {
                        if (col==0) {
                            Jacobian.set(row, col, -dfda());
                        }
                        else if (col == 1) {
                            Jacobian.set(row, col, -dfdb(xData[row], bp_est, D1, D2));
                        }
                        else {
                            Jacobian.set(row, col, -dfdbp(b_est, D2));
                        }
                    }
                    // estimate y with estimated parameters
                    yData_est[row] = piecewiseLinearFunc(xData[row], bp_est, a_est, b_est, D1, D2);

                    // calculate difference y_estimate and y_data
                    dError.set(row, 0, yData[row] - yData_est[row]);
                }

                // compute approximate Hessian
//                Hessian = Jacobian.transpose().mult(Jacobian);
                Hessian = Jacobian.transpose().mult(diagWEIGHT).mult(Jacobian);

                // compute initial error
                if (iter == 0) {
                    error = dError.dot(dError);
                }
            } // End of updateJacobian

            // Levenberg update rule : (lambda * Identity) + (Jt * J)
//            SimpleMatrix Hessian_lm = SimpleMatrix.identity(Nparams).scale(lambda).plus(Hessian);

            // Marquardt update rule : Construct diag(Jt * W * J) matrix for replacing Identity matrix in Levenberg's update rule
            double[] temp = new double[Nparams];
            temp[0] = Hessian.diag().get(0);
            temp[1] = Hessian.diag().get(1);
            temp[2] = Hessian.diag().get(2);
            SimpleMatrix Hessian_lm = SimpleMatrix.diag(temp).scale(lambda).plus(Hessian);

            /* *****************************************
             *  Stopping criteria 1:
             *    Breaks the loop if Hessian is singular
             *******************************************/
            if (Hessian_lm.conditionP2() > 1E15 || Hessian_lm.determinant() <= 1E-4) {
                break;
            }

//            SimpleMatrix params = Hessian_lm.invert().scale(-1).mult(Jacobian.transpose().mult(dError));
            SimpleMatrix params = Hessian_lm.invert().negative().mult(Jacobian.transpose().mult(diagWEIGHT).mult(dError));

            // update parameters
            double a_lm = a_est + params.get(0);
            double b_lm = b_est + params.get(1);
            double bp_lm = bp_est + params.get(2);

            /* ******************************************************************************
             *  Stopping criteria 2
             *    Tolerence check: Norm vector of difference of previous & current parameters
             ********************************************************************************/
            double diffParams = Math.pow(bp_lm-bp_est, 2) + Math.pow(a_lm-a_est, 2) + Math.pow(b_lm-b_est, 2);
            double diffNorm = Math.sqrt(diffParams);
            double Norm = Math.sqrt(Math.pow(bp_est, 2) + Math.pow(a_est, 2) + Math.pow(b_est, 2));
            if (diffNorm <= eps*(Norm+eps)) {
                break;
            }

            for (int row = 0; row < Ndata; row++) {
                if (xData[row] <= bp_lm) {
                    D1 = 1.0;
                    D2 = 0.0;
                }
                else {
                    D1 = 0.0;
                    D2 = 1.0;
                }

                yData_est_lm[row] = piecewiseLinearFunc(xData[row], bp_lm, a_lm, b_lm, D1, D2);
                dError_lm.set(row, 0, yData[row] - yData_est_lm[row]);
            }
            double error_lm = dError_lm.dot(dError_lm);

            // check if lambda needs to be decreased or increased
            if (error_lm < error) {
                lambda /= 10.0;

                bp_est = bp_lm;
                a_est = a_lm;
                b_est = b_lm;

                error = error_lm;
                updateJacobian = true;
            }
            else {
                updateJacobian = false;
                lambda *= 10.0;
            }
        } // End of Levenberg-Marquardt algorithm

        DataTable fittedPieceWise = new DataTable(Double.class, Double.class, Double.class, Double.class, Double.class);

        double deltaT = (xData[Ndata-1] - xData[0])/1000.0;
        for (double x = xData[0]; x <= xData[Ndata-1]; x += deltaT) {
            if (x <= bp_est) {
                D1 = 1.0;
                D2 = 0.0;
            }
            else {
                D1 = 0.0;
                D2 = 1.0;
            }
            double y = piecewiseLinearFunc(x, bp_est, a_est, b_est, D1, D2);
            fittedPieceWise.add(x, y, bp_est, a_est, b_est);
            if (deltaT == 0.0) {
                break;
            }
        }

        return fittedPieceWise;
    }

    private static Object[] SimpleLinearRegression(double[] x, double[] y, int Ndata) {

        double sumx = 0.0, sumy = 0.0;
        for (int n = 0; n < Ndata; n++) {
            sumx  += x[n];
            sumy  += y[n];
        }
        double xbar = sumx / Ndata;
        double ybar = sumy / Ndata;

        double xxbar = 0.0, xybar = 0.0;
        for (int i = 0; i < Ndata; i++) {
            xxbar += (x[i] - xbar) * (x[i] - xbar);
            xybar += (x[i] - xbar) * (y[i] - ybar);
        }
        double slope = xybar / xxbar;
        double intercept = ybar - slope * xbar;

        return new Object[]{slope, intercept};
    }



    private static double CDF(double x, double Amp, double mu, double sigma) {
        double z = (x-mu)/(sigma*Math.sqrt(2.0));
        return Amp * 0.5 * (1.0 + Erf.erf(z));
    }

    // TODO: check amplitude?
    private static double GaussFunc(double x, double Amp, double mu, double sigma) {
//        return 1.0/(Math.sqrt(2*Math.PI*Math.pow(sigma, 2))) * Math.exp(-Math.pow(x - mu, 2)/(2.0*Math.pow(sigma, 2)));
        return Amp*Math.exp(-Math.pow(x - mu, 2)/(2.0*Math.pow(sigma, 2)));
    }

    private static double DiffCDF(int param, double x, double Amp, double mu, double sigma) {
        double slope = 0;
        double eps = 1E-5;

        // Quadrature differentiation : 4 points interpolation
        if (param == 0) {
            slope = ( CDF(x, Amp-2*eps, mu, sigma)
                    - 8*CDF(x, Amp-eps, mu, sigma)
                    + 8*CDF(x, Amp+eps, mu, sigma)
                    - CDF(x, Amp+2*eps, mu, sigma) ) / (12*eps);
        }
        else if (param == 1) {
            slope = ( CDF(x, Amp, mu-2*eps, sigma)
                    - 8*CDF(x, Amp, mu-eps, sigma)
                    + 8*CDF(x, Amp, mu+eps, sigma)
                    - CDF(x, Amp, mu+2*eps, sigma) ) / (12*eps);
        }
        else if (param == 2) {
            slope = ( CDF(x, Amp, mu, sigma-2*eps)
                    - 8*CDF(x, Amp, mu, sigma-eps)
                    + 8*CDF(x, Amp, mu, sigma+eps)
                    - CDF(x, Amp, mu, sigma+2*eps) ) / (12*eps);
        }
        return slope;
    }

    public static Object[] lmCDF(double[] xData, double[] yData,
                                 int Ndata, int countZeros, int xCordMax,
                                 double max, double min) {

        // find an approximate inflection point
        double maxSlope = 0;
        int xCordatMaxSlope = 0;
        for (int i = 1; i < Ndata; i++) {
            if (xData[i] - xData[i-1] != 0.0) {
                double slope = (yData[i] - yData[i-1]) / (xData[i] - xData[i-1]);
                if (slope > maxSlope) {
                    maxSlope = slope;
                    xCordatMaxSlope = i-1;
                }
            }
        }

        // Initial guesses
        double Amp_est = max;
        double mu_est = xData[xCordatMaxSlope];
        double sigma_est = 0.21*(xData[Ndata-1] - xData[0]);

        // NO DATA!
//        if (countZeros == Ndata || (double)countZeros/Ndata > 0.9) {
//            mu_est = 0.0;
//            sigma_est = Double.NaN;
//        }
//        else if ((double)countZeros/Ndata > 0.5) {
//            A_est = max;
//            mu_est = xData[xCordMax];
//            sigma_est = 1.0;
//        }

        // LM configuration
        int Niter = 3000;
        int Nparams = 3;
        double eps = 1E-4;
        boolean updateJacobian = true;

        // LM variables
        double lambda = 1E-4;
        double error = 0.0;

        double[] yData_est = new double[Ndata];
        double[] yData_est_lm = new double[Ndata];

        SimpleMatrix Jacobian = new SimpleMatrix(Ndata, Nparams);
        SimpleMatrix dError = new SimpleMatrix(Ndata, 1);
        SimpleMatrix dError_lm = new SimpleMatrix(Ndata, 1);
        SimpleMatrix Hessian = new SimpleMatrix(Nparams, Nparams);

        // Initiate Levenberg-Marquardt algorithm
        for (int iter = 0; iter < Niter; iter++) {

            // updateJacobian
            if (updateJacobian) {

                Jacobian = new SimpleMatrix(Ndata, Nparams);

                // Populate Jacobian matrix
                for (int row = 0; row < Ndata; row++) {
                    for (int col = 0; col < Nparams; col++) {
                        Jacobian.set(row, col, -DiffCDF(col, xData[row], Amp_est, mu_est, sigma_est));
                    }

                    // estimate y with estimated parameters
                    yData_est[row] = CDF(xData[row], Amp_est, mu_est, sigma_est);

                    // calculate difference y_estimate and y_data
                    dError.set(row, 0, yData[row] - yData_est[row]);
                }

                // compute approximate Hessian
                Hessian = Jacobian.transpose().mult(Jacobian);

                // compute initial error
                if (iter == 0) {
                    error = dError.dot(dError);
                }
            } // End of updateJacobian

            SimpleMatrix Hessian_lm = SimpleMatrix.identity(Nparams).scale(lambda).plus(Hessian);

            /* *****************************************
             *  Stopping criteria 1:
             *    Breaks the loop if Hessian is singular
             *******************************************/
            if (Hessian_lm.conditionP2() > 1E15 || Hessian_lm.determinant() < 1E-6) {
                break;
            }

            SimpleMatrix params = Hessian_lm.invert().negative().mult(Jacobian.transpose().mult(dError));

            // update parameters
            double Amp_lm = Amp_est + params.get(0);
            double mu_lm = mu_est + params.get(1);
            double sigma_lm = sigma_est + params.get(2);

            /* ******************************************************************************
             *  Stopping criteria 2
             *    Tolerence check: Norm vector of difference of previous & current parameters
             ********************************************************************************/
            double diffParams = Math.pow(Amp_lm-Amp_est, 2) + Math.pow(mu_lm-mu_est, 2) + Math.pow(sigma_lm-sigma_est, 2);
            double diffNorm = Math.sqrt(diffParams);
            double Norm = Math.sqrt(Math.pow(Amp_est, 2) + Math.pow(mu_est, 2) + Math.pow(sigma_est, 2));
            if (diffNorm <= eps*(Norm+eps)) {
                break;
            }

            for (int row = 0; row < Ndata; row++) {
                yData_est_lm[row] = CDF(xData[row], Amp_lm, mu_lm, sigma_lm);
                dError_lm.set(row, 0, yData[row] - yData_est_lm[row]);
            }

            double error_lm = dError_lm.dot(dError_lm);

            // check if lambda needs to be decreased or increased
            if (error_lm < error) {
                updateJacobian = true;
                lambda /= 10.0;

                Amp_est = Amp_lm;
                mu_est = mu_lm;
                sigma_est = sigma_lm;

                error = error_lm;
            }
            else {
                updateJacobian = false;
                lambda *= 10.0;
            }
        } // End of Levenberg-Marquardt algorithm

        DataTable fittedGaussianCurve
                = new DataTable(Double.class, Double.class, Double.class, Double.class, Double.class);
        DataTable fittedCDFCurve = new DataTable(Double.class, Double.class);

        double dT = (xData[Ndata-1] - xData[0])/1000.0;
        for (double t = xData[0]; t <= xData[Ndata-1]; t += dT) {

            double yVals = GaussFunc(t, Amp_est, mu_est, sigma_est);

            if (yVals < 1E-50 && yVals >= 0) {
                yVals = 0.0;
            }

            double yVals2 = CDF(t, Amp_est, mu_est, sigma_est);

            fittedGaussianCurve.add(t, yVals, Amp_est, mu_est, sigma_est);
            fittedCDFCurve.add(t, yVals2);
        }

        return new Object[]{fittedGaussianCurve, fittedCDFCurve};
    }



    /* TODO: discuss further about log-normal fitting for Div0 prediction; this could be useful for estimating progress factor. */
//    public static double lognormFunc(double x, double A, double mu, double sigma) {
//        if (x == 0) {
//            return 0.0;
//        }
//       return (A/x) * Math.exp(-0.5 * (Math.pow(Math.log(x)-mu, 2))/Math.pow(sigma, 2));
//    }
//
//
//    private static double Diff3(int param, double x, double A, double mu, double sigma) {
//        double slope = 0;
//        double eps = 1E-5;
//
//        if (param == 0) {
//            slope = ( lognormFunc(x, A-2*eps, mu, sigma) - 8*lognormFunc(x, A-eps, mu, sigma)
//                    + 8*lognormFunc(x, A+eps, mu, sigma) - lognormFunc(x, A+2*eps, mu, sigma)  )/(12*eps);
//        }
//        else if (param == 1) {
//            slope = ( lognormFunc(x, A, mu-2*eps, sigma) - 8*lognormFunc(x, A, mu-eps, sigma)
//                    + 8*lognormFunc(x, A, mu+eps, sigma) - lognormFunc(x, A, mu+2*eps, sigma)  )/(12*eps);
//        }
//        else {
//            slope = ( lognormFunc(x, A, mu, sigma-2*eps) - 8*lognormFunc(x, A, mu, sigma-eps)
//                    + 8*lognormFunc(x, A, mu, sigma+eps) - lognormFunc(x, A, mu, sigma+2*eps)  )/(12*eps);
//        }
//
//        return slope;
//    }
//
//    public static DataTable lmLognormFit(double[] time, double[] gen0Res) {
//        int Ndata = time.length;
//
//        double xCordMax = 0;
//        double max = 0;
//        for (int i = 0; i < Ndata; i++) {
//
//            if (gen0Res[i] > max) {
//                max = gen0Res[i];
//                xCordMax = time[i];
//            }
//        }
//
//        /* Initial guesses */
//        double A_est = max;
//        double mu_est = xCordMax;
//        double sigma_est = 30.0;
//
//        /* LM configuration */
//        int Niter = 2000;
//        int Nparams = 3;
//        int updateJacobian = 1;
//
//        /* LM variables */
//        double lambda = 1E-3;
//        double error = 0;
//
//        double[] yData_est = new double[Ndata];
//        double[] yData_est_lm = new double[Ndata];
//
//        SimpleMatrix Jacobian = new SimpleMatrix(Ndata, Nparams);
//        SimpleMatrix dError = new SimpleMatrix(Ndata, 1);
//        SimpleMatrix dError_lm = new SimpleMatrix(Ndata, 1);
//        SimpleMatrix Hessian = new SimpleMatrix(Nparams, Nparams);
//
//        /* Initiate Levenberg-Marquardt algorithm */
//        for (int iter = 0; iter < Niter; iter++) {
//
//            /* updateJacobian*/
//            if (updateJacobian == 1) {
//
//                Jacobian = new SimpleMatrix(Ndata, Nparams);
//
//                /* Populate Jacobian matrix */
//                for (int row = 0; row < Ndata; row++) {
//                    for (int col = 0; col < Nparams; col++) {
//                        Jacobian.set(row, col, -Diff3(col, time[row], A_est, mu_est, sigma_est));
//                    }
//
//                    /* estimate y with estimated parameters */
//                    yData_est[row] = lognormFunc(time[row], A_est, mu_est, sigma_est);
//
//                    /* calculate difference y_estimate and y_data */
//                    dError.set(row, 0, gen0Res[row] - yData_est[row]);
//                }
//
//                /* compute Hessian */
//                Hessian = Jacobian.transpose().mult(Jacobian);
//
//                /* compute initial error */
//                if (iter == 0) {
//                    error = dError.dot(dError);
//                }
//            } // End of updateJacobian
//
//            SimpleMatrix Hessian_lm = SimpleMatrix.identity(Nparams).scale(lambda).plus(Hessian);
//
//            /* Breaks the loop if Hessian is singular */
//            if (Hessian_lm.determinant() == 0) {
//                break;
//            }
//
//            SimpleMatrix params = Hessian_lm.invert().scale(-1).mult(Jacobian.transpose().mult(dError));
//
//            /* update parameters */
//            double A_lm = A_est + params.get(0);
//            double mu_lm = mu_est + params.get(1);
//            double sigma_lm = sigma_est + params.get(2);
//
//            for (int row = 0; row < Ndata; row++) {
//                yData_est_lm[row] = lognormFunc(time[row], A_lm, mu_lm, sigma_lm);
//                dError_lm.set(row, 0, gen0Res[row] - yData_est_lm[row]);
//            }
//            double error_lm = dError_lm.dot(dError_lm);
//
//            /* check if lambda needs to be decreased or increased */
//            if (error_lm < error) {
//                lambda /= 1.2;
//
//                A_est = A_lm;
//                mu_est = mu_lm;
//                sigma_est = sigma_lm;
//
//                error = error_lm;
//                updateJacobian = 1;
//            }
//            else {
//                updateJacobian = 0;
//                lambda *= 1.2;
//            }
//        } // End of Levenberg-Marquardt algorithm
//
//
//        DataTable fittedLogNormCurve
//                = new DataTable(Double.class, Double.class, Double.class, Double.class, Double.class);
//
//        double deltaT = (time[Ndata-1] - time[0])/2000.0;
//        for (double x = time[0]; x <= time[Ndata-1]; x += deltaT) {
//            double y = lognormFunc(x, A_est, mu_est, sigma_est);
//            fittedLogNormCurve.add(x, y, A_est, mu_est, sigma_est);
//        }
//
//        return fittedLogNormCurve;
//    }
}
