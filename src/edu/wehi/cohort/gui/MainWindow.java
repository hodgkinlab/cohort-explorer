/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.wehi.cohort.gui;


import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.io.File;
import java.io.IOException;
import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;
import java.text.DecimalFormat;
//import java.awt.event.*;

import de.erichseifert.gral.ui.DrawablePanel;
import org.apache.commons.io.FilenameUtils;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.ejml.simple.SimpleMatrix;

import de.erichseifert.gral.ui.InteractivePanel;
import de.erichseifert.gral.graphics.Insets2D;
import de.erichseifert.gral.graphics.Dimension2D;
import de.erichseifert.gral.graphics.Location;
import de.erichseifert.gral.graphics.DrawableContainer;
import de.erichseifert.gral.graphics.layout.TableLayout;
import de.erichseifert.gral.data.DataSource;
import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.data.DataSeries;
import de.erichseifert.gral.data.statistics.Statistics;
import de.erichseifert.gral.plots.lines.LineRenderer;
import de.erichseifert.gral.plots.points.PointRenderer;
import de.erichseifert.gral.plots.axes.AxisRenderer;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.points.DefaultPointRenderer2D;

import edu.wehi.cohort.format.DataManager;
import edu.wehi.cohort.format.FitData;
import edu.wehi.cohort.io.ExportPlot;
import edu.wehi.cohort.io.ExportStat;

/**
 *
 * @author akan
 */
public class MainWindow extends JFrame {

    public static String version = "v2018.05.09";

    // TODO: make stop and pause buttons
    public static void main(String[] args) throws IOException {

        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | javax.swing.UnsupportedLookAndFeelException ex) {
                java.util.logging.Logger.getLogger(MainWindow.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }

        Filedroppo();
    }

    public MainWindow() {
//        buttonPanel = new ButtonPanel(totalPlotPanel, totalPlot);
//        getContentPane().add(buttonPanel, BorderLayout.SOUTH);
        this.setTitle("Cohort Explorer");
        this.setPreferredSize(new Dimension(900, 600));
        this.setLocationRelativeTo(null);
        this.setLocation(200, 200);
    }

    private XYPlot totalPlot;
    private InteractivePanel totalPlotPanel;

    // TODO: add more color options
//    final private Color[] COLORS = {
//            Color.BLUE,
//            Color.RED,
//            new Color(0f, 0.501960784f, 0f),
//            new Color(0.545098039f,0f,0.545098039f),
//            Color.ORANGE,
//            Color.BLACK,
//            Color.CYAN,
//            Color.BLUE.darker(),
//            Color.RED.darker(),
//            new Color(0f, 0.501960784f, 0f).darker(),
//            new Color(0.545098039f,0f,0.545098039f).darker(),
//            Color.ORANGE.darker(),
//            Color.BLACK.darker(),
//            Color.CYAN.darker()
//    };
    public static Color[] COLORS;

//    private final ButtonPanel buttonPanel;
//
//    static class ButtonPanel extends JPanel implements ActionListener {
//
//        private final InteractivePanel totalPlotPanel;
//        private final XYPlot totalPlot;
//
//        public ButtonPanel(InteractivePanel totalPlotPanel, XYPlot totalPlot) {
//                this.totalPlotPanel = totalPlotPanel;
//                this.totalPlot = totalPlot;
//
//                add(createButton("	", Color.BLUE));
//                add(createButton("Clear Drawing", null));
//        }
//
//        private JButton createButton(String text, Color background) {
//                JButton button = new JButton(text);
//                button.setBackground(background);
//                button.addActionListener(this);
//
//                return button;
//        }
//
//        @Override
//        public void actionPerformed(ActionEvent e) {
//                JButton button = (JButton) e.getSource();
//
//                if ("Clear Drawing".equals(e.getActionCommand())) {
//                        System.out.print("\"1\"");
//
//                        java.util.List<DataSource> dataList = totalPlot.getData();
//                        DataTable data = (DataTable) dataList.get(0);
//                        totalPlot.remove(data);
//
////                        DataTable currData = DataManager.loadData(5);
////                        totalPlot.add(currData);
//
////                        LineRenderer lines = new DefaultLineRenderer2D();
////                        totalPlot.setLineRenderers(currData, lines);
//
//                        //DataTable data = totalPlot.getData();
//                        //DataTable data = DataManager.loadData();
//                        //totalPlot.dataUpdated(data);
//
//                        totalPlotPanel.repaint();
//
////				java.util.List<PointRenderer> ptRenderers;
////				ptRenderers = totalPlot.getPointRenderers(data);
////				ptRenderers.get(0).setColor(color);
////
//                } else {
//                        System.out.print("\"2\"");
//
//                        java.util.List<DataSource> dataList = totalPlot.getData();
//                        DataTable data = (DataTable) dataList.get(0);
//
//                        LineRenderer lines = new DefaultLineRenderer2D();
//                        totalPlot.setLineRenderers(data, lines);
//
//                        totalPlotPanel.repaint();
//                }
//        }
//    }

    public static void Filedroppo() {

        JFrame frame = new JFrame("Cohort Explorer [" + version + "]");
        frame.setLayout(new FlowLayout());
        final JTextArea text = new JTextArea();

        JPanel progressPane = new JPanel();
        GridLayout layout = new GridLayout(0, 1, 0, 0);
        progressPane.setLayout(layout);
        JProgressBar progressBar = new JProgressBar(0, 100);
        JProgressBar progressBar2 = new JProgressBar(0, 100);
        progressBar.setStringPainted(true);
        progressBar2.setStringPainted(true);
        TitledBorder border = BorderFactory.createTitledBorder("Overall");
        TitledBorder border2 = BorderFactory.createTitledBorder("Current File");
//        border.setTitleJustification(TitledBorder.CENTER);
//        border2.setTitleJustification(TitledBorder.CENTER);
        progressBar.setBorder(border);
        progressBar2.setBorder(border2);
        progressPane.add(progressBar);
        progressPane.add(progressBar2);
        progressPane.setPreferredSize(new Dimension(600, 120));
        frame.add(progressPane);

        JScrollPane MainArea = new JScrollPane(text);
        MainArea.setPreferredSize(new Dimension(600, 648));
        frame.getContentPane().add(MainArea, BorderLayout.CENTER);
        text.append("This is Cohort Explorer " + version +"\n");
        text.append("Drag and drop 'FormatHLSEP2016' Excel files into this area...\n\n");

//        frame.setBounds(2000, 50, 600, 800);
        frame.setBounds(100, 100, 600, 800);
        frame.setResizable(false);
        frame.setDefaultCloseOperation(frame.EXIT_ON_CLOSE);
        frame.setVisible(true);

        new FileDrop(System.out, text, new FileDrop.Listener() {
            public void filesDropped(File[] files) {
//                JOptionPane.showMessageDialog(frame, "MEH");

                // TODO: create a warning popup/message when another file is dropped while it's computing
                // Put the task in a background thread -> makes the program responsive during calculation
                class Task extends SwingWorker<Void, Void>{

                    @Override
                    public Void doInBackground() throws Exception{
                        long startTime = System.nanoTime();

                        progressBar.setValue(0);
                        progressBar2.setValue(0);
                        boolean printList = true;

                        for (int i = 0; i < files.length; i++) {
                            try {
                                // Print list of dropped files
                                if(printList) {
                                    text.setText(files.length + " file(s) dropped! Analysing...\n");
                                    for (int fileIDX = 0; fileIDX < files.length; fileIDX++) {
                                        text.append("    ("+ (fileIDX+1) + ") " + files[fileIDX].getCanonicalPath() + "\n");
                                    }
                                    text.append("---------------------------------------------------------------------------------\n");
                                    printList = false;
                                }

                                // Show current working file on 'Overall' progress bar
                                Border currFile = BorderFactory.createTitledBorder("Overall \u21e2 " + files[i].getName());
                                progressBar.setBorder(currFile);

                                text.append("Working on...'" + files[i].getName() +"'");
                                File filez = new File(files[i].getCanonicalPath());
                                Dataset DATA = FormatHLSEP2016.importDataFromFile(filez);

                                // create N number of distinct colors
                                double nColors = DATA.numConditions;
                                COLORS = new Color[(int)nColors];

                                // hard-coded saturation & brightness level
                                float hue = 0.0f;
                                float saturation = 0.7f;
                                float brightness = 0.8f;

                                // user HSL color scheme for consistent color space
                                int index = 0;
                                while (index < nColors) {
                                    hue += 1. / nColors;
                                    Color c = Color.getHSBColor(hue, saturation, brightness);
                                    COLORS[index] = c;
                                    index++;
                                }

                                /* ***************************************************************************
                                 * Main plot routine begins here!
                                 * DataManager class creates an object that contains processed data.
                                 * Since GRAL requires data to be in DataTable/DataSeries class object,
                                 * 'processed data' is organised in List<DataSource>, including
                                 * fitted curves. This object also carries fitted parameters in double arrays.
                                 *****************************************************************************/
                                DataManager procDATA = new DataManager(DATA, text, progressBar2, files[i]);

                                /* ***************************************************************************
                                 * Create a new output folder for each input files
                                 * "inputFileName" Folder --inside--> figure Folder
                                 *****************************************************************************/
                                String baseName = FilenameUtils.getBaseName(files[i].getName());
                                File newFolder = new File(files[i].getParent() + "/" + baseName);
                                if (!newFolder.exists()) {
                                    newFolder.mkdir();
                                }
                                else {
                                    int fileCopyNum = 2;
                                    while (newFolder.exists()) {
                                        String Copy = String.format("_%03d", fileCopyNum++);
                                        newFolder = appendFileName(files[i], Copy);
                                        if (fileCopyNum == 100) {
                                            throw new IOException("number of file copies is too large");
                                        }
                                    }
                                    newFolder.mkdir();
                                }
                                File figFolder = new File(newFolder.getCanonicalPath() + "/" + "cohort_figs");
                                figFolder.mkdir();

                                /* ***************************************************************************
                                 * Show interactive plot windows
                                 * Export plots in figure Folder
                                 *****************************************************************************/
                                text.append("...Plotting...");
                                plotAll(files[i], figFolder, DATA, procDATA, progressBar2);
                                text.append("done\n");

                                /* ***************************************************************************
                                 * Export to excel files
                                 *****************************************************************************/
                                XSSFWorkbook wb = FormatHLSEP2016.exportDataForPrism(DATA, filez, newFolder);
                                ExportStat.toFiles(files[i], newFolder, DATA, procDATA, wb);

                                text.append("The file " + "'" + files[i].getCanonicalPath() + "'" + " is analysed. \n\n");

                                // Overall progress
                                double percent = 100*(i+1)/files.length;
                                progressBar.setValue((int)percent);

                            } catch (IOException e) {
                                text.append("\nBAD NEWS: " + e.getMessage() + "\n");
                                text.append("GOOD NEWS: you can try again ");
                                text.append("(after correcting your file or choosing another one) \n\n");
                                e.printStackTrace();

                                double percent = 100*(i+1)/files.length;
                                progressBar.setValue((int)percent);
                            }
                        }
                        Border currFile = BorderFactory.createTitledBorder("Overall");
                        progressBar.setBorder(currFile);

                        long estimatedTime = System.nanoTime() - startTime;
                        long hr = TimeUnit.NANOSECONDS.toHours(estimatedTime);
                        long min = TimeUnit.NANOSECONDS.toMinutes(estimatedTime - TimeUnit.HOURS.toNanos(hr));
                        long sec = TimeUnit.NANOSECONDS.toSeconds(estimatedTime - TimeUnit.HOURS.toNanos(hr) - TimeUnit.MINUTES.toNanos(min));
                        text.append(String.format("Total Elapsed Time - %d hrs : %d min : %d sec", hr, min, sec)+"\n");
                        text.append("\nDrag and drop another file onto the window to analyse \n \n");

                        return null;
                    }

                    // TODO: block all other incoming tasks while the program is running with pop-ups (asking user to terminate current task and start new one) - or wait for the program to finish and continue on the extra files
//                    @Override
//                    protected void done() { throw new RuntimeException("I want to produce a stack trace!"); }
                }
                Task task = new Task();
                task.execute();
            }
        });
    }

    /**
     * Creates a diagonal cross shape.
     *
     * @param armLength  the length of each 'arm'.
     * @param thickness  the thickness.
     *
     * @return A diagonal cross shape.
     */
    private static Shape createDiagonalCross(double armLength,
                                             double thickness) {
        final double SQRT2 = Math.pow(2.0, 0.5);
        final GeneralPath CROSS = new GeneralPath();

        CROSS.moveTo(-armLength - thickness, -armLength + thickness);
        CROSS.lineTo(-armLength + thickness, -armLength - thickness);
        CROSS.lineTo(0.0f, -thickness * SQRT2);
        CROSS.lineTo(armLength - thickness, -armLength - thickness);
        CROSS.lineTo(armLength + thickness, -armLength + thickness);
        CROSS.lineTo(thickness * SQRT2, 0.0f);
        CROSS.lineTo(armLength + thickness, armLength - thickness);
        CROSS.lineTo(armLength - thickness, armLength + thickness);
        CROSS.lineTo(0.0f, thickness * SQRT2);
        CROSS.lineTo(-armLength + thickness, armLength + thickness);
        CROSS.lineTo(-armLength - thickness, armLength - thickness);
        CROSS.lineTo(-thickness * SQRT2, 0.0f);
        CROSS.closePath();

        return CROSS;
    }

    private static File appendFileName(File file, String suffix) throws IOException {

        String fullFileName = file.getAbsolutePath();
        int strExtPoint = fullFileName.length() - 5;
        int strCutPoint = fullFileName.length() - 9;
        String fileExt = fullFileName.substring(strCutPoint, strExtPoint);
        fileExt = fileExt + suffix;
        fullFileName = fullFileName.substring(0, strCutPoint) + fileExt;
        file = new File(fullFileName);

        return file;
    }


    private static void plotAll(File inputFile,
                                File figFolder,
                                Dataset DATA,
                                DataManager procDATA,
                                JProgressBar progressBar) throws IOException {
        Border titledBorder = BorderFactory.createTitledBorder("Plotting...");
        progressBar.setBorder(titledBorder);
        progressBar.setValue(0);
        double totalProcess = 4;

        boolean[][] missingStat = DATA.missingStat;
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();

        // setup a main frame and tabs
        JFrame masterFRAME = new JFrame();
        masterFRAME.setSize(screenSize.width, screenSize.height);
        masterFRAME.setTitle(inputFile.getCanonicalPath());

        JTabbedPane masterTABPANE = new JTabbedPane();
        masterTABPANE.setLayout(new BorderLayout());
        masterTABPANE.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);

        // Plot [Cohort Vs. Gen] per time points
        titledBorder = BorderFactory.createTitledBorder("Plotting...[1/4 Panel]");
        progressBar.setBorder(titledBorder);

        MainWindow COHORTGEN = new MainWindow();
        COHORTGEN.PlotCohortGen(DATA, procDATA, missingStat, figFolder);
        COHORTGEN.setSize(900, 600);

        double currProcess = 100* (1.0/totalProcess);
        progressBar.setValue((int) currProcess);



        // Plot [Mean Division Number Vs. Time]
        // Plot [Arithmetically calculated Mean Division Number vs. Time]
        // Plot [Amplitude(i.e. total cohort number) Vs. Time]
        // Plot [Fitted Crash Plot]
        // Plot [Arithmetically calculated Total Cohort Number vs. Time]
        // [Arithmetic Crash plot]
        titledBorder = BorderFactory.createTitledBorder("Plotting...[2/4 Panel]");
        progressBar.setBorder(titledBorder);

        MainWindow fitMDNTIME = new MainWindow();
        XYPlot A3 = fitMDNTIME.PlotMeanDivTime(procDATA);
        MainWindow fitTOTCOHORT = new MainWindow();
        XYPlot B3 = fitTOTCOHORT.PlotAmplitude(procDATA);
        MainWindow fitCRASH = new MainWindow();
        XYPlot C3 = fitCRASH.PlotFittedCrash(procDATA);
        MainWindow arithMDNTIME = new MainWindow();
        XYPlot D3 = arithMDNTIME.PlotArithemticMeanDivTime(procDATA);
        MainWindow arithTOTCOHORT = new MainWindow();
        XYPlot E3 = arithTOTCOHORT.PlotArithemticToTCohort(procDATA);
        MainWindow arithCRASH = new MainWindow();
        XYPlot F3 = arithCRASH.PlotArithemticMDNCohort(procDATA);

        XYPlot[] FitVsArith = new XYPlot[]{A3, B3, C3, D3, E3, F3};
        ExportPlot.toSVG(FitVsArith, "3", figFolder);

        MainWindow exPLOTS = new MainWindow();
        DrawableContainer mainPlotContainer = new DrawableContainer(new TableLayout(3));
        for (XYPlot awesomePlot : FitVsArith) {
            awesomePlot.setBackground(Color.WHITE);
            mainPlotContainer.add(awesomePlot);
        }
        InteractivePanel NEW = new InteractivePanel(mainPlotContainer);
        ExportPlot.toPNG(mainPlotContainer, figFolder, "Fig3   MDN & Total Cohort & Crash");
        for (XYPlot awesomePlot : FitVsArith) {
            awesomePlot.setBackground(new Color(0.0f, 0.0f, 0.0f, 0));
        }
        exPLOTS.getContentPane().add(NEW);
        exPLOTS.pack();
        exPLOTS.setSize(900, 600);

        currProcess = 100* (2.0/totalProcess);
        progressBar.setValue((int) currProcess);



        // Plot [Sigma Vs. Time]
        // Plot [Mu Vs. Sigma]
        // Plot [gen0 cohort residue]
        titledBorder = BorderFactory.createTitledBorder("Plotting...[3/4 Panel]");
        progressBar.setBorder(titledBorder);

        MainWindow DFT = new MainWindow();
        XYPlot A4 = DFT.DistTFD(procDATA);
        MainWindow GEN0PREDICT = new MainWindow();
        XYPlot B4 = GEN0PREDICT.PlotGen0(procDATA);
        MainWindow SIGMATIME = new MainWindow();
        XYPlot C4 = SIGMATIME.PlotSigmaTime(procDATA);
        MainWindow MUSIGMA = new MainWindow();
        XYPlot D4 = MUSIGMA.PlotMuSigma(procDATA);

        XYPlot[] miscPlots = new XYPlot[]{A4, B4, C4, D4};
        ExportPlot.toSVG(miscPlots, "4", figFolder);

        MainWindow exPLOTS2 = new MainWindow();
        DrawableContainer exPlotContainer2 = new DrawableContainer(new TableLayout(2));
        for (XYPlot miscPlot : miscPlots) {
            miscPlot.setBackground(Color.WHITE);
            exPlotContainer2.add(miscPlot);
        }
        InteractivePanel NEW2 = new InteractivePanel(exPlotContainer2);
        ExportPlot.toPNG(exPlotContainer2, figFolder, "Fig4   Miscellaneous");
        for (XYPlot miscPlot : miscPlots) {
            miscPlot.setBackground(new Color(0.0f, 0.0f, 0.0f, 0));
        }
        exPLOTS2.getContentPane().add(NEW2);
        exPLOTS2.pack();
        exPLOTS2.setSize(900, 600);

        currProcess = 100* (3.0/totalProcess);
        progressBar.setValue((int) currProcess);



        // Plot [Total Cells Vs. Time]
        titledBorder = BorderFactory.createTitledBorder("Plotting...[4/4 Panel]");
        progressBar.setBorder(titledBorder);

        MainWindow TOTCELLS = new MainWindow();
        XYPlot AAA = TOTCELLS.PlotTotalCells(procDATA);
        XYPlot[] AA = new XYPlot[]{AAA};
        ExportPlot.toSVG(AA, "1", figFolder);
        AAA.setBackground(Color.WHITE);
        ExportPlot.toPNG(AA[0], figFolder, "Fig1   Total Cells vs. Time");
        AAA.setBackground(new Color(0.0f, 0.0f, 0.0f, 0));
        TOTCELLS.getContentPane().removeAll();
        TOTCELLS.PlotTotalCells(procDATA);
        TOTCELLS.setSize(900, 600);

        currProcess = 100* (4.0/totalProcess);
        progressBar.setValue((int) currProcess);

//        JScrollPane scrollArea1 = new JScrollPane(TOTCELLS.getContentPane());
        masterTABPANE.addTab("Fig1   Total Cell Number vs. Time", null, TOTCELLS.getContentPane());
//        masterTABPANE.addTab("Fig1   Total Cell Number vs. Time", scrollArea1);

//        JScrollPane scrollArea2 = new JScrollPane(COHORTGEN.getContentPane(), JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        masterTABPANE.addTab("Fig2   Cohort Number vs. Generation", null, COHORTGEN.getContentPane());
//        masterTABPANE.addTab("Fig2   Cohort Number vs. Generation", scrollArea2);

        masterTABPANE.addTab("Fig3   MDN/Total Cohort Number/Crash : Fitted & Arithmetic", null, exPLOTS.getContentPane());
        masterTABPANE.addTab("Fig4   Miscellaneous", null, exPLOTS2.getContentPane());
        masterFRAME.add(masterTABPANE, BorderLayout.CENTER);
        masterFRAME.setVisible(true);
    }

    /* **********************************************
     *  Figure 1 (Total Cells vs. Time):
     *    Basic total cell plot as a function of time
     ************************************************/
    private XYPlot PlotTotalCells(DataManager procDATA) {

        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 70.0,
                insetsBottom = 60.0,
                insetsRight = 180.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);
//        buttonPanel = new ButtonPanel(totalPlotPanel, totalPlot);

        /* Populate data & draw colored lines/datapoints per conditions */
        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.CellvTime.get(icnd));

            List<DataSource> dataList = totalPlot.getData();
            DataSeries tempData = (DataSeries) dataList.get(icnd);

            /* Set color for data points */
            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);

            /* Set error bars for data points */
            pts.setErrorColumnTop(2);
            pts.setErrorColumnBottom(2);
            pts.setErrorColor(COLORS[icnd]);
            pts.setErrorVisible(true);
            totalPlot.setPointRenderers(tempData, pts);

            LineRenderer lines = new DefaultLineRenderer2D();
            lines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(tempData, lines);
        }

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendVisible(true);
        totalPlot.getAxis(XYPlot.AXIS_Y).setMin(0.0);

        /* Styling axes */
        AxisRenderer rendererY = totalPlot.getAxisRenderer(XYPlot.AXIS_Y);
        DecimalFormat formatter = new DecimalFormat("#.#E0");
        rendererY.setTickLabelFormat(formatter);
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("time (hrs)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Total Number of Cells");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).setLabelDistance(2.0);

        totalPlot.getTitle().setText("A   Total Cells vs. Time");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    /* ***********************************************************************
     *  Figure 2 (Cohort vs. Gen):
     *    Place plots from Top to Bottom, Left to Right
     *    Number of panels would depend on the time points from the input data
     *    Each panel has following contents:
     *       1. Data points [square points + Standard Error Bars]
     *       2. Fitted discrete Gaussan curve
     *       3. Predicted Division 0 Cohort number [cross point + its value]
     *************************************************************************/
    private void PlotCohortGen(Dataset DATA, DataManager procDATA, boolean[][] missingStat, File newFolder) {

        XYPlot[] plotObjs = new XYPlot[]{};
        java.util.ArrayList<XYPlot> tempList = new ArrayList<>();
        Object[] plotLabels = new Object[]{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO", "AP", "AQ", "AR", "AS", "AT", "AU", "AV", "AW", "AX", "AY", "AZ"};

        DrawableContainer plotContainer = new DrawableContainer(new TableLayout(3));

        for (int itpt = 0; itpt < procDATA.numTimePoints; itpt++) {
            totalPlot = new XYPlot();

            /* Set size of the plot panel */
            double insetsTop = 5.0,
                    insetsLeft = 40.0,
                    insetsBottom = 60.0,
                    insetsRight = 150.0;
            totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
            totalPlot.setBackground(Color.WHITE);

            int DATACounter = 0;
            double max = 0.0;
//            double maxFit = 0.0;
            int IDX = 0;
            for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {

                if (!missingStat[itpt][icnd]) {

                    totalPlot.add(procDATA.CohortvTime.get(itpt).get(IDX));
                    totalPlot.add(procDATA.fittedGaussCurve.get(itpt).get(IDX));

                    double A = procDATA.getlmParamsGauss()[itpt][icnd][0];
                    double mu = procDATA.getlmParamsGauss()[itpt][icnd][1];
                    double sigma = procDATA.getlmParamsGauss()[itpt][icnd][2];

                    double predDiv0 = FitData.DGaussFunc(true, false, 0.0, A, mu, sigma);
                    DataTable predictedDiv0 = new DataTable(Double.class, Double.class);
                    predictedDiv0.add(0.0, predDiv0);
                    DataSeries PRED = new DataSeries("Pred", predictedDiv0);
                    totalPlot.add(PRED);
                    totalPlot.getLegend().remove(PRED);
                /* *************************************************************************************
                 * Load all the data inside the totalPlot including the fitted one
                 * So the data list is structured as following:
                 * Condition1:
                 *   {originalData = dataList[0], fittedData = dataList[1], predictedDiv0 = dataList[2]}
                 * Condition2:
                 *   {originalData = dataList[3], fittedData = dataList[4], predictedDiv0 = dataList5]}
                 * Condition3:
                 *   {originalData = dataList[6], fittedData = dataList[7], predictedDiv0 = dataList[8]}
                 * and so on
                 **************************************************************************************/
                    List<DataSource> dataList = totalPlot.getData();

                /* Split the original data & fitted data using counters */
                /* Section 1: original data */
                    DataSeries tempData = (DataSeries) dataList.get(DATACounter);
                    DATACounter += 3;

                /* Set color for data points */
                    PointRenderer DATApts = new DefaultPointRenderer2D();
                    DATApts.setShape(new Rectangle2D.Double(-3.0, -3.0, 6.0, 6.0));
                    DATApts.setColor(COLORS[icnd]);

                /* Set error bars for data points */
                    DATApts.setErrorColumnTop(2);
                    DATApts.setErrorColumnBottom(2);
                    DATApts.setErrorColor(COLORS[icnd]);
                    DATApts.setErrorVisible(true);
                    totalPlot.setPointRenderers(tempData, DATApts);

                /* Section 2: fitted one */
                    int FITCounter = IDX * 3 + 1;
                    DataSeries fittedTempData = (DataSeries) dataList.get(FITCounter);
                    LineRenderer fittedLines = new DefaultLineRenderer2D();
                    fittedLines.setColor(COLORS[icnd]);
                    totalPlot.setLineRenderers(fittedTempData, fittedLines);
                    totalPlot.setPointRenderers(fittedTempData, null);

                /* Section 3: render a point for predicted Division 0 */
                    int PREDCounter = IDX * 3 + 2;
                    DataSeries predictedDATA = (DataSeries) dataList.get(PREDCounter);
                    PointRenderer PREDpts = new DefaultPointRenderer2D();
                    PREDpts.setShape(createDiagonalCross(10.0, 0.8));
                    PREDpts.setColor(COLORS[icnd]);
                    PREDpts.setValueVisible(true);
                    PREDpts.setValueColor(COLORS[icnd]);
                    PREDpts.setValueFont(new Font("Latin Modern Math", Font.PLAIN, 13));
                    PREDpts.setValueFormat(new DecimalFormat("###"));
                    PREDpts.setValueAlignmentX(0.0);
//                PREDpts.setValueAlignmentY(0.5);
                    PREDpts.setValueDistance(1.5);
//                PREDpts.setValueRotation(45.0);
                    totalPlot.setPointRenderers(predictedDATA, PREDpts);

                /* Cut off predicted Div0 values if they are too large; focus on visualising Datapoints and fitted ones */
                    if (tempData.getColumn(1).getStatistics(Statistics.MAX) > max) {
                        max = tempData.getColumn(1).getStatistics(Statistics.MAX);
//                    totalPlot.getAxis(XYPlot.AXIS_Y).setMax(max);
                    }
//                if (fittedTempData.getColumn(1).getStatistics(Statistics.MAX) > maxFit) {
//                    maxFit = fittedTempData.getColumn(1).getStatistics(Statistics.MAX);
//                }

//                if (max > maxFit) {
//                    totalPlot.getAxis(XYPlot.AXIS_Y).setMax(max);
//                }
//                else if (maxFit > max) {
//                    totalPlot.getAxis(XYPlot.AXIS_Y).setMax(maxFit);
//                }
                    IDX++;
                }
            }
            totalPlot.getAxis(XYPlot.AXIS_Y).setMax(max);
            totalPlot.getAxis(XYPlot.AXIS_Y).setMin(0.0);
            totalPlot.getAxis(XYPlot.AXIS_X).setMin(-0.5);
            /* Plot styling */
            /* Put legend box outside of the plot area */
            totalPlot.setLegendLocation(Location.EAST);
            totalPlot.setLegendDistance(0.3);
            totalPlot.getLegend().setInsets(new Insets2D.Double(4, 4, 4, 4));
            totalPlot.getLegend().setGap(new Dimension2D.Double(2, 2));
            totalPlot.getLegend().setSymbolSize(new Dimension2D.Double(0.4,0.4));
            totalPlot.getLegend().setAlignmentX(0.0);
            totalPlot.setLegendVisible(true);

            /* Styling axes */
            AxisRenderer rendererY = totalPlot.getAxisRenderer(XYPlot.AXIS_Y);
            DecimalFormat formatter = new DecimalFormat("#");
            rendererY.setTickLabelFormat(formatter);
            totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("Generations");
            totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Cohort Number");
            totalPlot.getAxisRenderer(XYPlot.AXIS_Y).setLabelDistance(2.0);

            SimpleMatrix timeStamps = DATA.getTimeStamps();
            totalPlot.getTitle().setText(plotLabels[itpt] + "   " + timeStamps.get(itpt) +" (hrs)");
            totalPlot.getTitle().setAlignmentX(0.0);

            plotContainer.add(totalPlot);
            tempList.add(totalPlot);
        }
        ExportPlot.toSVG(tempList.toArray(plotObjs), "2", newFolder);
        totalPlotPanel = new InteractivePanel(plotContainer);

        ExportPlot.toPNG(plotContainer, newFolder, "Fig2   Cohort vs. Gen");
        for (XYPlot returnTransparency : tempList) {
            returnTransparency.setBackground(new Color(0.0f, 0.0f, 0.0f, 0));
        }

        getContentPane().add(totalPlotPanel);
        pack();
    }

    /* ************************************************
     *  Figure 3 (MDN & Total Cohort & Crash):
     *    3A - Fitted MDN vs. Time
     *    3B - Fitted Total Cohort Number vs. Time
     *    3C - Fitted Crash Plot
     *    3D - Arithmetic MDN vs. Time
     *    3E - Arithmetic Total Cohort Number vs. Time
     *    3F - Arithmetic Crash Plot
     **************************************************/
    private XYPlot PlotMeanDivTime(DataManager procDATA) {
        totalPlot = new XYPlot();
        
        /* Set size of the plot panel */
        double insetsTop = 20.0,
               insetsLeft = 60.0,
               insetsBottom = 60.0,
               insetsRight = 180.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);
        
        int oddCounter = 0;
        int evenCounter = 0;
        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.MDNvTime.get(icnd));
            totalPlot.add(procDATA.fittedPieceWise.get(icnd));
            
            List<DataSource> dataList = totalPlot.getData();

            /* Split the original data & fitted data using even/odd counters */
            DataSeries tempData = (DataSeries) dataList.get(evenCounter);
            evenCounter += 2;
            
            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);

            /* Set error bars for data points : 95% Confidence Interval from Bootstrap */
            pts.setErrorColumnTop(3);
            pts.setErrorColumnBottom(2);
            pts.setErrorColor(COLORS[icnd]);
            pts.setErrorVisible(true);
            totalPlot.setPointRenderers(tempData, pts);
          
            /* Plot lines for original datasets (per conditions) */
//            LineRenderer lines = new DefaultLineRenderer2D();
//            lines.setColor(COLORS[icnd]);
//            totalPlot.setLineRenderers(tempData, lines);
            
            oddCounter = icnd*2 + 1;
            DataSeries fittedTempData = (DataSeries) dataList.get(oddCounter);

            LineRenderer fittedLines = new DefaultLineRenderer2D();
            fittedLines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(fittedTempData, fittedLines);
            totalPlot.setPointRenderers(fittedTempData, null);
        }
   
        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendDistance(0.5);
        totalPlot.getLegend().setInsets(new Insets2D.Double(6, 6, 6, 6));
        totalPlot.getLegend().setSymbolSize(new Dimension2D.Double(0.9, 0.9));
        totalPlot.getLegend().setAlignmentX(0.0);
        totalPlot.setLegendVisible(true);
        
        /* Styling axes */
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("time (hrs)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Fitted MDN (\u00B5)");

        totalPlot.getTitle().setText("A   Fitted Mean Division Number vs. Time");
        totalPlot.getTitle().setAlignmentX(0.0);
        
        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    private XYPlot PlotAmplitude(DataManager procDATA) {
        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 110.0,
                insetsBottom = 60.0,
                insetsRight = 160.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);

        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.AmplitudevTime.get(icnd));

            List<DataSource> dataList = totalPlot.getData();
            DataSeries tempData = (DataSeries) dataList.get(icnd);

            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);
            totalPlot.setPointRenderers(tempData, pts);

            LineRenderer lines = new DefaultLineRenderer2D();
            lines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(tempData, lines);
        }
        totalPlot.getAxis(XYPlot.AXIS_Y).setMin(0);

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendDistance(0.5);
        totalPlot.setLegendVisible(true);

        // TODO: make it more general
        /* Check if maximum Y-value in the data is absurdly large */
//        if (totalPlot.getAxis(XYPlot.AXIS_Y).getMax().doubleValue() > 1E6) {
//            totalPlot.getAxis(XYPlot.AXIS_Y).setMax(1E4);
//        }

        /* Styling axes */
        AxisRenderer rendererY = totalPlot.getAxisRenderer(XYPlot.AXIS_Y);
        DecimalFormat formatter = new DecimalFormat("0");
        rendererY.setTickLabelFormat(formatter);
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("time (hrs)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Total Cohort Number (fitted A)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).setLabelDistance(2.1);

        totalPlot.getTitle().setText("B   Fitted Total Cohort Number vs. Time");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    private XYPlot PlotFittedCrash(DataManager procDATA) {
        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 45.0,
                insetsBottom = 60.0,
                insetsRight = 140.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);

        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.FittedCrash.get(icnd));

            List<DataSource> dataList = totalPlot.getData();
            DataSeries tempData = (DataSeries) dataList.get(icnd);

            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);
            totalPlot.setPointRenderers(tempData, pts);

            LineRenderer lines = new DefaultLineRenderer2D();
            lines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(tempData, lines);
        }
        totalPlot.getAxis(XYPlot.AXIS_Y).setMin(0);

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendDistance(0.5);
        totalPlot.setLegendVisible(true);

        // TODO: make it more general
        /* Make plots aesthetically visible! */
        /* Check if maximum Y-value in the data is absurdly large */
//        if (totalPlot.getAxis(XYPlot.AXIS_Y).getMax().doubleValue() > 1E6) {
//            totalPlot.getAxis(XYPlot.AXIS_Y).setMax(1E4);
//        }
        /* Check if minmum X-value in the data is negative */
        if (totalPlot.getAxis(XYPlot.AXIS_X).getMin().doubleValue() < 0) {
            totalPlot.getAxis(XYPlot.AXIS_X).setMin(0.0);
        }

        /* Styling axes */
        AxisRenderer rendererY = totalPlot.getAxisRenderer(XYPlot.AXIS_Y);
        DecimalFormat formatter = new DecimalFormat("0");
        rendererY.setTickLabelFormat(formatter);
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("Fitted MDN (\u00B5)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Total Cohort Number (fitted A)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).setLabelDistance(2.1);

        totalPlot.getTitle().setText("C   Fitted Crash Plot");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    private XYPlot PlotArithemticMeanDivTime(DataManager procDATA) {
        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 60.0,
                insetsBottom = 60.0,
                insetsRight = 180.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);

        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.ArithmeticMDNvTime.get(icnd));

            List<DataSource> dataList = totalPlot.getData();
            DataSeries tempData = (DataSeries) dataList.get(icnd);

            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);

            /* Set error bars for data points */
            pts.setErrorColumnTop(2);
            pts.setErrorColumnBottom(2);
            pts.setErrorColor(COLORS[icnd]);
            pts.setErrorVisible(true);
            totalPlot.setPointRenderers(tempData, pts);

            /* Plot lines for original datasets (per conditions) */
            LineRenderer lines = new DefaultLineRenderer2D();
            lines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(tempData, lines);
        }

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendDistance(0.5);
        totalPlot.setLegendVisible(true);

        /* Styling axes */
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("time (hrs)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("MDN");

        totalPlot.getTitle().setText("D   Arithmetic Mean Division Number vs. Time");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    private XYPlot PlotArithemticToTCohort(DataManager procDATA) {
        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 60.0,
                insetsBottom = 60.0,
                insetsRight = 160.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);

        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.ArithmeticCohortvTime.get(icnd));

            List<DataSource> dataList = totalPlot.getData();
            DataSeries tempData = (DataSeries) dataList.get(icnd);

            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);

            /* Set error bars for data points */
            pts.setErrorColumnTop(2);
            pts.setErrorColumnBottom(2);
            pts.setErrorColor(COLORS[icnd]);
            pts.setErrorVisible(true);
            totalPlot.setPointRenderers(tempData, pts);

            /* Plot lines for original datasets (per conditions) */
            LineRenderer lines = new DefaultLineRenderer2D();
            lines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(tempData, lines);
        }
        totalPlot.getAxis(XYPlot.AXIS_Y).setMin(0);

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendDistance(0.5);
        totalPlot.setLegendVisible(true);

        /* Styling axes */
        AxisRenderer rendererY = totalPlot.getAxisRenderer(XYPlot.AXIS_Y);
        DecimalFormat formatter = new DecimalFormat("0");
        rendererY.setTickLabelFormat(formatter);
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("time (hrs)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Total Cohort Number");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).setLabelDistance(2.1);

        totalPlot.getTitle().setText("E   Arithmetic Total Cohort Number vs. Time");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    private XYPlot PlotArithemticMDNCohort(DataManager procDATA) {
        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 45.0,
                insetsBottom = 60.0,
                insetsRight = 140.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);

        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.ArithmeticCrash.get(icnd));

            List<DataSource> dataList = totalPlot.getData();
            DataSeries tempData = (DataSeries) dataList.get(icnd);

            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);

            /* Set error bars for data points */
            pts.setErrorColumnTop(2);
            pts.setErrorColumnBottom(2);
            pts.setErrorColor(COLORS[icnd]);
            pts.setErrorVisible(true);
            totalPlot.setPointRenderers(tempData, pts);

            /* Plot lines for original datasets (per conditions) */
            LineRenderer lines = new DefaultLineRenderer2D();
            lines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(tempData, lines);
        }
        totalPlot.getAxis(XYPlot.AXIS_Y).setMin(0);

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendDistance(0.5);
        totalPlot.setLegendVisible(true);

        /* Styling axes */
        AxisRenderer rendererY = totalPlot.getAxisRenderer(XYPlot.AXIS_Y);
        DecimalFormat formatter = new DecimalFormat("0");
        rendererY.setTickLabelFormat(formatter);
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("MDN");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Total Cohort Number");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).setLabelDistance(2.1);

        totalPlot.getTitle().setText("F   Arithmetic Crash Plot");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    /* **********************************************************
     *  Figure 4 (Miscellaneous):
     *    4A - Distribution of time to first entry (fit with CDF)
     *    4B - Cohort [Data - Predicted] for Division 0
     *    4C - Variance vs. Time
     *    4D - Variance vs. MDN
     ************************************************************/
    private XYPlot DistTFD(DataManager procDATA) {

        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 68.0,
                insetsBottom = 60.0,
                insetsRight = 180.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);

        int dataIDX = 0;
        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.ArithmeticCohortProportion.get(icnd));
            totalPlot.add(procDATA.complementaryGaussCurve.get(icnd));
            totalPlot.add(procDATA.fittedCohortProportionCDF.get(icnd));

            List<DataSource> dataList = totalPlot.getData();

            /* Split the original data & fitted data using even/odd counters */
            DataSeries tempData = (DataSeries) dataList.get(dataIDX);
            dataIDX += 3;

            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);

            /* Set error bars for data points */
            pts.setErrorColumnTop(2);
            pts.setErrorColumnBottom(2);
            pts.setErrorColor(COLORS[icnd]);
            pts.setErrorVisible(true);
            totalPlot.setPointRenderers(tempData, pts);

            /* Plot lines for original datasets (per conditions) */
//            LineRenderer lines = new DefaultLineRenderer2D();
//            lines.setColor(COLORS[icnd]);
//            totalPlot.setLineRenderers(tempData, lines);

            int GaussFitIDX = icnd*3 + 1;
            DataSeries fittedGaussData = (DataSeries) dataList.get(GaussFitIDX);

            LineRenderer fittedGaussLines = new DefaultLineRenderer2D();
            fittedGaussLines.setColor(COLORS[icnd]);
            // dashed line for Gaussian curve
//            fittedGaussLines.setStroke(new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{9}, 0));
            // dotted line for Gaussian curve
            fittedGaussLines.setStroke(new BasicStroke(2, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[] {1, 2}, 0));
            totalPlot.setLineRenderers(fittedGaussData, fittedGaussLines);
            totalPlot.setPointRenderers(fittedGaussData, null);


            int CDFFitIDX = icnd*3 + 2;
            DataSeries fittedCDFData = (DataSeries) dataList.get(CDFFitIDX);

            LineRenderer fittedCDFLines = new DefaultLineRenderer2D();
            fittedCDFLines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(fittedCDFData, fittedCDFLines);
            totalPlot.setPointRenderers(fittedCDFData, null);

            totalPlot.getLegend().remove(fittedCDFData);
        }

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendDistance(0.5);
        totalPlot.getLegend().setInsets(new Insets2D.Double(6, 6, 6, 6));
        totalPlot.getLegend().setSymbolSize(new Dimension2D.Double(0.9, 0.9));
        totalPlot.getLegend().setAlignmentX(0.0);
        totalPlot.setLegendVisible(true);

        /* Styling axes */
        AxisRenderer rendererY = totalPlot.getAxisRenderer(XYPlot.AXIS_Y);
        DecimalFormat formatter = new DecimalFormat("0.0#");
        rendererY.setTickLabelFormat(formatter);
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("time (hrs)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Prop. of Total Cohort Number excluding Division 0");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).setLabelDistance(2.0);

        totalPlot.getTitle().setText("A   Distribution of time to first entry");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    private XYPlot PlotGen0(DataManager procDATA) {
        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 60.0,
                insetsBottom = 60.0,
                insetsRight = 180.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);

        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.gen0.get(icnd));
            List<DataSource> dataList = totalPlot.getData();

            DataSeries tempData = (DataSeries) dataList.get(icnd);

            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);
            totalPlot.setPointRenderers(tempData, pts);

            /* Set error bars for data points */
            pts.setErrorColumnTop(2);
            pts.setErrorColumnBottom(2);
            pts.setErrorColor(COLORS[icnd]);
            pts.setErrorVisible(true);
            totalPlot.setPointRenderers(tempData, pts);

            /* Plot lines for original datasets (per conditions) */
            LineRenderer lines = new DefaultLineRenderer2D();
            lines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(tempData, lines);
        }

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendVisible(true);

        // TODO: make it more general
//        if (totalPlot.getAxis(XYPlot.AXIS_Y).getMax().doubleValue() > 1E6) {
//            totalPlot.getAxis(XYPlot.AXIS_Y).setMax(1E4);
//        }
//        if (totalPlot.getAxis(XYPlot.AXIS_Y).getMin().doubleValue() < -1E6) {
//            totalPlot.getAxis(XYPlot.AXIS_Y).setMin(-5E3);
//        }

        /* Styling axes */
        AxisRenderer rendererY = totalPlot.getAxisRenderer(XYPlot.AXIS_Y);
        DecimalFormat formatter = new DecimalFormat("#.#E0");
        rendererY.setTickLabelFormat(formatter);
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("time (hrs)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Div0 Cohort Residue");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).setLabelDistance(2.0);

        totalPlot.getTitle().setText("B   Cohort [Data - Predicted] for Division 0");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    private XYPlot PlotSigmaTime(DataManager procDATA) {
        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 68.0,
                insetsBottom = 60.0,
                insetsRight = 180.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);

        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.SigmaTime.get(icnd));

            List<DataSource> dataList = totalPlot.getData();
            DataSeries tempData = (DataSeries) dataList.get(icnd);

            /* Set color for data points */
            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);

            /* Set error bars for data points */
            pts.setErrorColumnTop(2);
            pts.setErrorColumnBottom(2);
            pts.setErrorColor(COLORS[icnd]);
            pts.setErrorVisible(true);
            totalPlot.setPointRenderers(tempData, pts);

            LineRenderer lines = new DefaultLineRenderer2D();
            lines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(tempData, lines);
        }

        totalPlot.getAxis(XYPlot.AXIS_Y).setMin(0.0);

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendVisible(true);

        /* Styling axes */
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("time (hrs)");
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Variance, \u03C3\u00B2");

        totalPlot.getTitle().setText("C   Variance vs. Time");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }

    private XYPlot PlotMuSigma(DataManager procDATA) {
        totalPlot = new XYPlot();

        /* Set size of the plot panel */
        double insetsTop = 20.0,
                insetsLeft = 60.0,
                insetsBottom = 60.0,
                insetsRight = 180.0;
        totalPlot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
        totalPlotPanel = new InteractivePanel(totalPlot);

        for (int icnd = 0; icnd < procDATA.numConditions; icnd++) {
            totalPlot.add(procDATA.MuSigma.get(icnd));

            List<DataSource> dataList = totalPlot.getData();
            DataSeries tempData = (DataSeries) dataList.get(icnd);

            /* Set color for data points */
            PointRenderer pts = new DefaultPointRenderer2D();
            pts.setColor(COLORS[icnd]);

            /* Set error bars for data points */
            pts.setErrorColumnTop(2);
            pts.setErrorColumnBottom(2);
            pts.setErrorColor(COLORS[icnd]);
            pts.setErrorVisible(true);
            totalPlot.setPointRenderers(tempData, pts);

            LineRenderer lines = new DefaultLineRenderer2D();
            lines.setColor(COLORS[icnd]);
            totalPlot.setLineRenderers(tempData, lines);
        }
        totalPlot.getAxis(XYPlot.AXIS_Y).setMin(0.0);

        /* Plot styling */
        /* Put legend box outside of the plot area */
        totalPlot.setLegendLocation(Location.EAST);
        totalPlot.setLegendVisible(true);

        /* Styling axes */
        totalPlot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("Variance, \u03C3\u00B2");
        totalPlot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("MDN, \u00B5");

        totalPlot.getTitle().setText("D   Variance vs. MDN");
        totalPlot.getTitle().setAlignmentX(0.0);

        getContentPane().add(totalPlotPanel);
        pack();

        return totalPlot;
    }
}
