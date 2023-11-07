package edu.wehi.cohort.io;


import java.awt.*;
import java.text.DecimalFormat;
import java.util.Calendar;
import java.util.List;
import java.io.IOException;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.lang.Double;

import org.ejml.simple.SimpleMatrix;
import org.apache.commons.io.FilenameUtils;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.util.CellUtil;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFFont;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.util.CellRangeAddress;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import edu.wehi.cohort.gui.MainWindow;
import edu.wehi.cohort.format.DataManager;
import edu.wehi.cohort.gui.Dataset;

public class ExportStat {

    public static void toFiles(File inputFile, File dest, Dataset dataset, DataManager outputFittedData, XSSFWorkbook workbook) throws IOException {

        String baseNameinputFile = FilenameUtils.getBaseName(inputFile.getName());

        /* ***********************
         * PRINT META INFORMATION
         ************************/
        Calendar calendar = Calendar.getInstance();

        File metafile = new File(dest.getCanonicalPath() + "/" + "info.txt");
        Writer out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(metafile), StandardCharsets.UTF_8));

        // Print program version
        out.write("CohortExplorer " + MainWindow.version + "\r\n");

        // Print time stamp
        out.write("Created at : " + calendar.getTime().toString() + "\r\n");

        // Print input file info
        out.write("Input file : " + inputFile.getCanonicalPath() + "\r\n");

        // Print bootstrap iterations
        out.write("Bootstrap iterations per condition : " + outputFittedData.getBOOT_ITERATIONS() + "\r\n");

        Object[] plotLabels = new Object[]{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO", "AP", "AQ", "AR", "AS", "AT", "AU", "AV", "AW", "AX", "AY", "AZ"};
        // Print output files & tree
        out.write("Output files : \r\n\r\n");
        out.write("   " + dest.getName() + "\r\n");
        out.write("      ├── " + baseNameinputFile + "_4prism.xlsx" + "\r\n");
        out.write("      └── " + "cohort_figs" + "\r\n");
        out.write("           ├── " + "Fig1   Total Cells vs. Time.png" + "\r\n");
        out.write("           ├── " + "Fig2   Cohort vs. Gen.png" + "\r\n");
        out.write("           ├── " + "Fig3   MDN & Total Cohort & Crash.png" + "\r\n");
        out.write("           ├── " + "Fig4   Miscellaneous.png" + "\r\n");
        out.write("           ├── " + "1A   Total Cells vs. Time.svg" + "\r\n");
        for (int justLoop = 0; justLoop < dataset.numTimePoints; justLoop++) {
            out.write("           ├── 2" + plotLabels[justLoop] + "   " + dataset.getTimeStamps().get(justLoop) +" (hrs).svg \r\n");
        }
        out.write("           ├── " + "3A   Fitted Mean Division Number vs. Time.svg" + "\r\n");
        out.write("           ├── " + "3B   Fitted Total Cohort Number vs. Time.svg" + "\r\n");
        out.write("           ├── " + "3C   Fitted Crash Plot.svg" + "\r\n");
        out.write("           ├── " + "3D   Arithmetic Mean Division Number vs. Time.svg" + "\r\n");
        out.write("           ├── " + "3E   Arithmetic Total Cohort Number vs. Time.svg" + "\r\n");
        out.write("           ├── " + "3F   Arithmetic Crash Plot.svg" + "\r\n");
        out.write("           ├── " + "4A   Distribution of time to first entry.svg" + "\r\n");
        out.write("           ├── " + "4B   Cohort [Data - Predicted] for Division 0.svg" + "\r\n");
        out.write("           ├── " + "4C   Variance vs. Time.svg" + "\r\n");
        out.write("           └── " + "4D   Variance vs. MDN.svg" + "\r\n");

        out.flush();
        out.close();

        /* ***********************
         * PRINT SATAISTICS
         ************************/
//        XSSFWorkbook workbook = new XSSFWorkbook();

        exportCohortGenInfo(dataset, outputFittedData, workbook);
        exportMDNFitInfo(dataset, outputFittedData, workbook);
        exportCDFFitInfo(dataset, outputFittedData, workbook);

//        File file = new File(dest.getCanonicalPath() + "/" + baseNameinputFile + "_fitResults.xlsx");
//        FileOutputStream fileOutputStream = new FileOutputStream(file);
        File newFile1 = appendFileName(inputFile, "_4prism");
        String finalDestination = dest.getCanonicalPath() + "/" + newFile1.getName();
        File finalFile = new File(finalDestination);
        FileOutputStream outputStream = new FileOutputStream(finalFile);
        workbook.write(outputStream);
    }

    static File appendFileName(File file, String suffix) throws IOException {

        String fullFileName = file.getAbsolutePath();
        int strCutPoint = fullFileName.length() - 5;
        String fileExt = fullFileName.substring(strCutPoint);
        if (!fileExt.equals(".xlsx")) {
            throw new IOException();
        }

        fileExt = suffix + fileExt;
        fullFileName = fullFileName.substring(0, strCutPoint) + fileExt;
        file = new File(fullFileName);

        return file;
    }

    private static void exportCohortGenInfo(Dataset dataset, DataManager outputFittedData, XSSFWorkbook workbook) {
        XSSFSheet sheet = workbook.createSheet("fig2 cohort gen");

        XSSFCellStyle style = workbook.createCellStyle();
        XSSFFont font = workbook.createFont();
        font.setColor(IndexedColors.RED.getIndex());
        style.setFont(font);

        XSSFCellStyle style2 = workbook.createCellStyle();
        XSSFFont font2 = workbook.createFont();
        font2.setColor(IndexedColors.BLUE.getIndex());
        style2.setFont(font2);

        XSSFCellStyle styleBold = workbook.createCellStyle();
        XSSFFont fontBold = workbook.createFont();
        fontBold.setBold(true);
        styleBold.setFont(fontBold);

        DecimalFormat df = new DecimalFormat( "#0.0000");

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix timeStamps = dataset.getTimeStamps();
        double[][][] fittedVals = outputFittedData.getlmParamsGauss();
        double[][] lowMu = outputFittedData.getlowMu();
        double[][] highMu = outputFittedData.gethighMu();

        Row row1 = sheet.createRow(0);
        Cell cell1 = row1.createCell(0);
        cell1.setCellValue("fitted MDN (\u00B5)");
        cell1.setCellStyle(styleBold);

        row1 = sheet.createRow(1);
        cell1 = row1.createCell(0);
        cell1.setCellValue("time");
        cell1.setCellStyle(styleBold);


        double[][] lowSigma = outputFittedData.getlowSigma();
        double[][] highSigma = outputFittedData.gethighSigma();
        int rowOffset = dataset.numTimePoints + 3;

        Row row2 = sheet.createRow(rowOffset);
        Cell cell2 = row2.createCell(0);
        cell2.setCellValue("fitted standard deviation (\u03C3)");
        cell2.setCellStyle(styleBold);

        row2 = sheet.createRow(1 + rowOffset);
        cell2 = row2.createCell(0);
        cell2.setCellValue("time");
        cell2.setCellStyle(styleBold);


        double[][] lowA = outputFittedData.getlowA();
        double[][] highA = outputFittedData.gethighA();
        int rowOffset2 = 2*(dataset.numTimePoints + 3);

        Row row3 = sheet.createRow(rowOffset2);
        Cell cell3 = row3.createCell(0);
        cell3.setCellValue("fitted Amplitude (A)");
        cell3.setCellStyle(styleBold);

        row3 = sheet.createRow(1 + rowOffset2);
        cell3 = row3.createCell(0);
        cell3.setCellValue("time");
        cell3.setCellStyle(styleBold);

        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            cell1 = row1.createCell(3*icnd+1);
            cell1.setCellValue(currName);
            cell1.setCellStyle(styleBold);

            cell1 = row1.createCell(3*icnd+2);
            cell1.setCellValue("95% Conf Interval");
            cell1.setCellStyle(styleBold);
            CellUtil.setAlignment(cell1, workbook, CellStyle.ALIGN_CENTER);
            sheet.addMergedRegion(new CellRangeAddress(1, 1, 3*icnd+2, 3*icnd+3));

            cell2 = row2.createCell(3*icnd+1);
            cell2.setCellValue(currName);
            cell2.setCellStyle(styleBold);
            cell2 = row2.createCell(3*icnd+2);
            cell2.setCellValue("95% Conf Interval");
            cell2.setCellStyle(styleBold);
            CellUtil.setAlignment(cell2, workbook, CellStyle.ALIGN_CENTER);
            sheet.addMergedRegion(new CellRangeAddress(1 + rowOffset, 1 + rowOffset, 3*icnd+2, 3*icnd+3));

            cell3 = row3.createCell(3*icnd+1);
            cell3.setCellValue(currName);
            cell3.setCellStyle(styleBold);
            cell3 = row3.createCell(3*icnd+2);
            cell3.setCellValue("95% Conf Interval");
            cell3.setCellStyle(styleBold);
            CellUtil.setAlignment(cell3, workbook, CellStyle.ALIGN_CENTER);
            sheet.addMergedRegion(new CellRangeAddress(1 + rowOffset2, 1 + rowOffset2, 3*icnd+2, 3*icnd+3));
        }

        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row1 = sheet.createRow(itpt + 2);
            cell1 = row1.createCell(0);
            cell1.setCellValue(currTime);

            row2 = sheet.createRow(itpt + 2 + rowOffset);
            cell2 = row2.createCell(0);
            cell2.setCellValue(currTime);

            row3 = sheet.createRow(itpt + 2 + rowOffset2);
            cell3 = row3.createCell(0);
            cell3.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                cell1 = row1.createCell(3*icnd+1);
                cell1.setCellValue(Double.valueOf(df.format(fittedVals[itpt][icnd][1])));
                cell1 = row1.createCell(3*icnd+2);
                cell1.setCellValue(Double.valueOf(df.format(lowMu[itpt][icnd])));
                cell1.setCellStyle(style);
                cell1 = row1.createCell(3*icnd+3);
                cell1.setCellValue(Double.valueOf(df.format(highMu[itpt][icnd])));
                cell1.setCellStyle(style2);


                cell2 = row2.createCell(3*icnd+1);
                if (Double.isNaN(fittedVals[itpt][icnd][2])) {
                    cell2.setCellValue("NaN");
                }
                else {
                    cell2.setCellValue(Double.valueOf(df.format(fittedVals[itpt][icnd][2])));
                }

                if (Double.isNaN(lowSigma[itpt][icnd]) || Double.isNaN(highSigma[itpt][icnd])) {
                    cell2 = row2.createCell(3*icnd+2);
                    cell2.setCellValue("NaN");
                    cell2.setCellStyle(style);
                    cell2 = row2.createCell(3*icnd+3);
                    cell2.setCellValue("NaN");
                    cell2.setCellStyle(style2);
                }
                else {
                    cell2 = row2.createCell(3*icnd+2);
                    cell2.setCellValue(Double.valueOf(df.format(lowSigma[itpt][icnd])));
                    cell2.setCellStyle(style);
                    cell2 = row2.createCell(3*icnd+3);
                    cell2.setCellValue(Double.valueOf(df.format(highSigma[itpt][icnd])));
                    cell2.setCellStyle(style2);
                }

                cell3 = row3.createCell(3*icnd+1);
                cell3.setCellValue(Double.valueOf(df.format(fittedVals[itpt][icnd][0])));
                cell3 = row3.createCell(3*icnd+2);
                cell3.setCellValue(Double.valueOf(df.format(lowA[itpt][icnd])));
                cell3.setCellStyle(style);
                cell3 = row3.createCell(3*icnd+3);
                cell3.setCellValue(Double.valueOf(df.format(highA[itpt][icnd])));
                cell3.setCellStyle(style2);

            }
        }
    }

    private static void exportMDNFitInfo(Dataset dataset, DataManager outputFittedData, XSSFWorkbook workbook) {
        XSSFSheet sheet = workbook.createSheet("fig3A fitted MDN vs. Time");

        XSSFCellStyle style = workbook.createCellStyle();
        XSSFFont font = workbook.createFont();
        font.setColor(IndexedColors.RED.getIndex());
        style.setFont(font);

        XSSFCellStyle style2 = workbook.createCellStyle();
        XSSFFont font2 = workbook.createFont();
        font2.setColor(IndexedColors.BLUE.getIndex());
        style2.setFont(font2);

        XSSFCellStyle styleBold = workbook.createCellStyle();
        XSSFFont fontBold = workbook.createFont();
        fontBold.setBold(true);
        styleBold.setFont(fontBold);

        DecimalFormat df = new DecimalFormat( "#0.0000");
        List<String> condNames = dataset.getConditionNames();

        Row row1 = sheet.createRow(0);
        Cell cell1 = row1.createCell(0);
        cell1.setCellValue("[MDN vs. Time] fit results");
        cell1.setCellStyle(styleBold);

        row1 = sheet.createRow(1);
        cell1 = row1.createCell(0);
        cell1.setCellValue("condition");
        cell1.setCellStyle(styleBold);

        cell1 = row1.createCell(1);
        cell1.setCellValue("estimated time to enter first division (tFD)");
        sheet.autoSizeColumn(1);
        cell1.setCellStyle(styleBold);

        cell1 = row1.createCell(2);
        cell1.setCellValue("95% Conf Interval");
        cell1.setCellStyle(styleBold);
        CellUtil.setAlignment(cell1, workbook, CellStyle.ALIGN_CENTER);
        sheet.addMergedRegion(new CellRangeAddress(1, 1, 2, 3));

        cell1 = row1.createCell(4);
        cell1.setCellValue("estimated division destiny (DD)");
        sheet.autoSizeColumn(4);
        cell1.setCellStyle(styleBold);

        cell1 = row1.createCell(5);
        cell1.setCellValue("95% Conf Interval");
        cell1.setCellStyle(styleBold);
        CellUtil.setAlignment(cell1, workbook, CellStyle.ALIGN_CENTER);
        sheet.addMergedRegion(new CellRangeAddress(1, 1, 5, 6));

        cell1 = row1.createCell(7);
        cell1.setCellValue("estimated time to division destiny (tDD)");
        sheet.autoSizeColumn(7);
        cell1.setCellStyle(styleBold);

        cell1 = row1.createCell(8);
        cell1.setCellValue("95% Conf Interval");
        cell1.setCellStyle(styleBold);
        CellUtil.setAlignment(cell1, workbook, CellStyle.ALIGN_CENTER);
        sheet.addMergedRegion(new CellRangeAddress(1, 1, 8, 9));

        cell1 = row1.createCell(10);
        cell1.setCellValue("estimated time to enter from k to k+1 division (m)");
        sheet.autoSizeColumn(10);
        cell1.setCellStyle(styleBold);

        cell1 = row1.createCell(11);
        cell1.setCellValue("95% Conf Interval");
        cell1.setCellStyle(styleBold);
        CellUtil.setAlignment(cell1, workbook, CellStyle.ALIGN_CENTER);
        sheet.addMergedRegion(new CellRangeAddress(1, 1, 11, 12));

        double[][] fittedVals = outputFittedData.getlmParamsPieceWise();

        double[] lowDD = outputFittedData.getLowDD();
        double[] highDD = outputFittedData.getHighDD();

        double[] lowtDD = outputFittedData.getLowtDD();
        double[] hightDD = outputFittedData.getHightDD();

        double[] lowSlope = outputFittedData.getLowSlope();
        double[] highSlope = outputFittedData.getHighSlope();

        double[] lowtFirstDiv = outputFittedData.getlowtFirstDiv();
        double[] hightFirstDiv = outputFittedData.getHightFirstDiv();

        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {

            row1 = sheet.createRow(icnd + 2);
            cell1 = row1.createCell(0);
            cell1.setCellValue(condNames.get(icnd));

            cell1 = row1.createCell(1);
            cell1.setCellValue(Double.valueOf(df.format(fittedVals[icnd][3])));
            cell1 = row1.createCell(2);
            cell1.setCellValue(Double.valueOf(df.format(lowtFirstDiv[icnd])));
            cell1.setCellStyle(style);
            cell1 = row1.createCell(3);
            cell1.setCellValue(Double.valueOf(df.format(hightFirstDiv[icnd])));
            cell1.setCellStyle(style2);

            cell1 = row1.createCell(4);
            cell1.setCellValue(Double.valueOf(df.format(fittedVals[icnd][1])));
            cell1 = row1.createCell(5);
            cell1.setCellValue(Double.valueOf(df.format(lowDD[icnd])));
            cell1.setCellStyle(style);
            cell1 = row1.createCell(6);
            cell1.setCellValue(Double.valueOf(df.format(highDD[icnd])));
            cell1.setCellStyle(style2);

            cell1 = row1.createCell(7);
            cell1.setCellValue(Double.valueOf(df.format(fittedVals[icnd][0])));
            cell1 = row1.createCell(8);
            cell1.setCellValue(Double.valueOf(df.format(lowtDD[icnd])));
            cell1.setCellStyle(style);
            cell1 = row1.createCell(9);
            cell1.setCellValue(Double.valueOf(df.format(hightDD[icnd])));
            cell1.setCellStyle(style2);

            cell1 = row1.createCell(10);
            cell1.setCellValue(Double.valueOf(df.format(fittedVals[icnd][2])));
            cell1 = row1.createCell(11);
            cell1.setCellValue(Double.valueOf(df.format(lowSlope[icnd])));
            cell1.setCellStyle(style);
            cell1 = row1.createCell(12);
            cell1.setCellValue(Double.valueOf(df.format(highSlope[icnd])));
            cell1.setCellStyle(style2);
        }
    }

    private static void exportCDFFitInfo(Dataset dataset, DataManager outputFittedData, XSSFWorkbook workbook) {
        XSSFSheet sheet = workbook.createSheet("fig4A distribution of time to first entry");

        XSSFCellStyle style = workbook.createCellStyle();
        XSSFFont font = workbook.createFont();
        font.setColor(IndexedColors.RED.getIndex());
        style.setFont(font);

        XSSFCellStyle style2 = workbook.createCellStyle();
        XSSFFont font2 = workbook.createFont();
        font2.setColor(IndexedColors.BLUE.getIndex());
        style2.setFont(font2);

        XSSFCellStyle styleBold = workbook.createCellStyle();
        XSSFFont fontBold = workbook.createFont();
        fontBold.setBold(true);
        styleBold.setFont(fontBold);

        DecimalFormat df = new DecimalFormat( "#0.0000");
        List<String> condNames = dataset.getConditionNames();

        Row row1 = sheet.createRow(0);
        Cell cell1 = row1.createCell(0);
        cell1.setCellValue("[Distribution of time to first entry] CDF fit results");
        cell1.setCellStyle(styleBold);

        row1 = sheet.createRow(1);
        cell1 = row1.createCell(0);
        cell1.setCellValue("condition");
        cell1.setCellStyle(styleBold);

        cell1 = row1.createCell(1);
        cell1.setCellValue("fitted mean time to enter first division (\u00B5)");
        sheet.autoSizeColumn(1);
        cell1.setCellStyle(styleBold);

        cell1 = row1.createCell(2);
        cell1.setCellValue("95% Conf Interval");
        cell1.setCellStyle(styleBold);
        CellUtil.setAlignment(cell1, workbook, CellStyle.ALIGN_CENTER);
        sheet.addMergedRegion(new CellRangeAddress(1, 1, 2, 3));

        cell1 = row1.createCell(4);
        cell1.setCellValue("fitted standard deviation (\u03C3)");
        sheet.autoSizeColumn(4);
        cell1.setCellStyle(styleBold);

        cell1 = row1.createCell(5);
        cell1.setCellValue("95% Conf Interval");
        cell1.setCellStyle(styleBold);
        CellUtil.setAlignment(cell1, workbook, CellStyle.ALIGN_CENTER);
        sheet.addMergedRegion(new CellRangeAddress(1, 1, 5, 6));

        cell1 = row1.createCell(7);
        cell1.setCellValue("fitted Amplitude (A)");
        sheet.autoSizeColumn(7);
        cell1.setCellStyle(styleBold);

        cell1 = row1.createCell(8);
        cell1.setCellValue("95% Conf Interval");
        cell1.setCellStyle(styleBold);
        CellUtil.setAlignment(cell1, workbook, CellStyle.ALIGN_CENTER);
        sheet.addMergedRegion(new CellRangeAddress(1, 1, 8, 9));

        double[][] fittedVals = outputFittedData.getlmCDF();

        double[] lowMuCDF = outputFittedData.getLowMuCDF();
        double[] highMuCDF = outputFittedData.getHighMuCDF();

        double[] lowSigmaCDF = outputFittedData.getLowSigmaCDF();
        double[] highSigmaCDF = outputFittedData.getHighSigmaCDF();

        double[] lowAmpCDF = outputFittedData.getLowAmpCDF();
        double[] highAmpCDF = outputFittedData.getHighAmpCDF();

        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            row1 = sheet.createRow(icnd + 2);
            cell1 = row1.createCell(0);
            cell1.setCellValue(condNames.get(icnd));

            cell1 = row1.createCell(1);
            cell1.setCellValue(Double.valueOf(df.format(fittedVals[icnd][1])));
            cell1 = row1.createCell(2);
            cell1.setCellValue(Double.valueOf(df.format(lowMuCDF[icnd])));
            cell1.setCellStyle(style);
            cell1 = row1.createCell(3);
            cell1.setCellValue(Double.valueOf(df.format(highMuCDF[icnd])));
            cell1.setCellStyle(style2);

            cell1 = row1.createCell(4);
            cell1.setCellValue(Double.valueOf(df.format(fittedVals[icnd][2])));
            cell1 = row1.createCell(5);
            cell1.setCellValue(Double.valueOf(df.format(lowSigmaCDF[icnd])));
            cell1.setCellStyle(style);
            cell1 = row1.createCell(6);
            cell1.setCellValue(Double.valueOf(df.format(highSigmaCDF[icnd])));
            cell1.setCellStyle(style2);

            cell1 = row1.createCell(7);
            cell1.setCellValue(Double.valueOf(df.format(fittedVals[icnd][0])));
            cell1 = row1.createCell(8);
            cell1.setCellValue(Double.valueOf(df.format(lowAmpCDF[icnd])));
            cell1.setCellStyle(style);
            cell1 = row1.createCell(9);
            cell1.setCellValue(Double.valueOf(df.format(highAmpCDF[icnd])));
            cell1.setCellStyle(style2);
        }
    }
}
