package edu.wehi.cohort.gui;


import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.List;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.NoSuchElementException;

import org.apache.poi.ss.usermodel.*;
import org.ejml.simple.SimpleMatrix;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.apache.poi.xssf.usermodel.XSSFSheet;

/**
 * handles data from an Excel file in HL-SEP-2016 format
 *
 * @author akan
 */
public final class FormatHLSEP2016 {

    public static Dataset importDataFromFile(File file) throws IOException {
        class ExperimentSetup {

            int numTimePoints;
            int lastGeneration;
            int numConditions;
            int numReplicates;
            double numBeadsPut;
            double propBeadsGated;
            double expNumBeads;

            ExperimentSetup(XSSFSheet sheet) throws IOException {
                numTimePoints = 0;
                lastGeneration = 0;
                numConditions = 0;
                numReplicates = 0;
                numBeadsPut = 0;
                propBeadsGated = 0;
                expNumBeads = 0;

                Row row = sheet.getRow(2);
                Cell cell = row.getCell(1);
                if (cell == null) {
                    throw new IOException("number of time points not set");
                }
                else if (cell.getCellType() == Cell.CELL_TYPE_NUMERIC) {
                    double val = cell.getNumericCellValue();
                    if (val == Math.floor(val)) {
                        numTimePoints = (int) val;
                    }
                }

                row = sheet.getRow(3);
                cell = row.getCell(1);
                if (cell == null) {
                    throw new IOException("last generation not set");
                }
                else if (cell.getCellType() == Cell.CELL_TYPE_NUMERIC) {
                    double val = cell.getNumericCellValue();
                    if (val == Math.floor(val)) {
                        lastGeneration = (int) val;
                    }
                }

                row = sheet.getRow(4);
                cell = row.getCell(1);
                if (cell == null) {
                    throw new IOException("number of conditions not set");
                }
                else if (cell.getCellType() == Cell.CELL_TYPE_NUMERIC) {
                    double val = cell.getNumericCellValue();
                    if (val == Math.floor(val)) {
                        numConditions = (int) val;
                    }
                }

                row = sheet.getRow(5);
                cell = row.getCell(1);
                if (cell == null) {
                    throw new IOException("number of replicates not set");
                }
                else if (cell.getCellType() == Cell.CELL_TYPE_NUMERIC) {
                    double val = cell.getNumericCellValue();
                    if (val == Math.floor(val)) {
                        numReplicates = (int) val;
                    }
                }

                row = sheet.getRow(6);
                cell = row.getCell(1);
                if (cell == null) {
                    throw new IOException("number of beads not set");
                }
                else if (cell.getCellType() == Cell.CELL_TYPE_NUMERIC) {
                    double val = cell.getNumericCellValue();
                    if (val == Math.floor(val)) {
                        numBeadsPut = val;
                    }
                }

                row = sheet.getRow(7);
                cell = row.getCell(1);
                if (cell == null) {
                    throw new IOException("prop. beads gated not set");
                }
                else if (cell.getCellType() == Cell.CELL_TYPE_NUMERIC) {
                    double val = cell.getNumericCellValue();
                    if (val >= 0 && val <= 1) {
                        propBeadsGated = val;
                    }
                    else {
                        throw new IOException("prop. beads gated must be in range 0 and 1");
                    }
                }
                expNumBeads = propBeadsGated * numBeadsPut;
            }

            void testValidity() throws IOException {
                if (numTimePoints <= 0) {
                    throw new IOException("number of time points cannot be 0");
                }
                else if (numTimePoints == 1) {
                    throw new IOException("you must have more than 1 time point for piece-wise fitting" + "\n                 (In fact, you need at least 2 data points for every conditions)");
                }
                if (lastGeneration <= 0) {
                    throw new IOException("I noticed last generation is set 0! but why...?");
                }
                if (numConditions <= 0) {
                    throw new IOException("number of conditions is set 0, so there're no conditions to plot??? I'm confused.");
                }
//                else if(numConditions >= 15) {
//                    throw new IOException("Currently this program can only handle up to 14 conditions!" + "\n                 (This limit will be removed in next update, I promise)");
//                }
                if (numReplicates <= 0) {
                    throw new IOException("number of replicates is set 0, so is that mean there's no data?");
                }
                if (expNumBeads <= 0) {
                    throw new IOException("expected number of beads is 0. check your number of beads & prop. beads gated");
                }
            }
        }

        String baseErrorMessage = "problem with data importing";
        String debugMessage = baseErrorMessage;

        try (FileInputStream fis = new FileInputStream(file)) {
            XSSFWorkbook workbook = new XSSFWorkbook(fis);
            XSSFSheet sheet = workbook.getSheetAt(0);

            verifyTableFormat(sheet);

            ExperimentSetup exp = new ExperimentSetup(sheet);
            exp.testValidity();

            int numRows = (exp.lastGeneration + 1) * exp.numTimePoints;
            int numCols = exp.numConditions * exp.numReplicates;

            double[][][] totalCounts
                    = new double[exp.numTimePoints][exp.numConditions][exp.numReplicates];

            SimpleMatrix counts = new SimpleMatrix(numRows, numCols);
            SimpleMatrix cohorts = new SimpleMatrix(numRows, numCols);
            SimpleMatrix missing = new SimpleMatrix(numRows, numCols);
            missing = missing.plus(1); // initial mark all as missing

            SimpleMatrix timeStamps = new SimpleMatrix(exp.numTimePoints, 1);
            List<String> conditionNames = new ArrayList<>(exp.numConditions);

            // collect condition names
            for (int i = 0; i < exp.numConditions; i++) {
                Row row = sheet.getRow(11 + i);
                Cell cell = row.getCell(0, Row.CREATE_NULL_AS_BLANK);
                if (cell.getCellType() == Cell.CELL_TYPE_NUMERIC) {
                    String condName = String.valueOf(cell.getNumericCellValue());
                    conditionNames.add(condName);
                } else if (cell.getCellType() == Cell.CELL_TYPE_STRING) {
                    String condName = cell.getStringCellValue();
                    conditionNames.add(condName);
                } else if (cell.getCellType() == Cell.CELL_TYPE_BLANK) {
                    throw new IOException("there's an empty condition name. please specify one");
                }
            }

            // collect time stamps
            for (int i = 0; i < exp.numTimePoints; i++) {
                Row row = sheet.getRow(i * exp.numConditions * exp.numReplicates + 1);
                if (row == null) {
                    throw new IOException("problem with gathering time information at/after [time point: " + timeStamps.get(i-1) + "hrs]" + "\n                 perhaps you added extra time point which doesn't exist"+  "\n                 OR check your number of conditions & replicates" + "\n                 each time point must have equal size of data (put an empty row if no data)" );
                }
                Cell cell = row.getCell(3);
                if (cell == null || cell.getCellType() == Cell.CELL_TYPE_BLANK) {
                    throw new IOException("problem with gathering time information at/after [time point: " + timeStamps.get(i-1) + "hrs]" + "\n                 check your number of conditions & replicates" + "\n                 each time point must have equal size of data (put an empty row if no data)");
                }
                else {
                    CellType cellType = cell.getCellTypeEnum();

                    // check for cell type to avoid error in getting value.
                    double timeStamp;
                    if (cellType == CellType.NUMERIC) {
                        timeStamp = cell.getNumericCellValue();
                    } else if (cellType == CellType.STRING) {
                        timeStamp = Double.valueOf(cell.getStringCellValue());
                    }
                    else {
                        throw new IOException("problem with gathering time information at/after [time point: " + timeStamps.get(i-1) + "hrs]" + "\n                 make sure your cell type to be either numeric or string in the excel file");
                    }
                    timeStamps.set(i, timeStamp); // return -1 if any error
                }
            }

            // collect the data, compute cohort numbers
            int currRow = 1;
            for (int itpt = 0; itpt < exp.numTimePoints; itpt++) {
                for (int icnd = 0; icnd < exp.numConditions; icnd++) {
                    for (int irep = 0; irep < exp.numReplicates; irep++) {
                        debugMessage = String.format(
                                "%s: (%d,%d,%d) {%f,%s}",
                                baseErrorMessage, itpt, icnd, irep,
                                timeStamps.get(itpt), conditionNames.get(icnd));

                        Row row = sheet.getRow(currRow++);
                        if (row == null) {
                            continue;
                        }

                        Cell cell = row.getCell(5);
                        if (cell == null) {
                            continue;
                        }
                        if (cell.getCellType() == Cell.CELL_TYPE_NUMERIC) {
                            double recoveredBeads = cell.getNumericCellValue();

                            cell = row.getCell(6);
                            double v;
                            v = cell.getNumericCellValue();
                            v = v * exp.expNumBeads / recoveredBeads;
                            totalCounts[itpt][icnd][irep] = v;

                            for (int igen = 0; igen <= exp.lastGeneration; igen++) {
                                cell = row.getCell(7 + igen);

                                if (cell == null) {
                                    v = 0;
                                } else {
                                    v = cell.getNumericCellValue();
                                    v = v * exp.expNumBeads / recoveredBeads;
                                }

                                int dataRow = igen + itpt * (exp.lastGeneration + 1);
                                int dataCol = irep + icnd * (exp.numReplicates);

                                counts.set(dataRow, dataCol, v);
                                cohorts.set(dataRow, dataCol, v / Math.pow(2, igen));
                                missing.set(dataRow, dataCol, 0); // not missing
                            }
                        }
                    }
                }
            }

            return new Dataset(
                    file.getName(), totalCounts, counts, cohorts,
                    missing, timeStamps, conditionNames);

        } catch (FileNotFoundException | NoSuchElementException | NumberFormatException | IllegalStateException ex) {

            System.err.println(ex);
            throw new IOException(debugMessage);

        } catch (IOException e) {
            throw new IOException(e.getMessage());
        } catch (Exception ex) {
            System.err.println(ex);
            throw new IOException("unidentified problem with exporting data");
        }
    }

    public static XSSFWorkbook exportDataForPrism(Dataset dataset, File file, File destination) throws IOException {

        try {
            File newFile1 = appendFileName(file, "_4prism");
//            File newFile2 = appendFileName(file, "_ArithmeticResults");

//            int fileCopyNum = 2;
//            while (newFile.exists()) {
//                String suffix = String.format("_copy_%02d", fileCopyNum++);
//                newFile = appendFileName(file, suffix);
//                if (fileCopyNum == 100) {
//                    throw new IOException("number of file copies is too large");
//                }
//            }

            XSSFWorkbook workbook = new XSSFWorkbook();
            exportTotalCounts(dataset, workbook);
            exportCohortPerTime(dataset, workbook);
            exportCohortPerCondition(dataset, workbook);
            exportMeanDivNum(dataset, workbook);
            exportMeanDivNumreps(dataset, workbook);
            exportNormCohortNum(dataset, workbook);
            exportNormCohortNumreps(dataset, workbook);
            exporttotalcohortNum(dataset, workbook);
            exportTotalCohortNumreps(dataset, workbook);
            exportmdnvscohort(dataset, workbook);
            // export for CytonSolver - legacy data format
            exportTotalCellPerCond(dataset, workbook);

            String finalDestination = destination.getCanonicalPath() + "/" + newFile1.getName();
            File finalFile = new File(finalDestination);
            FileOutputStream outputStream = new FileOutputStream(finalFile);
            workbook.write(outputStream);

//            XSSFWorkbook workbook2 = new XSSFWorkbook();
//            exportMeanDivNum(dataset, workbook2);
//            exportMeanDivNumreps(dataset, workbook2);
//            exportNormCohortNum(dataset, workbook2);
//            exportNormCohortNumreps(dataset, workbook2);
//            exporttotalcohortNum(dataset, workbook2);
//            exportTotalCohortNumreps(dataset, workbook2);
//            exportmdnvscohort(dataset, workbook2);
//            exportTimeStamp(dataset, workbook);

//            String finalDestination2 = destination.getCanonicalPath() + "/" + newFile2.getName();
//            File finalFile2 = new File(finalDestination2);
//            FileOutputStream outputStream2 = new FileOutputStream(finalFile2);
//            workbook2.write(outputStream2);

            return workbook;

        } catch (IOException ex) {
            System.err.println(ex);
            throw ex;

        } catch (IndexOutOfBoundsException ex) {
            System.err.println(ex);
            throw new IOException("a data indexing problem");

        } catch (Exception ex) {
            System.err.println(ex);
            throw new IOException("unidentified problem with exporting data");
        }
    }

    static void verifyTableFormat(XSSFSheet sheet) throws IOException {
        Row row = sheet.getRow(0);
        Cell cell = row.getCell(1);

        if (cell.getCellType() == Cell.CELL_TYPE_STRING) {
            String formatLabel = cell.getStringCellValue();
            if (formatLabel.equals("HL-SEP-2016")) {
                return;
            }
        }

        throw new IOException("cannot verify table format");
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

    static void exportTotalCounts(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("tot cell num");

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix timeStamps = dataset.getTimeStamps();

        // header
        Row row = sheet.createRow(0);
        Cell cell = row.createCell(0);
        cell.setCellValue("time");

        int pos = 1;
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            for (int irep = 0; irep < dataset.numReplicates; irep++) {
                cell = row.createCell(pos++);
                cell.setCellValue(currName);
            }
        }

        // data
        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 1);
            cell = row.createCell(0);
            cell.setCellValue(currTime);

            pos = 1;
            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                for (int irep = 0; irep < dataset.numReplicates; irep++) {
                    cell = row.createCell(pos++);
                    if (dataset.missingStatreps[itpt][icnd][irep]==false) {
                        cell.setCellValue(dataset.totalNumbers[itpt][icnd][irep]);
                    }
                }
            }
        }
    }

    // export to Legacy data format (pure cell number information) - this can be used for an input format of Cyton Solver
    static void exportTotalCellPerCond(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("Input for CytonSolver");

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix missing = dataset.getMissing();
        SimpleMatrix timeStamps = dataset.getTimeStamps();
        SimpleMatrix counts = dataset.getCounts();

        // first row must be empty for legacy format
        int currRow = 0;
        Row row = sheet.createRow(currRow++);
        Cell cell = row.createCell(0);
        cell.setCellValue("This is an input worksheet of Cyton Solver. Please do NOT change tab name.");

        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            row = sheet.createRow(currRow++);
            cell = row.createCell(0);
            cell.setCellValue(condNames.get(icnd));

            int currCell = 1;
            for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
                double timeStamp = timeStamps.get(itpt);
                for (int irep = 0; irep < dataset.numReplicates; irep++) {

                    // check if there's at least one data point in the generation at that time point
                    int data_exist = 0;
                    for (int igen = 0; igen <= dataset.lastGeneration; igen++) {
                        int dataRow = igen + itpt * (dataset.lastGeneration + 1);
                        int dataCol = irep + icnd * dataset.numReplicates;
                        double isMissing = missing.get(dataRow, dataCol);

                        if (isMissing == 0) {
                            data_exist++;
                        }
                    }

                    // if there's at least one data point in the generation then print header
                    if (data_exist > 0) {
                        cell = row.createCell(currCell++);
                        cell.setCellValue(timeStamp);
                    }
                }
            }

            for (int igen = 0; igen <= dataset.lastGeneration; igen++) {
                row = sheet.createRow(currRow++);

                currCell = 1;
                for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
                    for (int irep = 0; irep < dataset.numReplicates; irep++) {
                        int dataRow = igen + itpt * (dataset.lastGeneration + 1);
                        int dataCol = irep + icnd * dataset.numReplicates;
                        double isMissing = missing.get(dataRow, dataCol);
                        double value = counts.get(dataRow, dataCol);

                        if (isMissing == 0) {
                            cell = row.createCell(currCell++);
                            cell.setCellValue(value);
                        }
                    }
                }
            }

            currRow++;
        }
    }

    static void exportCohortPerTime(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("cohort per TP");

        SimpleMatrix cohorts = dataset.getCohorts();
        SimpleMatrix missing = dataset.getMissing();
        SimpleMatrix timeStamps = dataset.getTimeStamps();
        List<String> conditionNames = dataset.getConditionNames();

        int currRow = 0;
        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            Row row = sheet.createRow(currRow++);
            Cell cell = row.createCell(0);
            cell.setCellValue(timeStamps.get(itpt));

            int currCell = 1;
            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                String condName = conditionNames.get(icnd);
                for (int irep = 0; irep < dataset.numReplicates; irep++) {
                    cell = row.createCell(currCell++);
                    cell.setCellValue(condName);
                }
            }

            for (int igen = 0; igen <= dataset.lastGeneration; igen++) {
                row = sheet.createRow(currRow++);

                currCell = 1;
                for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                    for (int irep = 0; irep < dataset.numReplicates; irep++) {
                        int dataRow = igen + itpt * (dataset.lastGeneration + 1);
                        int dataCol = irep + icnd * dataset.numReplicates;
                        double isMissing = missing.get(dataRow, dataCol);
                        double value = cohorts.get(dataRow, dataCol);

                        cell = row.createCell(currCell++);
                        if (isMissing == 0) {
                            cell.setCellValue(value);
                        }
                    }
                }
            }
        }
    }

    static void exportCohortPerCondition(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("cohort per cond");

        SimpleMatrix cohorts = dataset.getCohorts();
        SimpleMatrix missing = dataset.getMissing();
        SimpleMatrix timeStamps = dataset.getTimeStamps();
        List<String> conditionNames = dataset.getConditionNames();

        int currRow = 0;
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            Row row = sheet.createRow(currRow++);
            Cell cell = row.createCell(0);
            cell.setCellValue(conditionNames.get(icnd));

            int currCell = 1;
            for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
                double timeStamp = timeStamps.get(itpt);
                for (int irep = 0; irep < dataset.numReplicates; irep++) {
                    cell = row.createCell(currCell++);
                    cell.setCellValue(timeStamp);
                }
            }

            for (int igen = 0; igen <= dataset.lastGeneration; igen++) {
                row = sheet.createRow(currRow++);

                currCell = 1;
                for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
                    for (int irep = 0; irep < dataset.numReplicates; irep++) {
                        int dataRow = igen + itpt * (dataset.lastGeneration + 1);
                        int dataCol = irep + icnd * dataset.numReplicates;
                        double isMissing = missing.get(dataRow, dataCol);
                        double value = cohorts.get(dataRow, dataCol);

                        cell = row.createCell(currCell++);
                        if (isMissing == 0) {
                            cell.setCellValue(value);
                        }
                    }
                }
            }
        }
    }

    static void exporttotalcohortNum(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("Total Cohort Number");
        Row row;
        Cell cell;

        int rowOffset = 0;
        int colOffset = 0;

        //Title
        row = sheet.createRow(0);
        cell = row.createCell(0);
        cell.setCellValue("Total Cohort Number");

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix timeStamps = dataset.getTimeStamps();

        //Time
        row = sheet.createRow(1);
        cell = row.createCell(0);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            cell = row.createCell(icnd + 1);
            cell.setCellValue(currName);
        }

        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 2);
            cell = row.createCell(0);
            cell.setCellValue(currTime);
            //cohort numbers
            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                if (!dataset.missingStat[itpt][icnd]) {
                    cell = row.createCell(icnd + 1);
                    cell.setCellValue(dataset.totCohortnotnorm[itpt][icnd]);
                }
            }
        }

        // total cohort number excluding first t.p.
        rowOffset = dataset.numTimePoints + 3;
        colOffset = 0;

        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("total cohort number (excl. first t.p.)");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            cell = row.createCell(icnd + 1);
            cell.setCellValue(currName);
        }

        for (int itpt = 1; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 1 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {

                // NEW ONE (wihtout if statement)
                cell = row.createCell(icnd + 1 + colOffset);
                cell.setCellValue(dataset.totCohortnotnorm[itpt][icnd]);

                // OLD ONE
//                if (!dataset.missingEx0[itpt][icnd]) {
//                    cell = row.createCell(icnd + 1 + colOffset);
//                    cell.setCellValue(dataset.totCohortnotnorm[itpt][icnd]);
//                }
            }
        }

    }

    static void exportNormCohortNum(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("Norm Total Cohort Number");

        Row row;
        Cell cell;

        int rowOffset = 0;
        int colOffset = 0;

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix timeStamps = dataset.getTimeStamps();

        // total cohort number
        // title
        row = sheet.createRow(0);
        cell = row.createCell(0);
        cell.setCellValue("Normalised Total Cohort Number");

        // header
        row = sheet.createRow(1);
        cell = row.createCell(0);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            cell = row.createCell(icnd + 1);
            cell.setCellValue(currName);
        }
        //Set cohort numbers
        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 2);
            cell = row.createCell(0);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                if (!dataset.missingStat[itpt][icnd]) {
                    cell = row.createCell(icnd + 1);
                    cell.setCellValue(dataset.totCohort[itpt][icnd]);
                }
            }
        }

        // normalised total cohort number excluding first t.p.
        rowOffset = dataset.numTimePoints + 3;
        colOffset = 0;

        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("normalised total cohort number (excl. first t.p.)");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            cell = row.createCell(icnd + 1);
            cell.setCellValue(currName);
        }

        for (int itpt = 1; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 1 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {

                // NEW ONE (without if statement)
                cell = row.createCell(icnd + 1 + colOffset);
                cell.setCellValue(dataset.totExcTP0[itpt][icnd]);

                // OLD ONE
//                if (!dataset.missingEx0[itpt][icnd]) {
//                    cell = row.createCell(icnd + 1 + colOffset);
//                    cell.setCellValue(dataset.totExcTP0[itpt][icnd]);
//                }
            }
        }
    }

    static void exportTotalCohortNumreps(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("Total Cohort Num reps");

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix timeStamps = dataset.getTimeStamps();

        int rowOffset = 0;
        int colOffset = 0;

        Row row;
        Cell cell;

        // mean division number
        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("Total Cohort Number reps");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            for (int irep = 0; irep < dataset.numReplicates; irep++) {
                cell = row.createCell(icnd * dataset.numReplicates + irep + 1);
                cell.setCellValue(currName);
            }
        }

        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 2 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                for (int irep = 0; irep < dataset.numReplicates; irep++) {
                    if (!dataset.missingStatreps[itpt][icnd][irep]) {
                        cell = row.createCell(icnd * dataset.numReplicates + 1 + colOffset + irep);
                    if (!dataset.missingStatreps[itpt][icnd][irep] && !dataset.missingStatreps[0][icnd][irep]) {
                        cell.setCellValue(dataset.cohortSumreps[itpt][icnd][irep]);
                    }
                    }
                }
            }
        }
        // mean division number, excluding gen. 0
        rowOffset = dataset.numTimePoints + 3;
        colOffset = 0;

        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("Total Cohort Number (excl. first t.p.)");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            for (int irep = 0; irep < dataset.numReplicates; irep++) {
                cell = row.createCell(icnd * dataset.numReplicates + irep + 1);
                cell.setCellValue(currName);
            }
        }

        for (int itpt = 1; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 1 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                for (int irep = 0; irep < dataset.numReplicates; irep++) {

                    // NEW ONE (without missingEx0reps flag)
                    cell = row.createCell(icnd * dataset.numReplicates + 1 + colOffset + irep);
                    if (!dataset.missingStatreps[itpt][icnd][irep] && !dataset.missingStatreps[1][icnd][irep]) {
                        cell.setCellValue(dataset.cohortSumreps[itpt][icnd][irep]);
                    }

                    // OLD ONE
//                    if (!dataset.missingEx0reps[itpt][icnd][irep]) {
//                        cell = row.createCell(icnd * 3 + 1 + colOffset + irep);
//                        if (!dataset.missingStatreps[itpt][icnd][irep] && !dataset.missingStatreps[1][icnd][irep]) {
//                            cell.setCellValue(dataset.cohortSumreps[itpt][icnd][irep]);
//                        }
//                    }
                }
            }
        }

    }

    static void exportNormCohortNumreps(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("Norm Total Cohort Num reps");

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix timeStamps = dataset.getTimeStamps();

        int rowOffset = 0;
        int colOffset = 0;

        Row row;
        Cell cell;

        // mean division number
        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("Norm Total Cohort Number");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            for (int irep = 0; irep < dataset.numReplicates; irep++) {
                cell = row.createCell(icnd * dataset.numReplicates + irep + 1);
                cell.setCellValue(currName);
            }
        }

        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 2 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                for (int irep = 0; irep < dataset.numReplicates; irep++) {
                    cell = row.createCell(icnd * dataset.numReplicates + 1 + colOffset + irep);
                    if (dataset.missingStatreps[itpt][icnd][irep] ==false && dataset.missingStatreps[0][icnd][irep]==false) {
                        cell.setCellValue(dataset.totCohorteach[itpt][icnd][irep]);
                    }
                }
            }
        }
        // mean division number, excluding gen. 0
        rowOffset = dataset.numTimePoints + 3;
        colOffset = 0;

        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("Norm Total Cohort Number (excl. first t.p.)");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            for (int irep = 0; irep < dataset.numReplicates; irep++) {
                cell = row.createCell(icnd * dataset.numReplicates + irep + 1);
                cell.setCellValue(currName);
            }
        }

        for (int itpt = 1; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 1 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                for (int irep = 0; irep < dataset.numReplicates; irep++) {
                    
                    cell = row.createCell(icnd * dataset.numReplicates + 1 + colOffset + irep);
                    if (dataset.missingStatreps[itpt][icnd][irep] ==false && dataset.missingStatreps[1][icnd][irep]==false) {
                        cell.setCellValue(dataset.totexcTP0each[itpt][icnd][irep]);
                    }

                }
            }
        }

    }

    static void exportmdnvscohort(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("Mean Division Number vs Total Cohort Number");

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix timeStamps = dataset.getTimeStamps();

        int rowOffset = 0;
        int colOffset = 0;

        Row row;
        Cell cell;

        // mean division number
        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("Mean Division number vs Total cohort Number");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            for (int i = 0; i < (dataset.numReplicates + 1); i++) {
                if (i == 0) {
                    cell = row.createCell(icnd * (dataset.numReplicates + 1) + i);
                    cell.setCellValue("MDN");
                } else {
                    String currName = condNames.get(icnd);
                    cell = row.createCell(icnd * (dataset.numReplicates + 1) + i);
                    cell.setCellValue(currName);
                }
            }
        }

        int irep=0;

        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {

            row = sheet.createRow(itpt + 2 + rowOffset);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                irep=0;
                for (int i = 0; i < dataset.numReplicates + 1; i++) {
                    if (i == 0) {
                        if (!dataset.missingStat[itpt][icnd]) {
                            cell = row.createCell(icnd * (dataset.numReplicates + 1) + colOffset);
                            cell.setCellValue(dataset.meanDivNum[itpt][icnd]);
                        }
                    } else {
                        cell = row.createCell(icnd * (dataset.numReplicates + 1) + i + colOffset);
                        if (dataset.missingStatreps[itpt][icnd][irep]==false)
                    {
                        cell.setCellValue(dataset.cohortSumreps[itpt][icnd][irep]);
                    }
                        irep=irep+1;
                    }

                }
            }
        }
        // mean division number, excluding gen. 0
        rowOffset = dataset.numTimePoints + 3;
        colOffset = 0;

        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("mean division number (gen 0 excluded) vs total cohort number");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);

        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            for (int i = 0; i < (dataset.numReplicates + 1); i++) {
                if (i == 0) {
                    cell = row.createCell(icnd * (dataset.numReplicates + 1) + i);
                    cell.setCellValue("MDN");
                } else {
                    String currName = condNames.get(icnd);
                    cell = row.createCell(icnd * (dataset.numReplicates + 1) + i);
                    cell.setCellValue(currName);
                }
            }
        }


        for (int itpt = 1; itpt < dataset.numTimePoints; itpt++) {

            row = sheet.createRow(itpt + 2 + rowOffset-1);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                irep=0;
                for (int i = 0; i < dataset.numReplicates + 1; i++) {
                    if (i == 0) {
                        if (!dataset.missingStat[itpt][icnd]) {
                            cell = row.createCell(icnd * (dataset.numReplicates + 1) + colOffset);
//                            cell.setCellValue(dataset.meanDivNum[itpt][icnd]);
                            cell.setCellValue(dataset.mdnExcGen0[itpt][icnd]);
                        }
                    } else {
                        cell = row.createCell(icnd * (dataset.numReplicates + 1) + i + colOffset);
                        if (dataset.missingStatreps[itpt][icnd][irep]==false)
                    {
                        cell.setCellValue(dataset.cohortSumreps[itpt][icnd][irep]);
                    }
                        irep=irep+1;

                    }

                }
            }
        }
    }

    static void exportMeanDivNum(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("Mean Division Number");

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix timeStamps = dataset.getTimeStamps();

        int rowOffset = 0;
        int colOffset = 0;

        Row row;
        Cell cell;

        // mean division number
        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("mean division number");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            cell = row.createCell(icnd + 1);
            cell.setCellValue(currName);
        }

        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 2 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                if (!dataset.missingStat[itpt][icnd]) {
                    cell = row.createCell(icnd + 1 + colOffset);
                    cell.setCellValue(dataset.meanDivNum[itpt][icnd]);
                }
            }
        }

        // mean division number, excluding gen. 0
        rowOffset = dataset.numTimePoints + 3;
        colOffset = 0;

        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("mean division number (gen 0 excluded)");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            cell = row.createCell(icnd + 1);
            cell.setCellValue(currName);
        }

        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 2 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                if (!dataset.missingEx0[itpt][icnd]) {
                    cell = row.createCell(icnd + 1 + colOffset);
                    cell.setCellValue(dataset.mdnExcGen0[itpt][icnd]);
                }
            }
        }
    }

    static void exportMeanDivNumreps(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("Mean Division Number reps");

        List<String> condNames = dataset.getConditionNames();
        SimpleMatrix timeStamps = dataset.getTimeStamps();

        int rowOffset = 0;
        int colOffset = 0;

        Row row;
        Cell cell;

        // mean division number
        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("mean division number");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            for (int irep = 0; irep < dataset.numReplicates; irep++) {
                cell = row.createCell(icnd * dataset.numReplicates + irep + 1);
                cell.setCellValue(currName);
            }
        }

        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 2 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                for (int irep = 0; irep < dataset.numReplicates; irep++) {
                    if (!dataset.missingStatreps[itpt][icnd][irep]) {
                        cell = row.createCell(icnd * dataset.numReplicates + 1 + colOffset + irep);
                        if (!dataset.missingStatreps[itpt][icnd][irep]) {
                            cell.setCellValue(dataset.meanDiveach[itpt][icnd][irep]);
                        }
                    }
                }
            }
        }
        // mean division number, excluding gen. 0
        rowOffset = dataset.numTimePoints + 3;
        colOffset = 0;

        // title
        row = sheet.createRow(0 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("mean division number (gen 0 excluded)");

        // header
        row = sheet.createRow(1 + rowOffset);
        cell = row.createCell(0 + colOffset);
        cell.setCellValue("time");
        for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
            String currName = condNames.get(icnd);
            for (int irep = 0; irep < dataset.numReplicates; irep++) {
                cell = row.createCell(icnd * dataset.numReplicates + irep + 1);
                cell.setCellValue(currName);
            }
        }

        for (int itpt = 0; itpt < dataset.numTimePoints; itpt++) {
            double currTime = timeStamps.get(itpt);

            row = sheet.createRow(itpt + 2 + rowOffset);
            cell = row.createCell(0 + colOffset);
            cell.setCellValue(currTime);

            for (int icnd = 0; icnd < dataset.numConditions; icnd++) {
                for (int irep = 0; irep < dataset.numReplicates; irep++) {
                    if (!dataset.missingEx0reps[itpt][icnd][irep]) {
                        cell = row.createCell(icnd * dataset.numReplicates + 1 + colOffset + irep);
                        if (!dataset.missingStatreps[itpt][icnd][irep]) {
                            cell.setCellValue(dataset.mdnExcGen0reps[itpt][icnd][irep]);
                        }
                    }
                }
            }
        }

    }

    static void exportTimeStamp(Dataset dataset, XSSFWorkbook workbook) throws IOException {
        XSSFSheet sheet = workbook.createSheet("version");

        //Setting Time and Date, Time stamping
        CellStyle cellStyle = workbook.createCellStyle();
        CreationHelper createHelper = workbook.getCreationHelper();
        short dateFormat = createHelper.createDataFormat().getFormat("h:mm dd-MM-yyyy");
        cellStyle.setDataFormat(dateFormat);
        Row row = sheet.createRow(0);
        Cell cell = row.createCell(1);
        cell.setCellValue(Calendar.getInstance());
        cell.setCellStyle(cellStyle);
        row.createCell(0).setCellValue("The Spreadsheet was created at:");
        sheet.autoSizeColumn(0);
    }
}
