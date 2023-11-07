package edu.wehi.cohort.io;


import java.awt.Dimension;
import java.awt.Toolkit;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import de.erichseifert.gral.graphics.Insets2D;
import de.erichseifert.gral.graphics.layout.TableLayout;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.graphics.DrawableContainer;
import de.erichseifert.gral.io.plots.DrawableWriter;
import de.erichseifert.gral.io.plots.DrawableWriterFactory;

public class ExportPlot {

    public static void toPNG(DrawableContainer plot, File dest, String fileName) {
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();

        double xReduce = 0;
        double X = 0;
        double Y = 0;

        /* check if container has more than 1 plot panel */
        if (plot.getDrawables().toArray().length > 2) {
            xReduce = screenSize.width*0.06;
            X = screenSize.width - xReduce;
            Y = screenSize.height*X/screenSize.width;
        }
        else {
            xReduce = screenSize.width*0.3;
            X = screenSize.width - xReduce;
            Y = screenSize.height*X/screenSize.width;
        }

        DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/png");

        try {
            String destination = dest.getCanonicalPath() + "/";
            writer.write(plot, new FileOutputStream(destination + fileName +".png"), X, Y);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void toSVG(XYPlot[] plotObjs, String figNum, File dest) {

        for (XYPlot obj : plotObjs) {

            DrawableContainer WRAPPER = new DrawableContainer(new TableLayout(1));
            obj.getPlotArea().setBorderStroke(null);
            WRAPPER.add(obj);
            Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();

            double perc = 0.4;
            double xReduce = screenSize.width * perc;
            double X = screenSize.width - xReduce;
            double Y = screenSize.height * X / screenSize.width;

            WRAPPER.setInsets(new Insets2D.Double(0.0, 50.0, 0.0, 50.0));

            DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/svg+xml");

            // Strange bug in printing PDF/EPS format
//            DrawableWriter writer = DrawableWriterFactory.getInstance().get("application/pdf");
//            DrawableWriter writer = DrawableWriterFactory.getInstance().get("application/postscript");

            try {
                String destination = dest.getCanonicalPath() + "/";
                writer.write(WRAPPER, new FileOutputStream(destination + figNum + obj.getTitle().getText() + ".svg"), 0, 0, X, Y);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
